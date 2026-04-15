# -*- coding: iso-8859-1 -*-
from cyvcf2 import VCF, Writer

IN_VCF  = "/PFAD ZUR DATEI/HBeeID_272_markers_Amel45_multisample_biallelic_noIndel.vcf"
OUT_VCF = "/PFAD ZUR DATEI/HBeeID_272_markers_FILTERED_Amel45.vcf"

# --- PARAMETER & SCHWELLENWERTE ---
MQ_MIN         = 35
QD_MIN         = 2.0
GQ_MIN         = 20
DP_MIN_DEFAULT = 15
DP_MIN_B1_1    = 15
MISSING        = -2147483648

# Spezifische Max-DP-Werte. Samples, die nicht drin stehen oder 'None' haben, bekommen kein Limit.
MAX_DP = {
    "B1-1_sorted":     None, 
    "B13018-3_sorted": None,
    "B13018-4_sorted": None,
    "B20-6_sorted":    None,
    "B26-3_sorted":    None,
    "B34-4_sorted":    None,
    "B40-6_sorted":    None,
    "B57-1_sorted":    None,
    "B6-1_sorted":     None,
    "B6-3_sorted":     None,
}

# --- HILFSFUNKTION ---
def safe_format(variant, tag):
    """
    Sicheres Abrufen von FORMAT-Tags. F ngt cyvcf2 KeyError ab, 
    falls der Tag im Header (wie z.B. DP) fehlt.
    """
    try:
        return variant.format(tag)
    except KeyError:
        return None

# --- HAUPTSKRIPT ---
vcf     = VCF(IN_VCF)
samples = vcf.samples
print("Samples: " + str(samples), flush=True)

w = Writer(OUT_VCF, vcf)
n_total = n_global_fail = n_written = 0

for variant in vcf:
    n_total += 1
    fail_global = False

    # 1. Global: MQ >= 35 (Wie von dir vorgegeben: ">=" 35)
    mq = variant.INFO.get("MQ")
    if mq is not None and mq < MQ_MIN:
        fail_global = True

    # 2. Global: MQBZ > -4 (Wenn NICHT strikt gr  er als -4, dann filtern)
    if not fail_global:
        mqbz = variant.INFO.get("MQBZ")
        if mqbz is not None and not (mqbz > -4.0):
            fail_global = True

    # 3. Global: -5 < RPBZ < 5 (Wenn NICHT strikt zwischen -5 und 5, dann filtern)
    if not fail_global:
        rpbz = variant.INFO.get("RPBZ")
        if rpbz is not None and not (-5.0 < rpbz < 5.0):
            fail_global = True

    # 4. Global: QD >= 2.0 (Berechnet aus globaler QUAL / globaler DP)
    if not fail_global:
        qual    = variant.QUAL
        info_dp = variant.INFO.get("DP")
        if qual is not None and info_dp is not None and info_dp > 0:
            if (qual / info_dp) < QD_MIN:
                fail_global = True

    gt = variant.genotypes
    
    # Wenn ein globaler Filter greift: Setze alle Proben dieser Zeile auf "fehlend" (./.)
    if fail_global:
        n_global_fail += 1
        for i in range(len(gt)):
            gt[i][0], gt[i][1], gt[i][2] = -1, -1, False
        variant.genotypes = gt
        w.write_record(variant)
        n_written += 1
        continue

    # ==========================================
    # LOKALE FILTER (Anwendung pro Biene)
    # ==========================================
    pl_data = safe_format(variant, "PL")
    gq_data = safe_format(variant, "GQ")
    dp_data = safe_format(variant, "DP")

    for i, sname in enumerate(samples):
        gq = None
        if gq_data is not None:
            v = int(gq_data[i][0])
            if v != MISSING:
                gq = v
                
        # Wenn kein GQ, aber PL vorhanden -> GQ berechnen
        if gq is None and pl_data is not None:
            pls = [int(p) for p in pl_data[i] if int(p) != MISSING]
            if len(pls) >= 2:
                sp = sorted(pls)
                gq = min(sp[1] - sp[0] if sp[0] == 0 else sp[0], 99)

        dp = None
        if dp_data is not None:
            v = int(dp_data[i][0])
            if v != MISSING:
                dp = v

        # Individuelle DP Min-Zuweisung pr fen (B1-1 vs Default)
        dp_min = DP_MIN_B1_1 if sname == "B1-1_sorted" else DP_MIN_DEFAULT
        dp_max = MAX_DP.get(sname)

        # Pr fen, ob lokaler Filter triggert
        set_missing = (
            (gq is not None and gq < GQ_MIN) or
            (dp is not None and dp < dp_min) or
            (dp_max is not None and dp is not None and dp > dp_max)
        )
        
        # Lokalen Genotyp auf ./. setzen (-1, -1, False), falls Filter greift
        if set_missing:
            gt[i][0], gt[i][1], gt[i][2] = -1, -1, False

    variant.genotypes = gt
    w.write_record(variant)
    n_written += 1

w.close()
vcf.close()

print("Gesamt:          " + str(n_total))
print("Global fail:     " + str(n_global_fail))
print("Geschrieben:     " + str(n_written))
print("Output: " + OUT_VCF)
