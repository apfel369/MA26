# The script was used to filter the multisample VCF, which serves as the basis for the GT workflow of unsupervised learning methods. 

#!/usr/bin/env bash
# =============================================================================
# Multi-Sample VCF Genotype-Level Filter – Step 01
# Filterkriterien:
#   FORMAT/DP  >= 5 und <= 2x Median(DP) [PER SAMPLE]
#   MQ         >= 40
#   FORMAT/GQ  >= 20
#   QD (custom) = QUAL / SUM(FORMAT/DP aller non-homref/non-missing) >= 2
#   MQBZ        > -4   (nur wenn Feld vorhanden)
#   RPBZ        >= -5 und <= 5  (nur wenn Feld vorhanden)
#
# Strategie: Genotypen, die Grenzen nicht bestehen -> ./. (missing) gesetzt
#            Site-level INFO-Filter (MQ, MQBZ, RPBZ, QD) -> soft-filter
# =============================================================================

set -euo pipefail

# -- Plugin-Pfad setzen --------------------------------------------------------
export BCFTOOLS_PLUGINS="/PFAD/bcftools"

if [[ ! -f "${BCFTOOLS_PLUGINS}/setGT.so" ]]; then
    echo "[ERROR] setGT.so nicht gefunden in: ${BCFTOOLS_PLUGINS}"
    exit 1
fi
echo "[INFO] BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS}"

# -- Pfade ---------------------------------------------------------------------
INPUT_DIR="/PFAD/"
INPUT_VCF="${INPUT_DIR}/multisample-snvs-only_expanded.sorted.vcf.gz"
OUT_DIR="${INPUT_DIR}/filtered-genotype-level_step-01"
OUT_VCF="${OUT_DIR}/multisample-genotype-filtered_step-01_korrigiert_expanded_complete.vcf.gz"

THREADS=18

mkdir -p "${OUT_DIR}"


# -- Schritt 1: DP-Grenzen aus TXT einlesen -----------------------------------
DP_LIMITS_TXT="${INPUT_DIR}/dp_filtergrenzen.txt"

if [[ ! -f "${DP_LIMITS_TXT}" ]]; then
    echo "[ERROR] DP-Grenzen-TXT nicht gefunden: ${DP_LIMITS_TXT}"
    exit 1
fi

echo "[INFO] Lese DP-Grenzen aus: ${DP_LIMITS_TXT}"

mapfile -t SAMPLES < <(bcftools query -l "${INPUT_VCF}")
declare -A SAMPLE_MAX_DP

while read -r SAMPLE MIN_DP MAX_DP MEDIAN_DP; do
    # Kommentarzeilen überspringen
    [[ "${SAMPLE}" == \#* ]] && continue
    SAMPLE_MAX_DP["${SAMPLE}"]="${MAX_DP}"
    echo "[INFO]   ${SAMPLE}: Median DP = ${MEDIAN_DP} -> MAX_DP = ${MAX_DP}"
done < "${DP_LIMITS_TXT}"

# Sicherheitscheck: Alle Samples der VCF müssen in der TXT stehen
for S in "${SAMPLES[@]}"; do
    if [[ -z "${SAMPLE_MAX_DP[${S}]+x}" ]]; then
        echo "[ERROR] Sample '${S}' fehlt in der DP-Grenzen-TXT!"
        exit 1
    fi
done

# -- Schritt 2: Custom-QD annotieren ------------------------------------------
echo "[INFO] Berechne und annotiere Custom-QD ..."

PYTHON_QD_SCRIPT=$(mktemp /tmp/calc_qd_XXXXXX.py)

cat > "${PYTHON_QD_SCRIPT}" << 'PYEOF'
import sys

for line in sys.stdin:
    line = line.rstrip('\n')
    if line.startswith('##INFO=<ID=CUSTOM_QD'): continue
    if line.startswith('#CHROM'):
        print('##INFO=<ID=CUSTOM_QD,Number=1,Type=Float,Description="Custom QD: QUAL/SUM(FMT/DP non-homref non-missing samples)">')
        print(line); continue
    if line.startswith('#'): print(line); continue
    fields = line.split('\t')
    if len(fields) < 10: print(line); continue
    try: qual = float(fields[5]) if fields[5] not in ('.','') else None
    except ValueError: qual = None
    fmt = fields[8].split(':')
    dp_idx = fmt.index('DP') if 'DP' in fmt else None
    gt_idx = fmt.index('GT') if 'GT' in fmt else None
    dp_sum = 0; valid = False
    if dp_idx is not None and gt_idx is not None and qual is not None:
        for sample in fields[9:]:
            sf = sample.split(':')
            gt = sf[gt_idx] if gt_idx < len(sf) else '.'
            gt_alleles = gt.replace('|','/').split('/')
            if '.' in gt_alleles: continue
            if all(a=='0' for a in gt_alleles): continue
            try:
                dp_val = int(sf[dp_idx]) if dp_idx < len(sf) and sf[dp_idx]!='.' else 0
                dp_sum += dp_val; valid = True
            except ValueError: pass
    if valid and dp_sum > 0:
        cqd = round(qual/dp_sum, 4)
        info = fields[7]
        if info=='.': fields[7]='CUSTOM_QD='+str(cqd)
        elif 'CUSTOM_QD=' in info:
            parts=[p for p in info.split(';') if not p.startswith('CUSTOM_QD=')]
            parts.append('CUSTOM_QD='+str(cqd)); fields[7]=';'.join(parts)
        else: fields[7]=info+';CUSTOM_QD='+str(cqd)
    print('\t'.join(fields))
PYEOF

bcftools view --threads "${THREADS}" "${INPUT_VCF}" \
    | python3 "${PYTHON_QD_SCRIPT}" \
    | bgzip -@ "${THREADS}" \
    > "${OUT_DIR}/tmp_qd_annotated.vcf.gz"

bcftools index --threads "${THREADS}" --tbi "${OUT_DIR}/tmp_qd_annotated.vcf.gz"
rm -f "${PYTHON_QD_SCRIPT}"

# -- Schritt 3: Genotypen auf missing setzen (DP per-sample + GQ) -------------
# Da bcftools +setGT keinen per-Sample-MAX_DP unterstuetzt, wird ein
# Python-Skript genutzt, das die sample-spezifischen Grenzen direkt anwendet.
echo "[INFO] Setze Genotypen auf missing (per-sample DP-Grenzen + GQ) ..."

PYTHON_GT_SCRIPT=$(mktemp /tmp/mask_gt_XXXXXX.py)

# Sample-MAX_DP als Leerzeichen-getrennte Liste uebergeben:
# Format: SAMPLE1=MAX1 SAMPLE2=MAX2 ...
MAX_DP_ARG=""
for SAMPLE in "${SAMPLES[@]}"; do
    MAX_DP_ARG="${MAX_DP_ARG} ${SAMPLE}=${SAMPLE_MAX_DP[${SAMPLE}]}"
done

cat > "${PYTHON_GT_SCRIPT}" << PYEOF
import sys

# Per-Sample MAX_DP aus Kommandozeilenargumenten einlesen
# Aufruf: python3 skript.py SAMPLE1=MAX1 SAMPLE2=MAX2 ...
sample_max_dp = {}
for arg in sys.argv[1:]:
    name, val = arg.split('=')
    sample_max_dp[name] = int(val)

sample_order = []  # wird beim ersten #CHROM-Header gesetzt

for line in sys.stdin:
    line = line.rstrip('\n')

    if line.startswith('#CHROM'):
        fields = line.split('\t')
        sample_order = fields[9:]
        print(line)
        continue
    if line.startswith('#'):
        print(line)
        continue

    fields = line.split('\t')
    if len(fields) < 10:
        print(line)
        continue

    fmt = fields[8].split(':')
    gt_idx  = fmt.index('GT') if 'GT' in fmt else None
    dp_idx  = fmt.index('DP') if 'DP' in fmt else None
    gq_idx  = fmt.index('GQ') if 'GQ' in fmt else None

    new_samples = []
    for i, sample_field in enumerate(fields[9:]):
        sample_name = sample_order[i] if i < len(sample_order) else None
        max_dp = sample_max_dp.get(sample_name, 999999)

        sf = sample_field.split(':')

        # Genotyp auslesen
        if gt_idx is None or gt_idx >= len(sf):
            new_samples.append(sample_field)
            continue

        gt = sf[gt_idx]
        # bereits missing -> unberuehrt lassen
        if '.' in gt.replace('|','/').split('/'):
            new_samples.append(sample_field)
            continue

        # DP pruefen
        dp_fail = False
        if dp_idx is not None and dp_idx < len(sf) and sf[dp_idx] != '.':
            try:
                dp = int(sf[dp_idx])
                if dp < 5 or dp > max_dp:
                    dp_fail = True
            except ValueError:
                dp_fail = True
        else:
            dp_fail = True  # kein DP vorhanden -> als fehlend behandeln

        # GQ pruefen
        gq_fail = False
        if gq_idx is not None and gq_idx < len(sf) and sf[gq_idx] != '.':
            try:
                gq = int(sf[gq_idx])
                if gq < 20:
                    gq_fail = True
            except ValueError:
                gq_fail = True
        else:
            gq_fail = True  # kein GQ vorhanden -> als fehlend behandeln

        if dp_fail or gq_fail:
            # GT auf missing setzen, Rest des FORMAT-Feldes beibehalten
            sf[gt_idx] = './.'
            new_samples.append(':'.join(sf))
        else:
            new_samples.append(sample_field)

    fields[9:] = new_samples
    print('\t'.join(fields))
PYEOF

bcftools view --threads "${THREADS}" "${OUT_DIR}/tmp_qd_annotated.vcf.gz" \
    | python3 "${PYTHON_GT_SCRIPT}" ${MAX_DP_ARG} \
    | bgzip -@ "${THREADS}" \
    > "${OUT_DIR}/tmp_gt_masked.vcf.gz"

bcftools index --threads "${THREADS}" --tbi "${OUT_DIR}/tmp_gt_masked.vcf.gz"
rm -f "${PYTHON_GT_SCRIPT}"

# -- Schritt 4: Site-Level Filter (MQ, MQBZ, RPBZ, CUSTOM_QD) -----------------
echo "[INFO] Wende Site-Level Filter an ..."

bcftools filter \
    --threads "${THREADS}" \
    --soft-filter 'GenotypeLevelFilter_Step01' \
    --mode '+' \
    --exclude 'MQ < 40 || (INFO/MQBZ != "." && INFO/MQBZ <= -4) || (INFO/RPBZ != "." && (INFO/RPBZ < -5 || INFO/RPBZ > 5)) || (INFO/CUSTOM_QD != "." && INFO/CUSTOM_QD < 2)' \
    "${OUT_DIR}/tmp_gt_masked.vcf.gz" \
    -Oz -o "${OUT_VCF}"

bcftools index --threads "${THREADS}" --tbi "${OUT_VCF}"

# -- Aufraeumen ----------------------------------------------------------------
rm -f "${OUT_DIR}"/tmp_qd_annotated.vcf.gz* \
      "${OUT_DIR}"/tmp_gt_masked.vcf.gz*

echo "[DONE] Ausgabe: ${OUT_VCF}"

# -- Statistik -----------------------------------------------------------------
echo "[INFO] Filter-Statistik:"
bcftools stats --threads "${THREADS}" "${OUT_VCF}" \
    | grep -E '^SN|^FILTER'

# -- nur PASS-Sites behalten -----------------------------------------------------------------
bcftools view -f PASS \
multisample-genotype-filtered_step-01_korrigiert_expanded_complete.vcf.gz \
-Oz -o multisample-genotype-filtered_step-01_korrigiert_expanded_complete_pass.vcf.gz
bcftools index multisample-genotype-filtered_step-01_korrigiert_expanded_complete_pass.vcf.gz
