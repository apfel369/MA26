#!/bin/bash

# =========================================================
# Shell-Wrapper: InDels entfernen (auf NA setzen)
# =========================================================

echo "Pruefe auf InDels und setze Eintraege auf fehlend (./.)..."

python - << 'EOF'
from cyvcf2 import VCF, Writer

# --- DATEIPFADE ---
IN_VCF  = "/PFAD/v2-HBeeID_272_markers_POS_MATCH.vcf"
OUT_VCF = "/PPFAD/v2-HBeeID_272_markers_POS_MATCH_noInDels.vcf"

vcf = VCF(IN_VCF)
w = Writer(OUT_VCF, vcf)

n_total = 0
n_indels = 0

for variant in vcf:
    n_total += 1
    
    # 1. Pruefen, ob ein InDel vorliegt
    is_indel = variant.is_indel
    
    if not is_indel:
        # Manuelle Pruefung der String-Laengen zur Sicherheit
        if len(variant.REF) != 1:
            is_indel = True
        for alt in variant.ALT:
            if alt and alt != '.' and len(alt) != 1:
                is_indel = True

    # 2. Wenn InDel: Alle Genotypen auf -1 (NA) setzen
    if is_indel:
        n_indels += 1
        gt = variant.genotypes
        for i in range(len(gt)):
            # -1 bedeutet in cyvcf2 "fehlend" (./.)
            gt[i][0], gt[i][1], gt[i][2] = -1, -1, False
        variant.genotypes = gt
        
    w.write_record(variant)

w.close()
vcf.close()

print("Fertig!")
print("Gesamt verarbeitet: " + str(n_total))
print("Davon InDels auf NA (./.) gesetzt: " + str(n_indels))
print("Output gespeichert unter: " + OUT_VCF)

EOF

echo "Skript erfolgreich beendet."