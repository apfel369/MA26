# -*- coding: utf-8 -*-
import csv
from cyvcf2 import VCF

# --- DATEIPFADE ---
VCF_FILE = "/PFAD/v2-HBeeID_272_markers_POS_MATCH_noInDels.vcf"
SNP_ORDER_FILE = "/PFAD/excel_snps.txt"
OUT_CSV = "/PFAD/v2-HBeeID_Input.csv"

# ==========================================
# 1. Gewünschte SNP-Reihenfolge einlesen
# ==========================================
snp_order = []
with open(SNP_ORDER_FILE, 'r') as f:
    for line in f:
        snp_id = line.strip()
        if snp_id:  # Leere Zeilen überspringen
            snp_order.append(snp_id)

print(f"Lese {len(snp_order)} SNP-IDs aus der Textdatei ein...")

# ==========================================
# 2. VCF auslesen und Genotypen übersetzen
# ==========================================
vcf = VCF(VCF_FILE)
samples = vcf.samples

# Ein Dictionary vorbereiten: sample_name -> {snp_id: genotyp_code}
genotype_data = {sample: {} for sample in samples}

for variant in vcf:
    snp_id = variant.ID
    
    # Wenn die Variante keine ID hat, überspringen wir sie
    if snp_id is None or snp_id == '.':
        continue

    for i, sample in enumerate(samples):
        # cyvcf2 liefert Genotypen als Array: [Allel1, Allel2, Phasing-Flag]
        # -1 bedeutet fehlend (z. B. ./. durch deine vorherige Filterung)
        a1, a2, _ = variant.genotypes[i]

        if a1 == -1 or a2 == -1:
            g_code = "NA"
        elif a1 == 0 and a2 == 0:     # 0/0 -> homozygot Referenz
            g_code = 0
        elif (a1 == 0 and a2 == 1) or (a1 == 1 and a2 == 0): # 0/1 oder 1/0 -> heterozygot
            g_code = 1
        elif a1 == 1 and a2 == 1:     # 1/1 -> homozygot Alternativ
            g_code = 2
        else:
            # Fallback für unerwartete Genotypen (z. B. multiallelisch 1/2)
            g_code = "NA"

        genotype_data[sample][snp_id] = g_code

vcf.close()

# ==========================================
# 3. CSV-Datei schreiben (im HBeeID-Format)
# ==========================================
with open(OUT_CSV, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)

    # Kopfzeile: Erste Zelle ist "Sample_ID", danach exakt die SNPs aus der txt
    header = ["Sample_ID"] + snp_order
    writer.writerow(header)

    # Datenzeilen: Pro Biene eine Zeile
    for sample in samples:
        row = [sample]
        for snp_id in snp_order:
            # Wir rufen den Genotyp ab. Falls der SNP in der VCF-Datei komplett
            # fehlt (z.B. gelöscht beim Liftover), wird stattdessen "NA" eingetragen.
            genotype = genotype_data[sample].get(snp_id, "NA")
            row.append(genotype)
        
        writer.writerow(row)

print(f"Fertig! Die Datei wurde gespeichert unter:\n{OUT_CSV}")