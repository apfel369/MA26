# -*- coding: utf-8 -*-
import csv
from cyvcf2 import VCF

# --- DATEIPFADE ---
VCF_FILE = "/PFAD ZUR DATEI/HBeeID_272_markers_FILTERED_Amel45.vcf"
SNP_ORDER_FILE = "/PFAD ZUR DATEI/excel_snps.txt"
OUT_CSV = "/PFAD ZUR DATEI/HBeeID_Input-amel45.csv"

# ==========================================
# 1. SNP-Reihenfolge einlesen + POS-Index aufbauen
# ==========================================
snp_order = []
# pos_to_snpid: { position_als_int -> snp_id }
pos_to_snpid = {}

with open(SNP_ORDER_FILE, 'r') as f:
    for line in f:
        snp_id = line.strip()
        if not snp_id:
            continue
        snp_order.append(snp_id)
        # Letzte Zahl nach dem letzten Unterstrich = Position
        pos = int(snp_id.rsplit('_', 1)[-1])
        pos_to_snpid[pos] = snp_id

print(f"Eingelesene SNPs: {len(snp_order)}")
print(f"Positionen im Index: {len(pos_to_snpid)}")

# ==========================================
# 2. VCF auslesen und Genotypen  bersetzen
# ==========================================
vcf = VCF(VCF_FILE)
samples = vcf.samples
genotype_data = {sample: {} for sample in samples}

matched = 0
skipped = 0

for variant in vcf:
    pos = variant.POS  # 1-basierte Position in cyvcf2

    if pos not in pos_to_snpid:
        skipped += 1
        continue

    snp_id = pos_to_snpid[pos]
    matched += 1

    for i, sample in enumerate(samples):
        a1, a2, _ = variant.genotypes[i]

        if a1 == -1 or a2 == -1:
            g_code = "NA"
        elif a1 == 0 and a2 == 0:
            g_code = 0
        elif (a1 == 0 and a2 == 1) or (a1 == 1 and a2 == 0):
            g_code = 1
        elif a1 == 1 and a2 == 1:
            g_code = 2
        else:
            g_code = "NA"  # z.B. multiallelisch

        genotype_data[sample][snp_id] = g_code

vcf.close()

print(f"Gematchte Varianten: {matched}")
print(f"Nicht gematchte Varianten (Position nicht in txt): {skipped}")

# ==========================================
# 3. CSV-Datei schreiben (im HBeeID-Format)
# ==========================================
with open(OUT_CSV, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    header = ["Sample_ID"] + snp_order
    writer.writerow(header)

    for sample in samples:
        row = [sample]
        for snp_id in snp_order:
            row.append(genotype_data[sample].get(snp_id, "NA"))
        writer.writerow(row)

print(f"\nFertig! Gespeichert unter:\n{OUT_CSV}")