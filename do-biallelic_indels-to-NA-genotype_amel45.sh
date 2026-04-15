#!/bin/bash

# ==========================================
# Post-Processing: Split multiallelics + InDels maskieren
# ==========================================

REF_FASTA="/PFAD ZUR DATEI/Amel_4.5/GCF_000002195.4_Amel_4.5_genomic.fna"
IN_VCF="/PFAD ZUR DATEI/Amel_4.5/HBeeID-multi-sample-vcf-272markers_Amel45/HBeeID_272_markers_Amel45_multisample.vcf"
OUT_VCF="/PFAD ZUR DATEI/Amel_4.5/HBeeID-multi-sample-vcf-272markers_Amel45/HBeeID_272_markers_Amel45_multisample_biallelic_noIndel.vcf.gz"
THREADS=8

echo "Starte Post-Processing..."

# Schritt 1: Multiallelische Sites auf mehrere Zeilen aufteilen (biallelisch)
# Schritt 2: InDel-Positionen ? Genotypen auf ./. setzen
bcftools norm \
    --threads "$THREADS" \
    --multiallelics -any \
    --fasta-ref "$REF_FASTA" \
    --output-type u \
    "$IN_VCF" | \
bcftools +setGT \
    --threads "$THREADS" \
    -Oz -o "$OUT_VCF" \
    -- \
    -t q \
    -n . \
    -i 'TYPE="indel"'

# Index erstellen
bcftools index --threads "$THREADS" -t "$OUT_VCF"

echo "Fertig: $OUT_VCF"
