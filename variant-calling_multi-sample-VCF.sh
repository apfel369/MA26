# This script was used to generate multisample VCFs. 
# It was primarily used to generate the VCF that was used in the GT workflow for unsupervised learning methods.

#!/bin/bash

# 1. Verzeichnisse & Referenz festlegen
MY_BAM_DIR="/PFAD/"
REF_BAM_DIR="/PFAD/dedup_bam"
OUTPUT_DIR="/PFAD/"
REF_GENOME="/PFAD/GCF_003254395.2_Amel_HAv3.1_genomic.fna"

mkdir -p "$OUTPUT_DIR"

THREADS=22

BAM_LIST="$OUTPUT_DIR/all_samples.bam.list"
FINAL_VCF="$OUTPUT_DIR/multisample.norm.split.vcf.gz"
LOG_FILE="$OUTPUT_DIR/multisample_variant_calling.log"

# --- 1. BAI-Indices pruefen / erstellen -------------------------------------
echo "Pruefe BAM-Indizes..."
for BAM in "$MY_BAM_DIR"/*_sorted_amel_dedup.bam "$REF_BAM_DIR"/*_dedup.bam; do
    [ -e "$BAM" ] || continue
    if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then
        echo "  Indexiere: $(basename $BAM)"
        samtools index -@ "$THREADS" "$BAM"
    fi
done

# --- 2. BAM-Liste erstellen -------------------------------------------------
ls "$MY_BAM_DIR"/*_sorted_amel_dedup.bam > "$BAM_LIST"
ls "$REF_BAM_DIR"/*_dedup.bam            >> "$BAM_LIST"

echo "BAM-Dateien in der Liste:"
cat "$BAM_LIST"
echo "Gesamtanzahl: $(wc -l < "$BAM_LIST") Samples"

# --- 3. Multi-Sample Variant Calling ----------------------------------------
if [ -s "$FINAL_VCF" ]; then
    echo "Fertige Multi-Sample-VCF existiert bereits. Ueberspringe."
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starte Multi-Sample Variant Calling..."

    bcftools mpileup \
        --threads "$THREADS" \
        -q 30 \
        -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
        -f "$REF_GENOME" \
        -b "$BAM_LIST" \
        -Ou 2> "$LOG_FILE" | \
    bcftools call \
        --threads "$THREADS" \
        -mv \
        -a GQ,GP \
        -Ou 2>> "$LOG_FILE" | \
    bcftools norm \
        --threads "$THREADS" \
        -f "$REF_GENOME" \
        -m -any \
        -Oz -o "$FINAL_VCF" 2>> "$LOG_FILE"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Variant Calling abgeschlossen."
fi

# --- 4. Indexieren ----------------------------------------------------------
echo "Indexiere finale VCF..."
bcftools index --threads "$THREADS" "$FINAL_VCF"

# --- 5. Statistiken ---------------------------------------------------------
echo "Generiere Statistiken..."
bcftools stats "$FINAL_VCF" > "$OUTPUT_DIR/multisample_vcf_stats.txt"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fertig! Multi-Sample-VCF: $FINAL_VCF"
