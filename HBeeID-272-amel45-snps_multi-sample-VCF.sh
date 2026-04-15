#!/bin/bash

# ==========================================
# 1. Pfade, Dateien und Threads definieren
# ==========================================
THREADS=8  # <--- HIER die Anzahl der gewünschten Threads eintragen

REF_FASTA="/PFAD ZUR DATEI/Amel_4.5/GCF_000002195.4_Amel_4.5_genomic.fna"
BAM_DIR="/PFAD ZUR DATEI/Amel_4.5/mapped_data/bam_files/dedup-picard_bam"
MARKER_FILE="/PFAD ZUR DATEI/HBeeID_SNP_data/HBeeID_SNPs_CHROM-POS.txt"
OUT_DIR="/PFAD ZUR DATEI/Amel_4.5/HBeeID-multi-sample-vcf-272markers_Amel45"

# Name der finalen VCF-Datei (entsprechend Amel_4.5 benannt)
OUT_VCF="${OUT_DIR}/HBeeID_272_markers_Amel45_multisample.vcf"
BAM_LIST="${OUT_DIR}/bam_list.txt"

# ==========================================
# 2. Sicherheitsprüfungen
# ==========================================
echo "Starte Targeted Joint Calling (Amel_4.5) mit $THREADS Threads..."
echo "---------------------------------------------------------"

if [ ! -f "$REF_FASTA" ]; then
    echo "Fehler: Referenzgenom nicht gefunden: $REF_FASTA"
    exit 1
fi

if [ ! -f "$MARKER_FILE" ]; then
    echo "Fehler: Marker-Datei nicht gefunden: $MARKER_FILE"
    exit 1
fi

# Ausgabeverzeichnis erstellen, falls nicht vorhanden
mkdir -p "$OUT_DIR"

# Liste aller passenden BAM-Dateien erstellen (geändertes Namensmuster!)
ls "${BAM_DIR}"/*_sorted_amel45_dedup.bam > "$BAM_LIST" 2>/dev/null

# Prüfen, ob BAMs gefunden wurden
BAM_COUNT=$(wc -l < "$BAM_LIST")
if [ "$BAM_COUNT" -eq 0 ]; then
    echo "Fehler: Keine BAM-Dateien mit dem Muster *_sorted_amel45_dedup.bam gefunden!"
    exit 1
fi

echo "Gefundene BAM-Dateien: $BAM_COUNT"
echo "Referenzgenom: $(basename "$REF_FASTA")"
echo "Marker-Datei: $(basename "$MARKER_FILE")"
echo "Ausgabe-VCF: $(basename "$OUT_VCF")"
echo "---------------------------------------------------------"

# ==========================================
# 3. Das Variant Calling (mpileup + call)
# ==========================================
echo "Führe mpileup und call aus. Bitte warten..."

bcftools mpileup \
    --threads "$THREADS" \
    -f "$REF_FASTA" \
    -T "$MARKER_FILE" \
    -a FORMAT/DP,FORMAT/AD \
    -Ou \
    -b "$BAM_LIST" | \
bcftools call \
    --threads "$THREADS" \
    -m -O v -o "$OUT_VCF"

echo "---------------------------------------------------------"
echo "Fertig! Die Multi-Sample VCF für Amel_4.5 wurde erfolgreich erstellt."
echo "Zu finden unter: $OUT_VCF"