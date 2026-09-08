# The script works in the same way for variant calling of the Buckfast samples and SRA samples from UpSet

#!/bin/bash

# 1. Verzeichnisse & Referenz festlegen
INPUT_DIR="/PFAD/"
OUTPUT_DIR="/PFAD/"
REF_GENOME="/PFAD/Amel_HAv3/GCF_003254395.2_Amel_HAv3.1_genomic.fna"

# Ausgabeverzeichnis erstellen
mkdir -p "$OUTPUT_DIR"

# Anzahl der Threads (passe dies an deinen Server an)
THREADS=15

# 2. Schleife über alle dedup BAM-Dateien
for BAM_FILE in "$INPUT_DIR"/*_dedup.bam; do

    # Verhindern, dass die Schleife bei leerem Ordner stolpert
    [ -e "$BAM_FILE" ] || { echo "Keine _dedup.bam Dateien gefunden!"; break; }

    # Probenname extrahieren (z.B. SRR30718962)
    filename=$(basename "$BAM_FILE")
    SAMPLE_NAME="${filename%_dedup.bam}"

    FINAL_VCF="$OUTPUT_DIR/${SAMPLE_NAME}.norm.split.vcf.gz"
    LOG_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_variant_calling.log"

    # Überspringen, falls die fertige VCF schon existiert
    if [ -s "$FINAL_VCF" ]; then
        echo "Fertige VCF für $SAMPLE_NAME existiert bereits. Überspringe."
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starte Variant Calling & Splitting für $SAMPLE_NAME..."

    # 3. Die magische All-in-One Pipeline
    # mpileup (-Ou = unkomprimierter BCF-Output für schnelle Pipe)
    # call    (-Ou = unkomprimierter BCF-Output für schnelle Pipe)
    # norm    (-Oz = komprimierter VCF.gz Output als finales Format)
    bcftools mpileup --threads "$THREADS" -q 30 -Ou -f "$REF_GENOME" "$BAM_FILE" 2> "$LOG_FILE" | \
    bcftools call --threads "$THREADS" -mv -Ou 2>> "$LOG_FILE" | \
    bcftools norm --threads "$THREADS" -f "$REF_GENOME" -m -both -Oz -o "$FINAL_VCF" 2>> "$LOG_FILE"

    # 4. Finales Indexieren (.csi oder .tbi)
    echo "Indexiere finale VCF-Datei..."
    bcftools index --threads "$THREADS" "$FINAL_VCF"

    # 5. Statistiken generieren
    echo "Generiere Statistiken..."
    bcftools stats "$FINAL_VCF" > "$OUTPUT_DIR/${SAMPLE_NAME}_vcf_stats.txt"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Probe $SAMPLE_NAME komplett fertig!"
    echo "---------------------------------------------------"
done

echo "Alle Variant-Calling-Prozesse abgeschlossen."
