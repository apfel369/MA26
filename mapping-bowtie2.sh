# The script is designed for the expanded dataset, which includes Buckfast, Monticola, and the samples from the CNP paper. 
# Since Buckfast samples have already been mapped and processed using the analog method, the script is tailored to Monticola and CNP samples.

#!/bin/bash

# 1. Pfade festlegen
INPUT_DIR="/PFAD ZUR DATEI/"
OUTPUT_DIR="/PFAD ZUR DATEI/"
INDEX="/PFAD/"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/bam_files"

THREADS=22

# 2. Schleife ueber alle Forward-Dateien (_1.trimmed.fastq.gz)
for R1 in "$INPUT_DIR"/*_1.trimmed.fastq.gz; do

    [ -e "$R1" ] || { echo "Keine _1.trimmed.fastq.gz Dateien gefunden im Ordner!"; break; }

    R2="${R1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}"

    if [[ ! -f "$R2" ]]; then
        echo "WARNUNG: Kein passendes R2 fuer $R1 gefunden. Ueberspringe."
        continue
    fi

    SAMPLE_NAME=$(basename "$R1" _1.trimmed.fastq.gz)
    SORTED_BAM_FILE="$OUTPUT_DIR/bam_files/${SAMPLE_NAME}.sorted.bam"
    LOG_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_bowtie2.log"

    # 3. Ueberspringen, falls fertiges BAM bereits existiert
    if [ -s "$SORTED_BAM_FILE" ]; then
        echo "Fertige BAM-Datei fuer $SAMPLE_NAME existiert bereits. Ueberspringe Mapping."
        continue
    fi

    echo "Starte Mapping fuer Probe $SAMPLE_NAME..."

    # 4. Datensatz-Zuweisung anhand des Sample-Namens
    if [[ "$SAMPLE_NAME" == ana_tur* || "$SAMPLE_NAME" == car_aut_hun* || \
          "$SAMPLE_NAME" == mel_irl* || "$SAMPLE_NAME" == rut_mlt* ]]; then
        DATASET="CNP0001986"
        PLATFORM_MODEL="HiSeqXTen"
        PLATFORM_UNIT="CNGBdb"
    elif [[ "$SAMPLE_NAME" == mon_* ]]; then
        DATASET="PRJNA357367"
        PLATFORM_MODEL="HiSeq2500"
        PLATFORM_UNIT="SRA"
    else
        echo "WARNUNG: Unbekanntes Sample-Praefix fuer '$SAMPLE_NAME'. RG-Zuweisung fehlgeschlagen."
        exit 1
    fi

    # 5. Mapping-Pipeline
    bowtie2 -p "$THREADS" -x "$INDEX" -1 "$R1" -2 "$R2" \
        --rg-id "${SAMPLE_NAME}_${DATASET}" \
        --rg "SM:$SAMPLE_NAME" \
        --rg "PL:ILLUMINA" \
        --rg "LB:$SAMPLE_NAME" \
        --rg "PU:$PLATFORM_UNIT" \
        --rg "PM:$PLATFORM_MODEL" \
        2> "$LOG_FILE" | \
        samtools view -bS -@ "$THREADS" - | \
        samtools sort -@ "$THREADS" -o "$SORTED_BAM_FILE"

    echo "Mapping und Sortierung fuer $SAMPLE_NAME abgeschlossen."

    # 6. BAM-Datei indexieren
    echo "Erstelle Index (.bai) fuer $SAMPLE_NAME..."
    samtools index -@ "$THREADS" "$SORTED_BAM_FILE"

    echo "Pipeline fuer $SAMPLE_NAME komplett fertig!"
    echo "---------------------------------------------------"
done
