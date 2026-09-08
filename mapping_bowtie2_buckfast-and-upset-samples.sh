# The script was used in the same way for the Buckfast samples as well as for all SRA samples in the UpSet plot, with the paths adjusted accordingly in each case. 

#!/bin/bash

# 1. Pfade festlegen (angepasst an die HDD04 SRA-Daten)
INPUT_DIR="/PFAD/"
OUTPUT_DIR="/PFAD/"
INDEX="/PFAD/"

# Verzeichnisse erstellen
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/bam_files"

# Anzahl der Threads pro Bowtie2/Samtools-Aufruf
THREADS=20

# 2. Schleife über alle SRA-Forward-Dateien (_1.trimmed.fastq.gz)
for R1 in "$INPUT_DIR"/*_1.trimmed.fastq.gz; do
    
    # Prüfen, ob Dateien gefunden wurden
    [ -e "$R1" ] || { echo "Keine _1.trimmed.fastq.gz Dateien gefunden im Ordner!"; break; }

    # Ableiten des entsprechenden R2-Dateinamens
    R2="${R1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}"

    # Extrahieren des Probenamens (z.B. SRR30718962)
    SAMPLE_NAME=$(basename "$R1" _1.trimmed.fastq.gz)

    # Ausgabedateien definieren
    SORTED_BAM_FILE="$OUTPUT_DIR/bam_files/${SAMPLE_NAME}.sorted.bam"
    LOG_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_bowtie2.log"

    # 3. Überspringen, falls das fertige BAM schon existiert
    if [ -s "$SORTED_BAM_FILE" ]; then
        echo "Fertige BAM-Datei für $SAMPLE_NAME existiert bereits. Überspringe Mapping."
        continue
    fi

    echo "Starte Mapping, Konvertierung und Sortierung für Probe $SAMPLE_NAME..."

    # 4. Der magische Pipeline-Befehl
    # bowtie2 mappt -> gibt SAM aus -> piped an samtools view (macht BAM daraus) -> piped an samtools sort
    # 2> "$LOG_FILE" fängt die Mapping-Statistiken ab (Prozent der gemappten Reads etc.)
    bowtie2 -p "$THREADS" -x "$INDEX" -1 "$R1" -2 "$R2" 2> "$LOG_FILE" | \
    samtools view -bS -@ "$THREADS" - | \
    samtools sort -@ "$THREADS" -o "$SORTED_BAM_FILE"

    echo "Mapping und Sortierung für $SAMPLE_NAME abgeschlossen."

    # 5. BAM-Datei indexieren (.bai Datei erstellen)
    # Das ist zwingend nötig, wenn du die Datei später in IGV anschauen oder weiter filtern willst.
    echo "Erstelle Index (.bai) für $SAMPLE_NAME..."
    samtools index -@ "$THREADS" "$SORTED_BAM_FILE"

    echo "Pipeline für $SAMPLE_NAME komplett fertig!"
    echo "---------------------------------------------------"
done
