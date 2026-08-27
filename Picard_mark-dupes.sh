# The script was used for all Picard MarkDuplicates runs. 
# The only changes made were to adjust paths and directories as needed for each dataset. 

#!/bin/bash

# 1. Pfade festlegen
input_dir="/PFAD/"
rg_dir="$input_dir/"
dedup_dir="$input_dir/"

PICARD="/PFAD/picard"

# Ausgabe-Verzeichnisse erstellen
mkdir -p "$rg_dir"
mkdir -p "$dedup_dir"

# 2. Durch alle sortierten BAM-Dateien iterieren
for bam_file in "$input_dir"/*.sorted.bam; do

    [ -e "$bam_file" ] || { echo "Keine .sorted.bam Dateien gefunden!"; break; }

    filename=$(basename "$bam_file")
    sample_name="${filename%.sorted.bam}"

    echo "Verarbeite Probe: $sample_name"

    # Standardmaessig geht das Original-BAM in den Dedup-Schritt
    input_for_dedup="$bam_file"

    # 3. Schritt 1: Pruefen ob Read Groups (@RG) vorhanden sind
    if samtools view -H "$bam_file" | grep -q "^@RG"; then
        echo "-> Read Group (@RG) ist bereits vorhanden. Ueberspringe AddOrReplaceReadGroups."
    else
        echo "-> Keine Read Group gefunden. Fuege Read Groups hinzu..."

        rg_bam="$rg_dir/${sample_name}_RG.bam"

        _JAVA_OPTIONS="-Xmx16G" "$PICARD" AddOrReplaceReadGroups \
            I="$bam_file" \
            O="$rg_bam" \
            RGID="${sample_name}_ref" \
            RGLB="$sample_name" \
            RGPL=ILLUMINA \
            RGSM="$sample_name"

        input_for_dedup="$rg_bam"
    fi

    # 4. Schritt 2: Duplikate markieren (MarkDuplicates)
    final_bam="$dedup_dir/${sample_name}_dedup.bam"

    # Ueberspringen falls Dedup-BAM bereits existiert
    if [ -s "$final_bam" ]; then
        echo "-> Dedup-BAM existiert bereits. Ueberspringe MarkDuplicates."
        continue
    fi

    echo "-> Markiere Duplikate..."

    _JAVA_OPTIONS="-Xmx16G" "$PICARD" MarkDuplicates \
        I="$input_for_dedup" \
        O="$final_bam" \
        M="$dedup_dir/${sample_name}_dedup.metrics.txt"

    # 5. Schritt 3: Indexieren der finalen BAM-Datei
    echo "-> Erstelle Index (.bai) fuer $sample_name..."
    samtools index -@ 15 "$final_bam"

    echo "Probe $sample_name komplett verarbeitet!"
    echo "---------------------------------------------------"
done
