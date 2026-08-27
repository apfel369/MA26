#!/bin/bash

# Verzeichnisse festlegen
input_dir="/PFAD ZUR DATEI/"
output_dir="/PFAD ZUR DATEI/"

# Sicherstellen, dass das Ausgabeverzeichnis existiert
mkdir -p "${output_dir}"

# Durch alle Forward-Dateien (_1.fastq.gz) iterieren
for forward_read in "${input_dir}"/_1.fastq.gz; do
    
    # Prüfen, ob überhaupt Dateien gefunden wurden (verhindert Fehler bei leerem Ordner)
    [ -e "$forward_read" ] || { echo "Keine _1.fastq.gz Dateien gefunden!"; break; }

    # Reverse-Read-Datei bestimmen (_1 durch _2 ersetzen)
    reverse_read="${forward_read/_1.fastq.gz/_2.fastq.gz}"

    # Prüfen, ob das Reverse-Read existiert
    if [[ ! -f "$reverse_read" ]]; then
        echo "WARNUNG: Kein passendes Reverse-Read (_2) für $forward_read gefunden!"
        continue
    fi

    # Basisname extrahieren 
    base_name=$(basename "${forward_read}" _1.fastq.gz)

    # Ausgabedateien festlegen
    output_forward="${output_dir}/${base_name}_1.trimmed.fastq.gz"
    output_reverse="${output_dir}/${base_name}_2.trimmed.fastq.gz"
    html_report="${output_dir}/${base_name}_fastp_report.html"
    json_report="${output_dir}/${base_name}_fastp_report.json"

    echo "Starte fastp für ${base_name}..."

    # fastp mit Multi-Threading ausführen (Hinweis: fastp limitiert max Threads auf 16)
    fastp -i "${forward_read}" \
          -I "${reverse_read}" \
          -o "${output_forward}" \
          -O "${output_reverse}" \
          -w 16 \
          -h "${html_report}" \
          -j "${json_report}"

    echo "Fastp abgeschlossen für ${base_name}."
    echo "Report generiert: ${html_report}"
    echo "---------------------------------------------------"
done