# The script works in the same way for filtering variants in the Buckfast and SRA samples from UpSet, as well as for filtering single-sample VCFs.

#!/bin/bash

# 1. Ein- und Ausgabeverzeichnisse festlegen
INPUT_DIR="/PFAD/"
OUTPUT_DIR="/PFAD/"

mkdir -p "$OUTPUT_DIR"

# Feste Grenzen definieren
MIN_DEPTH=15
FLOOR_LIMIT=40

echo "Start der Filter-Pipeline..."

# 2. Schleife über die gezipten VCF-Dateien starten
for file in "$INPUT_DIR"/*.vcf.gz; do

    # Check, ob Dateien existieren
    [ -e "$file" ] || { echo "Keine .vcf.gz Dateien gefunden!"; break; }

    # Dateinamen bereinigen (entfernt die Endung für einen sauberen Output-Namen)
    filename=$(basename "$file" .vcf.gz)
    # Falls noch ".sorted" im Namen ist, nehmen wir das auch weg
    name="${filename%.sorted}"

    echo "=========================================================="
    echo "Verarbeite: $name"

    # --- A: Dynamische Max-Depth Berechnung ---
    # bcftools query liest die DP-Werte aus, awk berechnet den Median
    median=$(bcftools query -f '%INFO/DP\n' "$file" | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
    
    # Fallback, falls kein Median gefunden wurde
    if [ -z "$median" ]; then median=0; fi
    
    calc_max=$(( median * 2 ))
    
    if [ "$calc_max" -lt "$FLOOR_LIMIT" ]; then
        final_max=$FLOOR_LIMIT
        echo "  -> Median ($median) ist niedrig. Setze Max-Depth auf: $final_max"
    else
        final_max=$calc_max
        echo "  -> Median ($median) ist stark. Setze Max-Depth auf: $final_max"
    fi

    # Ausgabedatei definieren
    OUTPUT_FILE="$OUTPUT_DIR/${name}_filtered_HighQuality_v2.vcf.gz"

    # --- B: Filtern (bcftools + Perl-Trick) ---
    echo "  -> Filtere Werte und berechne GQ aus PL..."

# 1. bcftools view: Filtert alles außer GQ (DP, MQ, RPBZ, MQBZ, QD berechnet als QUAL/DP)
#    Fehlende Bias-Werte (MQBZ, RPBZ) werden als neutral behandelt (Site wird behalten)
# 2. Perl: Berechnet GQ aus PL und druckt nur Zeilen mit GQ >= 20
# 3. bgzip: Komprimiert das Ergebnis direkt
    bcftools view -i "
        DP >= $MIN_DEPTH && DP <= $final_max &&
        MQ >= 40 &&
        (MQBZ = \".\" || MQBZ > -4) &&
        (RPBZ = \".\" || (RPBZ > -5 && RPBZ < 5)) &&
        QUAL/DP >= 2
    " "$file" | \
    perl -ne '
    # Header immer ausgeben
    if (/^#/) { print; next; }

    # Zeile am Tabulator splitten
    @cols = split(/\t/);
    $format = $cols[8];
    $sample = $cols[9];
    
    # Wo steht "PL" im FORMAT-Feld?
    @tags = split(/:/, $format);
    $pl_idx = -1;
    for ($i=0; $i<@tags; $i++) { if ($tags[$i] eq "PL") { $pl_idx=$i; last; } }
    
    # Wenn PL gefunden wurde:
    if ($pl_idx != -1) {
        @vals = split(/:/, $sample);
        $pl_string = $vals[$pl_idx];
        
        # PL Werte in Zahlen umwandeln und sortieren
        @pl_nums = split(/,/, $pl_string);
        @sorted = sort { $a <=> $b } @pl_nums;
        
        # GQ ist der zweitkleinste Wert. Drucken, wenn >= 20
        if (scalar(@sorted) > 1 && $sorted[1] >= 20) {
            print;
        }
    }
    ' | bgzip -c > "$OUTPUT_FILE"

    # --- C: Finale Indexierung ---
    echo "  -> Erstelle Index für die gefilterte Datei..."
    bcftools index "$OUTPUT_FILE"
    
    echo "  -> Abgeschlossen für $name!"
done

echo "=========================================================="
echo "Komplette Filterung erfolgreich beendet. Daten sind bereit für UpSet-Plots!"
