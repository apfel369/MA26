# The script works in the same way for filtering variants in the Buckfast and SRA samples from UpSet, as well as for filtering single-sample VCFs.

# 1. Ausgabeverzeichnis erstellen
# mkdir -p vcf_biall_filt_DP-MQ-MQBZ-RPBZ-QD-GQ
# Feste Grenzen definieren
MIN_DEPTH=15
FLOOR_LIMIT=40
# 2. Schleife starten
for file in *_sorted_amel_dedup.norm.split.vcf; do
    name=$(basename "$file" _sorted_amel_dedup.norm.split.vcf)
    echo "=========================================================="
    echo "Verarbeite: $name"
    # --- A: Dynamische Max-Depth Berechnung ---
    median=$(bcftools query -f '%INFO/DP\n' "$file" | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
    if [ -z "$median" ]; then median=0; fi
    calc_max=$(( median * 2 ))
    if [ "$calc_max" -lt "$FLOOR_LIMIT" ]; then
        final_max=$FLOOR_LIMIT
        echo "  -> Median ($median) niedrig. Max-Depth: $final_max"
    else
        final_max=$calc_max
        echo "  -> Median ($median) okay. Max-Depth: $final_max"
    fi
    # --- B: Filtern (Der Trick mit Perl) ---
    # 1. bcftools view: Filtert alles außer GQ (DP, MQ, RPBZ, QD...)
    # 2. Perl One-Liner: Berechnet GQ aus PL und filtert >= 20
    bcftools view -i "DP >= $MIN_DEPTH && DP <= $final_max && MQ >= 40 && MQBZ > -4 && RPBZ > -5 && RPBZ < 5 && QUAL/DP >= 2" "$file" | \
    perl -ne '
    # Header immer ausgeben
    if (/^#/) { print; next; }
    # Zeile splitten (Tabs)
    @cols = split(/\t/);
    $format = $cols[8];
    $sample = $cols[9];
     # Wo steht "PL" im FORMAT-Feld? (z.B. GT:PL:AD -> Index 1)
    @tags = split(/:/, $format);
    $pl_idx = -1;
    for ($i=0; $i<@tags; $i++) { if ($tags[$i] eq "PL") { $pl_idx=$i; last; } }
  # Wenn PL gefunden wurde:
    if ($pl_idx != -1) {
        @vals = split(/:/, $sample);
        $pl_string = $vals[$pl_idx];
      # PL Werte (z.B. "82,1,0") in Zahlen umwandeln und sortieren
        @pl_nums = split(/,/, $pl_string);
        @sorted = sort { $a <=> $b } @pl_nums;
        # GQ ist der zweitkleinste Wert (Index 1, da Index 0 der kleinste ist)
        # Wir drucken die Zeile nur, wenn GQ >= 20
        if (scalar(@sorted) > 1 && $sorted[1] >= 20) {
            print;
        }
    }
    ' > "vcf_biall_filt_DP-MQ-MQBZ-RPBZ-QD-GQ/${name}_filtered_DP-MQ-MQBZ-RPBZ-QD-GQ.vcf"
    echo "  -> Filterung abgeschlossen."
done
echo "=========================================================="
echo "Fertig! Inklusive manueller GQ >= 20 Berechnung."
