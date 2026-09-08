# Liftover zu Amel_4.5 zurück: 

# 1. Pfade als Variablen definieren
DIR="/PFAD/"
IN_VCF="$DIR/v2-HBeeID_272_markers_FILTERED_Amelhav31.vcf"
OUT_VCF="$DIR/v2-HBeeID_272_markers_FILTERED_excel_lifted_Amel4.5.vcf"
REJECT_VCF="$DIR/v2-HBeeID_272_markers_FILTERED_excel_REJECTED.vcf"
 
CHAIN="/PFAD/GCF_000002195.4_Amel_4.5_genomic.fna"
 
# 2. Den eigentlichen bcftools liftover Befehl ausführen
bcftools +liftover "$IN_VCF" -o "$OUT_VCF" -- -c "$CHAIN" -f "$TARGET_FASTA" --reject "$REJECT_VCF"

# Jetzt neue VCF machen, die nur die Einträge besitzt wo die ID (letzte Zahl) mit der POS matchen
IN_VCF="/PFAD/v2-HBeeID_272_markers_FILTERED_excel_lifted_Amel4.5.vcf"
OUT_VCF="/PFAD/v2-HBeeID_272_markers_POS_MATCH.vcf"
 
awk -F'\t' '
/^#/ { print $0; next } 
{
    # ID am Unterstrich aufsplitten
    n = split($3, id_parts, "_");
    
    # Pruefen, ob der letzte Teil der ID mit der POS uebereinstimmt
    if (id_parts[n] == $2) {
        print $0;
    }
}' "$IN_VCF" > "$OUT_VCF"
