#!/usr/bin/env bash
BCF_DIR="/home/Drives/HDD04_06T_SDF/Buckfast_08.2024/AmelHAv3.1/bcftool_vcf_outputv2/"
OUTFILE="/home/Drives/HDD04_06T_SDF/Buckfast_08.2024/AmelHAv3.1/bcftool_vcf_outputv2/all_stats_combined.txt"

echo "# bcftools stats results for all samples" > "$OUTFILE"

for f in ${BCF_DIR}/*.bcf; do
  base=$(basename "$f" .bcf)
  echo "" >> "$OUTFILE"
  echo "=== Sample: $base ===" >> "$OUTFILE"
  bcftools stats "$f" >> "$OUTFILE"
done





