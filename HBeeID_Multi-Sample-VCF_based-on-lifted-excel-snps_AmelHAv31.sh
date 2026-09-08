#!/bin/bash
# =============================================================================
# HBeeID Multi-Sample VCF aus gelifteter BED (AmelHAv3.1) # =============================================================================

set -euo pipefail

# --- Pfade -------------------------------------------------------------------
REF="/PFAD/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
BED="/PFAD/excel-HBeeID_SNPs-geliftet-to-AmelHAv31.bed"
BAM_DIR="/PFAD/"
OUT_DIR="/PFAD/"
OUT_VCF="${OUT_DIR}/multi-sample-vcf_based-on-excel_amelhav31.vcf"
THREADS=8

# --- Ausgabeverzeichnis erstellen --------------------------------------------
mkdir -p "${OUT_DIR}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Ausgabeverzeichnis erstellt: ${OUT_DIR}"

# --- BED sortieren -----------------------------------------------------------
BED_SORTED="${OUT_DIR}/regions_sorted_uniq.bed"
sort -k1,1 -k2,2n "${BED}" | uniq > "${BED_SORTED}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] BED sortiert: $(wc -l < "${BED_SORTED}") Eintraege"

# --- BAM-Dateien sammeln -----------------------------------------------------
BAM_FILES=$(ls "${BAM_DIR}"/*_sorted_amel_dedup.bam 2>/dev/null)
BAM_COUNT=$(echo "${BAM_FILES}" | wc -l)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Gefundene BAM-Dateien: ${BAM_COUNT}"
echo "${BAM_FILES}" | while read f; do echo "  -> $(basename $f)"; done

if [ "${BAM_COUNT}" -eq 0 ]; then
    echo "FEHLER: Keine BAM-Dateien gefunden in ${BAM_DIR}"
    exit 1
fi

# --- Variant Calling (angepasst: FORMAT/DP und GQ explizit anfordern) --------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starte bcftools mpileup | bcftools call ..."

bcftools mpileup \
    --fasta-ref "${REF}" \
    --regions-file "${BED_SORTED}" \
    --annotate FORMAT/DP \
    --threads "${THREADS}" \
    ${BAM_FILES} \
  | bcftools call \
    --multiallelic-caller \
    --format-fields GQ \
    --output-type u \
  | bcftools norm \
    --multiallelics -any \
    --fasta-ref "${REF}" \
    --output-type v \
    --output "${OUT_VCF}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] VCF erstellt: ${OUT_VCF}"

# --- Anzahl Eintraege zaehlen ------------------------------------------------
TOTAL=$(grep -v '^#' "${OUT_VCF}" | wc -l)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Anzahl Varianten: ${TOTAL}"

# --- Stats -------------------------------------------------------------------
bcftools stats --threads "${THREADS}" "${OUT_VCF}" \
    > "${OUT_DIR}/multi-sample-vcf_based-on-excel_amelhav31_stats.txt"

grep "^SN" "${OUT_DIR}/multi-sample-vcf_based-on-excel_amelhav31_stats.txt" \
    | grep -E "number of (SNPs|MNPs|indels|multiallelic|samples|records)"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fertig."