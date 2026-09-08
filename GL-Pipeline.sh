#!/bin/bash

export PATH="${HOME}/.cargo/bin:${PATH}"
export LD_LIBRARY_PATH=/PFAD/angsd_env/lib:$LD_LIBRARY_PATH


set -euo pipefail

# ------------------------------------------------------------------------------
# PIPELINE: ANGSD -> ngsLD (LD-Pruning) -> PCAngsd
# Apis mellifera WGS - Populationsstrukturanalyse (run-a3)
# ------------------------------------------------------------------------------

# --- Pfade --------------------------------------------------------------------
REF_BAM_DIR="/PFAD/"
MY_BAM_DIR="/PFAD/"
REF_GENOME="/PFAD/GCF_003254395.2_Amel_HAv3.1_genomic.fna"

OUT_DIR="/PFAD/run-a3"

mkdir -p "${OUT_DIR}"

BAM_LIST="${OUT_DIR}/all_samples.bamlist"
OUT_PREFIX="${OUT_DIR}/all_samples_angsd"
PCANGSD_OUT="${OUT_DIR}/pcangsd_result_run-a3"
LOG="${OUT_DIR}/pipeline_full.log"
THREADS=23

# --- LD-Pruning Parameter (ngsLD) ---------------------------------------------
LD_DIR="${OUT_DIR}/ld_per_chrom"
MAX_KB_DIST=1     # <-- deinen echten Wert aus der Decay-Kurve eintragen
MIN_WEIGHT=0.04    # <-- fuer prune_graph-Schritt
N_PARALLEL_JOBS=4
THREADS_PER_JOB=5  # 4 Jobs x 5 Threads = 20 Threads gesamt, passt zu 22 Kernen

# --- R-Check als Datei schreiben (kein PCA-Plot mehr) -------------------------
R_CHECK="${OUT_DIR}/check_packages.R"

cat > "${R_CHECK}" << 'REOF'
pkgs <- character(0)
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
    cat(sprintf("[INFO] Installiere fehlende R-Pakete: %s\n", paste(missing, collapse=", ")))
    install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
    still_missing <- missing[!sapply(missing, requireNamespace, quietly = TRUE)]
    if (length(still_missing) > 0) {
        stop(sprintf("R-Pakete konnten nicht installiert werden: %s", paste(still_missing, collapse=", ")))
    }
}
cat("[OK] Alle R-Pakete verfuegbar.\n")
REOF

# --- Hilfsfunktionen ----------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }
die() { log "FEHLER: $*"; exit 1; }

skip_if_exists() {
    local file="$1"
    local label="$2"
    if [[ -f "${file}" && -s "${file}" ]]; then
        log "  [SKIP] ${label} - Ausgabe existiert bereits: ${file}"
        return 0
    fi
    return 1
}

# --- SCHRITT 0: Tool-Checks ---------------------------------------------------
log "----------------------------------------"
log "SCHRITT 0: Pruefe benoetigte Tools..."
log "----------------------------------------"


for tool in angsd pcangsd samtools Rscript ngsLD prune_graph; do
    if command -v "${tool}" &>/dev/null; then
        log "  [OK]   ${tool} - $(command -v ${tool})"
    else
        die "${tool} nicht gefunden! Bitte installieren oder PATH pruefen."
    fi
done

log "Pruefe R-Pakete..."
Rscript "${R_CHECK}" 2>&1 | tee -a "${LOG}"
log "Tool-Checks abgeschlossen."

# --- SCHRITT 1: BAM-Liste erstellen (mon_*_MF_*.sorted.bam ausschliessen) -----
log "----------------------------------------"
log "SCHRITT 1: BAM-Liste erstellen (Ausschluss mon_*_MF_*.sorted.bam)..."
log "----------------------------------------"

if skip_if_exists "${BAM_LIST}" "BAM-Liste"; then
    log "  Verwende bestehende Liste mit $(wc -l < "${BAM_LIST}") Samples."
else
    find "${REF_BAM_DIR}" -name "*.bam" ! -name "mon_*_MF_*_dedup.bam" | sort >  "${BAM_LIST}"
    find "${MY_BAM_DIR}"  -name "*.bam" ! -name "mon_*_MF_*_dedup.bam" | sort >> "${BAM_LIST}"
    log "  BAM-Liste erstellt: $(wc -l < "${BAM_LIST}") Samples"
fi
log "  Samples:"
cat "${BAM_LIST}" | tee -a "${LOG}"

# --- SCHRITT 2: BAM-Sortierung & Index pruefen --------------------------------
log "----------------------------------------"
log "SCHRITT 2: Pruefe BAM-Sortierung und Index..."
log "----------------------------------------"

while read -r bam; do
    [[ -z "${bam}" ]] && continue

    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        log "  [WARN] Kein Index fuer $(basename ${bam}) - erstelle Index..."
        samtools index "${bam}" 2>&1 | tee -a "${LOG}"
    else
        log "  [OK]   Index vorhanden: $(basename ${bam})"
    fi

    sort_order=$(samtools view -H "${bam}" | grep "^@HD" | grep -o "SO:[a-z]*" || echo "SO:unknown")
    if [[ "${sort_order}" != "SO:coordinate" ]]; then
        sorted_bam="${bam%.bam}.sorted.bam"
        if skip_if_exists "${sorted_bam}" "Sortierung $(basename ${bam})"; then
            sed -i "s|${bam}|${sorted_bam}|g" "${BAM_LIST}"
        else
            log "  [WARN] $(basename ${bam}) nicht koordinatensortiert (${sort_order}) - sortiere..."
            samtools sort -@ "${THREADS}" -o "${sorted_bam}" "${bam}" 2>&1 | tee -a "${LOG}"
            samtools index "${sorted_bam}" 2>&1 | tee -a "${LOG}"
            sed -i "s|${bam}|${sorted_bam}|g" "${BAM_LIST}"
            log "  [OK]   Sortierung abgeschlossen: $(basename ${sorted_bam})"
        fi
    else
        log "  [OK]   Sortierung OK: $(basename ${bam}) (${sort_order})"
    fi
done < "${BAM_LIST}"

# --- SCHRITT 3: ANGSD - Beagle fuer ngsLD/PCAngsd -----------------------------
log "----------------------------------------"
log "SCHRITT 3: ANGSD - Genotype Likelihoods (Beagle)..."
log "----------------------------------------"

BEAGLE="${OUT_PREFIX}.beagle.gz"

if skip_if_exists "${BEAGLE}" "ANGSD Beagle-Output"; then
    log "  Ueberspringe ANGSD-Lauf."
else
    log "  Starte ANGSD mit ${THREADS} Threads..."
    angsd \
        -bam          "${BAM_LIST}" \
        -ref          "${REF_GENOME}" \
        -out          "${OUT_PREFIX}" \
        -nThreads     "${THREADS}" \
        -GL           1 \
        -doGlf        2 \
        -doMajorMinor 1 \
        -doMaf        1 \
        -SNP_pval     1e-6 \
        -minMapQ      40 \
        -minQ         20 \
        -minInd       57 \
        -remove_bads  1 \
        2>&1 | tee -a "${LOG}"
    log "  ANGSD abgeschlossen."
    N_SITES=$(zcat "${BEAGLE}" | tail -n +2 | wc -l)
    log "  Beagle-File enthaelt ${N_SITES} SNP-Sites."
fi

# --- SCHRITT 4: ngsLD - LD-Pruning --------------------------------------------
log "----------------------------------------"
log "SCHRITT 4: ngsLD - LD-Berechnung & Pruning..."
log "----------------------------------------"

PRUNED_BEAGLE="${OUT_PREFIX}.pruned.beagle.gz"
UNLINKED_POS="${OUT_DIR}/all_samples_unlinked.pos"
COMBINED_LD="${OUT_DIR}/all_samples_combined.ld"

if skip_if_exists "${PRUNED_BEAGLE}" "Gepruntes Beagle-File"; then
    log "  Ueberspringe ngsLD/Pruning."
else
    N_IND=$(wc -l < "${BAM_LIST}")
    log "  N_IND = ${N_IND}"

    if skip_if_exists "${COMBINED_LD}" "ngsLD-Berechnung (combined.ld)"; then
        log "  Ueberspringe ngsLD-Berechnung, verwende bestehende ${COMBINED_LD}."
    else
        CHROM_LIST=$(zcat "${BEAGLE}" | tail -n +2 | cut -f1 | \
            sed -E 's/^(.*)_[0-9]+$/\1/' | sort -u)

        mkdir -p "${LD_DIR}"

        for CHR in ${CHROM_LIST}; do
            (
            CHR_LD="${LD_DIR}/${CHR}.ld"
            if [[ -f "${CHR_LD}" && -s "${CHR_LD}" ]]; then
                log "  [SKIP] ${CHR} - .ld existiert bereits."
            else
                zcat "${BEAGLE}" | \
                    awk -v chr="${CHR}" -F'\t' 'NR==1 || $1 ~ "^"chr"_"' | \
                    gzip > "${LD_DIR}/${CHR}.beagle.gz"

                N_SITES_CHR=$(zcat "${LD_DIR}/${CHR}.beagle.gz" | tail -n +2 | wc -l)

                zcat "${LD_DIR}/${CHR}.beagle.gz" | tail -n +2 | cut -f1 | \
                        sed -E 's/^(.+)_([0-9]+)$/\1\t\2/' > "${LD_DIR}/${CHR}.pos"
                ngsLD \
                    --geno "${LD_DIR}/${CHR}.beagle.gz" \
                    --pos "${LD_DIR}/${CHR}.pos" \
                    --probs \
                    --n_ind "${N_IND}" \
                    --n_sites "${N_SITES_CHR}" \
                    --max_kb_dist "${MAX_KB_DIST}" \
                    --ignore_miss_data \
                    --min_maf 0.01 \
                    --rnd_sample 0.3 \
                    --seed 1 \
                    --n_threads "${THREADS_PER_JOB}" \
                    --out "${CHR_LD}"

                log "  Fertig: ${CHR}"
            fi
            ) &

            while [ "$(jobs -r | wc -l)" -ge "${N_PARALLEL_JOBS}" ]; do
                wait -n
            done
        done

        wait
        log "  Alle Chromosomen fertig."

        FIRST=1
        > "${COMBINED_LD}"
        for f in "${LD_DIR}"/*.ld; do
            if [[ ${FIRST} -eq 1 ]]; then
                cat "${f}" >> "${COMBINED_LD}"
                FIRST=0
            else
                tail -n +2 "${f}" >> "${COMBINED_LD}"
            fi
        done
    fi

    if skip_if_exists "${UNLINKED_POS}" "prune_graph-Ergebnis"; then
        log "  Ueberspringe prune_graph, verwende bestehende ${UNLINKED_POS}."
    else

        prune_graph \
            --header \
            --in "${COMBINED_LD}" \
            --weight "r2_ExpG" \
            --filter "dist <= ${MAX_KB_DIST}000 && r2_ExpG >= ${MIN_WEIGHT}" \
            --out "${UNLINKED_POS}"

        log "  Pruning abgeschlossen: ${UNLINKED_POS}"

    fi

    log "  Filtere Beagle-File auf ungelinkte Sites (via sort+join)..."

    { zcat "${BEAGLE}" || true; } | head -1 > "${OUT_DIR}/pruned_header.tmp"

    sed -E 's/:/_/' "${UNLINKED_POS}" | LC_ALL=C sort -k1,1 > "${OUT_DIR}/unlinked_sorted.tmp"

    zcat "${BEAGLE}" | tail -n +2 | LC_ALL=C sort -k1,1 -S 4G --parallel="${THREADS}" -T "${OUT_DIR}" \
        > "${OUT_DIR}/beagle_sorted.tmp"

    LC_ALL=C join -t $'\t' -1 1 -2 1 "${OUT_DIR}/unlinked_sorted.tmp" "${OUT_DIR}/beagle_sorted.tmp" \
        > "${OUT_DIR}/pruned_body.tmp"

    cat "${OUT_DIR}/pruned_header.tmp" "${OUT_DIR}/pruned_body.tmp" | gzip > "${PRUNED_BEAGLE}"
    rm -f "${OUT_DIR}/pruned_header.tmp" "${OUT_DIR}/pruned_body.tmp" \
          "${OUT_DIR}/unlinked_sorted.tmp" "${OUT_DIR}/beagle_sorted.tmp"

    N_SITES_PRUNED=$(zcat "${PRUNED_BEAGLE}" | tail -n +2 | wc -l)
    log "  Gepruntes Beagle-File enthaelt ${N_SITES_PRUNED} SNP-Sites."
fi

# --- SCHRITT 5: PCAngsd -------------------------------------------------------
log "----------------------------------------"
log "SCHRITT 5: PCAngsd - PCA & Admixture (auf LD-geprunten Sites)..."
log "----------------------------------------"

COV_FILE="${PCANGSD_OUT}.cov"

if skip_if_exists "${COV_FILE}" "PCAngsd Kovarianzmatrix"; then
    log "  Ueberspringe PCAngsd-Lauf."
else
    log "  Starte PCAngsd..."
    pcangsd \
        -b "${PRUNED_BEAGLE}" \
        -o "${PCANGSD_OUT}" \
        -t 22 \
        --selection \
        --admix \
        2>&1 | tee -a "${LOG}"
    log "  PCAngsd abgeschlossen."
    log "  Ausgabedateien:"
    ls -lh "${OUT_DIR}"/pcangsd_result_run-a3* 2>/dev/null | tee -a "${LOG}"
fi

# --- Abschlusszusammenfassung -------------------------------------------------
log "----------------------------------------"
log "PIPELINE ABGESCHLOSSEN"
log "----------------------------------------"
log "Ausgabeverzeichnis: ${OUT_DIR}"
log "Dateien:"
ls -lh "${OUT_DIR}" 2>/dev/null | tee -a "${LOG}"
log "Log gespeichert unter: ${LOG}"
