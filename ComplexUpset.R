# ==============================================================================
# UpSet Plot - RAM-minimal, farbige Punkte via upset_query, Set Sizes links
# Basis: v3-DegreeSorted, erweitert um SRR30718962 + SRR30719121
# Visuell angepasst an Referenzplot (grafik-2.jpg)
# ==============================================================================

library(data.table)
library(ggplot2)
library(ComplexUpset)
library(scales)

output_dir  <- "/PFAD/"
output_file <- paste0(output_dir, "/v3-ComplexUpset.svg")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

input_dir <- "/PFAD/"

vcf_files_dir <- list.files(path = input_dir, pattern = "_only-snps_sorted\\.vcf\\.gz$", full.names = TRUE)
if (length(vcf_files_dir) == 0) stop("Keine VCF-Dateien im input_dir gefunden!")

vcf_files_extra <- c(
  "/PFAD/SRR30718962-biall-snvs-only_filtered_HighQuality.vcf.gz",
  "/PFAD/SRR30719121-biall-snvs-only_filtered_HighQuality.vcf.gz"
)

extra_sample_names <- gsub("-biall.*\\.vcf\\.gz$", "", basename(vcf_files_extra))
options(scipen = 999)

# ==============================================================================
# Farbzuweisung
# ==============================================================================
probe_colors <- c(
  "B1-1"        = "#5C4E8A",
  "B13018-3"    = "#7A9FD4",
  "B13018-4"    = "#6ECFE8",
  "B20-6"       = "#6DDBB8",
  "B26-3"       = "#A8D97F",
  "B34-4"       = "#D4C96A",
  "B40-6"       = "#E8B86D",
  "B57-1"       = "#E08050",
  "B6-1"        = "#C05A45",
  "B6-3"        = "#8B4040",
  "SRR30718962" = "#C97BB2",
  "SRR30719121" = "#E8A8D0"
)

# ==============================================================================
# SCHRITT 1: VCF einlesen
# ==============================================================================
cat("Lese VCF-Dateien ein...\n")
variant_list <- list()
sample_names <- c()

for (vcf in vcf_files_dir) {
  sample_name  <- gsub("_only-snps_sorted\\.vcf\\.gz$", "", basename(vcf))
  sample_names <- c(sample_names, sample_name)
  cat("  -> Probe (dir):", sample_name, "\n")
  cmd      <- paste0("zcat ", vcf, " | grep -v '^#' | awk '{print $1\"_\"$2\"_\"$5\"_\"$4}'")
  variants <- fread(cmd = cmd, header = FALSE, col.names = "ID")
  variant_list[[sample_name]] <- variants$ID
  rm(variants); gc()
}

for (i in seq_along(vcf_files_extra)) {
  vcf          <- vcf_files_extra[i]
  sample_name  <- extra_sample_names[i]
  sample_names <- c(sample_names, sample_name)
  cat("  -> Probe (extra):", sample_name, "\n")
  cmd      <- paste0("zcat ", vcf, " | grep -v '^#' | awk '{print $1\"_\"$2\"_\"$5\"_\"$4}'")
  variants <- fread(cmd = cmd, header = FALSE, col.names = "ID")
  variant_list[[sample_name]] <- variants$ID
  rm(variants); gc()
}

# ==============================================================================
# SCHRITT 2: Aggregieren
# ==============================================================================
cat("Erstelle Binary-Matrix und aggregiere...\n")
dt_long <- rbindlist(lapply(names(variant_list), function(s) {
  data.table(VariantID = variant_list[[s]], Sample = s, Presence = TRUE)
}))
rm(variant_list); gc()

upset_wide <- dcast(dt_long, VariantID ~ Sample, value.var = "Presence", fill = FALSE)
rm(dt_long); gc()

upset_wide[, VariantID := NULL]

all_patterns <- upset_wide[, .(intersect_size = .N), by = sample_names]
all_patterns[, degree := rowSums(.SD), .SDcols = sample_names]

rm(upset_wide); gc()

# ==============================================================================
# Set Sizes (wahre Variantenzahlen pro Probe)
# ==============================================================================
cat("Berechne wahre Set Sizes...\n")
true_set_sizes <- setNames(
  sapply(sample_names, function(s) {
    sum(all_patterns$intersect_size[as.logical(all_patterns[[s]])])
  }),
  sample_names
)
cat("Set Sizes:\n")
print(true_set_sizes)
sample_names_sorted <- names(sort(true_set_sizes, decreasing = FALSE))

# ==============================================================================
# inclusive_union und ratio
# ==============================================================================
cat("Berechne inclusive_union und ratio...\n")
pattern_matrix_full <- as.matrix(all_patterns[, ..sample_names])

inclusive_union_vec <- sapply(seq_len(nrow(pattern_matrix_full)), function(i) {
  active_cols <- which(pattern_matrix_full[i, ])
  if (length(active_cols) == 0L) return(0L)
  any_active <- rowSums(pattern_matrix_full[, active_cols, drop = FALSE]) > 0
  sum(all_patterns$intersect_size[any_active])
})

all_patterns[, inclusive_union := as.numeric(inclusive_union_vec)]
all_patterns[, ratio := as.numeric(intersect_size) / inclusive_union]

rm(pattern_matrix_full, inclusive_union_vec); gc()


# ==============================================================================
# AUSGABE: Top 5 Degree-2 Schnittmengen
# ==============================================================================
cat("\n=== Top 5 Schnittmengen Degree 4 ===\n")

#deg2 <- all_patterns[degree == 2][order(-intersect_size)][1:min(5, .N)]
deg4 <- all_patterns[degree == 4][order(-intersect_size)][1:min(5, .N)]

for (i in seq_len(nrow(deg4))) {
  proben  <- sample_names[as.logical(as.matrix(deg4[i, ..sample_names]))]
  isect   <- deg4$intersect_size[i]
  union   <- deg4$inclusive_union[i]
  ratio   <- deg4$ratio[i]

  cat(sprintf(
    "\n[%d] %s\n    Schnittmenge : %s SNVs\n    Union        : %s SNVs\n    Anteil       : %.2f %%\n",    i,
    paste(proben, collapse = "  &  "),
    format(isect, big.mark = ".", decimal.mark = ",", scientific = FALSE),
    format(union, big.mark = ".", decimal.mark = ",", scientific = FALSE),
    ratio * 100
  ))
}

stop("Diagnose abgeschlossen - kein Plot erstellt.")


# --- DIAGNOSE ---
cat("\n=== Top 3 Intersektionen je Degree (ab Degree 6) ===\n")
for (d in 6:length(sample_names)) {
  sub <- all_patterns[degree == d][order(-intersect_size)][1:min(3, .N)]
  if (nrow(sub) == 0) next
  cat(sprintf("\nDegree %d (%d Muster gesamt):\n", d, all_patterns[degree == d, .N]))
  for (i in seq_len(nrow(sub))) {
    aktive <- paste(sample_names[as.logical(as.matrix(sub[i, ..sample_names]))], collapse = " | ")
    cat(sprintf("  [%d] intersect_size = %d  |  ratio = %.4f  |  Proben: %s\n",
                i, sub$intersect_size[i], sub$ratio[i], aktive))
  }
}

# ==============================================================================
# SCHRITT 3: Ziel-Schnittmengen selektieren
# ==============================================================================
cat("Selektiere Ziel-Muster...\n")
target_idx <- integer(0)

idx1 <- all_patterns[degree == 1, which = TRUE]
target_idx <- c(target_idx, idx1[order(-all_patterns$intersect_size[idx1])])

max_degree <- length(sample_names)
for (d in 2:(max_degree - 1)) {
  idx_d <- all_patterns[degree == d, which = TRUE]
  if (length(idx_d) > 0) {
    idx_d <- idx_d[order(-all_patterns$intersect_size[idx_d])]
    target_idx <- c(target_idx, head(idx_d, 3))
  }
}

idx_max <- all_patterns[degree == max_degree, which = TRUE]
if (length(idx_max) > 0) {
  target_idx <- c(target_idx, idx_max[order(-all_patterns$intersect_size[idx_max])])
}

target_patterns <- all_patterns[target_idx]
rm(all_patterns); gc()

# ==============================================================================
# SCHRITT 4: upset_df
# ==============================================================================
cat("Erstelle kompakten upset_df...\n")

upset_df <- as.data.frame(target_patterns[, ..sample_names])
for (s in sample_names) upset_df[[s]] <- as.logical(upset_df[[s]])

upset_df$intersect_size    <- target_patterns$intersect_size
upset_df$precomputed_union <- target_patterns$inclusive_union
upset_df$precomputed_ratio <- target_patterns$ratio

gc()

intersect_list <- lapply(seq_len(nrow(target_patterns)), function(i) {
  sample_names[as.logical(as.matrix(target_patterns[i, ..sample_names]))]
})

# ==============================================================================
# SCHRITT 5: Plot
# ==============================================================================
cat("Erstelle Plot...\n")

balken_breite <- 0.65

# Set Sizes als statischer Datensatz mit echten Variantenzahlen
ss_df <- data.frame(
  x    = factor(sample_names_sorted, levels = sample_names_sorted),
  y    = as.numeric(true_set_sizes[sample_names_sorted]),
  fill = probe_colors[sample_names_sorted]
)

color_queries <- lapply(sample_names_sorted, function(s) {
  upset_query(
    set             = s,
    fill            = probe_colors[s],
    color           = probe_colors[s],
    only_components = c("intersections_matrix")
  )
})

p <- upset(
  upset_df,
  intersect          = sample_names_sorted,
  intersections      = intersect_list,
  sort_sets          = FALSE,
  sort_intersections = FALSE,
  name               = "Proben",
  height_ratio       = 0.6,
  encode_sets        = FALSE,

  queries = color_queries,

  matrix = intersection_matrix(
    geom    = geom_point(size = 2.5),
    segment = geom_segment(linewidth = 0.5)
  ),

  set_sizes = (
    upset_set_size(
      geom = geom_col(
        data        = ss_df,
        mapping     = aes(x = x, y = y, fill = x),
        width       = balken_breite,
        stat        = "identity",
        inherit.aes = FALSE
      ),
      position             = "left",
      filter_intersections = FALSE
    )
    + scale_fill_manual(values = setNames(ss_df$fill, ss_df$x), guide = "none")
    + geom_text(
        data        = ss_df,
        mapping     = aes(
          x     = x,
          y     = y * 0.85,
          label = format(round(y / 1000, 1), nsmall = 1, decimal.mark = ",")
        ),
        angle       = 90,
        hjust       = 1.0,
        size        = 3.0,
        color       = "white",
        inherit.aes = FALSE
      )
    + scale_y_continuous(
        name   = "Anzahl SNVs [10\u00B3]",
        labels = function(x) format(x / 1000, big.mark = ".", decimal.mark = ","),
        trans  = "reverse",
        expand = expansion(mult = c(0, 0.05))
      )
    + theme(
        axis.title.x = element_text(size = 11),
        axis.text.x  = element_text(size = 9),
        axis.text.y  = element_text(size = 9)
      )
  ),

  base_annotations = list(
    'relative Schnittmenge [%]' = (
      ggplot(mapping = aes(x = intersection))
      + stat_summary(
          mapping = aes(y = precomputed_ratio),
          fun     = "mean",
          geom    = "bar",
          fill    = "grey40",
          width   = balken_breite
        )
      + stat_summary(
          mapping = aes(
            y     = precomputed_ratio,
            label = paste0(format(round(after_stat(y) * 100, 1), nsmall = 1, decimal.mark = ","), "%")
          ),
          fun     = "mean",
          geom    = "text",
          vjust   = -0.4,
          size    = 3.0
        )
      + geom_text(
          mapping = aes(
            y     = 0.008,
            label = format(precomputed_union, big.mark = ".", decimal.mark = ",", scientific = FALSE)
          ),
          angle   = 90,
          hjust   = 0,
          size    = 2.5,
          color   = "white"
        )

      + scale_y_continuous(
    name   = "relative Schnittmenge [%]",
    labels = percent_format(accuracy = 0.1, decimal.mark = ","),
    limits = c(0, 0.30),                    # <-- fest bis 30%
    breaks = seq(0, 0.30, by = 0.10),       # Striche bei 0%, 10%, 20%, 30%
    expand = expansion(mult = c(0, 0.02))   # kaum Puffer oben
  )

      + theme(
          axis.title.y = element_text(size = 11),
          axis.text.y  = element_text(size = 9),
          axis.text.x  = element_blank()
        )
    )
  ),

  annotations = list(
    'absolute Schnittmenge [10\u00B3]' = (
      ggplot(mapping = aes(x = intersection))
      + stat_summary(
          mapping = aes(y = intersect_size),
          fun     = "mean",
          geom    = "bar",
          fill    = "grey30",
          width   = balken_breite
        )
      + stat_summary(
          mapping = aes(
            y     = intersect_size,
            label = format(round(after_stat(y) / 1000, 1), nsmall = 1, decimal.mark = ",")
          ),
          fun     = "mean",
          geom    = "text",
          vjust   = -0.5,
          size    = 3.2
        )
      + scale_y_continuous(
          name   = "absolute Schnittmenge [10\u00B3]",
          labels = function(x) format(x / 1000, big.mark = ".", decimal.mark = ","),
          expand = expansion(mult = c(0, 0.15))
        )
      + theme(
          axis.title.y = element_text(size = 11),
          axis.text.y  = element_text(size = 9),
          axis.text.x  = element_blank()
        )
    )
  ),

  themes = upset_modify_themes(list(
    'intersections_matrix' = theme(
        axis.text.y  = element_text(size = 10),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
    ),
    'overall_sizes' = theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y  = element_text(size = 9),
        axis.title.x = element_text(size = 11)
    )
  ))
)

# ==============================================================================
# SCHRITT 6: Export
# ==============================================================================
cat("Speichere als SVG...\n")

svg(filename = output_file, width = 16, height = 10, bg = "transparent")
print(p)
dev.off()

message("Fertig! Plot gespeichert unter:\n", output_file)
