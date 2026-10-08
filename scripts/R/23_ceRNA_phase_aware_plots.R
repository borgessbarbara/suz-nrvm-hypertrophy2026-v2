#!/usr/bin/env Rscript
# =============================================================================
# R23 — Figure 6A: phase-aware ceRNA coherence trajectories
#
# PORTED, NOT REWRITTEN.
# Source: archive/.../31_phase_aware_ceRNA_coherence_plots.R (445 lines)
# Produced by scripts/dev/port_from_archive.sh, which applies exactly five
# classes of change and passes everything else through unmodified:
#   1. this preamble (optparse CLI, lib/ sourcing, base_dir/out_dir/dpi)
#   2. hard-coded C:/Users/sebau/... base_dir and out_dir commented out
#   3. base = "C:/Users/sebau/..." default arguments -> base = base_dir
#   4. browser-download duplicates "name(1).csv" removed from candidate lists
#   5. the /mnt/data container search path removed
#
# Diff this file against the archived original to audit the port:
#   diff archive/.../31_phase_aware_ceRNA_coherence_plots.R scripts/R/21_ceRNA_phase_aware_plots.R
# =============================================================================

suppressPackageStartupMessages(library(optparse))

.here <- dirname(normalizePath(sub("^--file=", "",
          commandArgs(trailingOnly = FALSE)[grepl("^--file=",
          commandArgs(trailingOnly = FALSE))][1]), winslash = "/"))
source(file.path(.here, "lib", "io.R"))
source(file.path(.here, "lib", "plotting.R"))

.opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--in-dir", type = "character", default = "results",
                help = "root directory searched for inputs [default %default]"),
    make_option("--out-dir", type = "character", default = "results/R23_ceRNA_phase_aware",
                help = "output directory [default %default]"),
    make_option("--dpi", type = "integer", default = 600,
                help = "raster figure resolution [default %default]")),
  description = "\nR21: Figure 6A: phase-aware ceRNA coherence trajectories\nPorted from 31_phase_aware_ceRNA_coherence_plots.R."))

base_dir <- normalizePath(.opt$`in-dir`, winslash = "/", mustWork = FALSE)
out_dir  <- ensure_out_dir(.opt$`out-dir`)
png_dpi  <- .opt$dpi
dpi_png  <- .opt$dpi
message("[R23] in-dir : ", base_dir)
message("[R23] out-dir: ", out_dir)

# ----------------------------- ported body ----------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================
# Phase-aware expression coherence of filtered ceRNA triplets
#
# Main question:
#   How do interaction-supported ceRNA triplets distribute across
#   Early, Sustained, and Late ncRNA temporal programs?
#
# Strict inclusion rules:
#   - IntaRNA MFE <= -20 kcal/mol
#   - coherent ncRNA-miRNA-mRNA direction
#   - miRNA-mRNA support already present in the integrated triplet table
#
# Recommended main plot:
#   Four panels retain both temporal program and sampling state:
#     Early (6 h)
#     Sustained (6 h)
#     Sustained (24 h)
#     Late (24 h)
#
# This avoids pooling the two sustained measurements into one value.
# circRNA trajectories remain exploratory because their exact variants
# and adjusted differential-expression support are more limited.
# ============================================================

out_dir <- "Figure6A_phase_aware_ceRNA_coherence"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mfe_cutoff <- -20
include_circRNA <- TRUE

# -----------------------------
# Helpers
# -----------------------------
find_input <- function(exact_names, regex_pattern = NULL, required = TRUE) {
  for (nm in exact_names) {
    if (!is.na(nm) && file.exists(nm)) return(nm)
  }
  if (!is.null(regex_pattern)) {
    hits <- list.files(".", recursive = TRUE, full.names = TRUE,
                       pattern = regex_pattern, ignore.case = TRUE)
    if (length(hits) > 0) {
      info <- file.info(hits)
      return(hits[order(info$mtime, decreasing = TRUE)][1])
    }
  }
  if (required) stop("Required input not found: ", paste(exact_names, collapse = ", "))
  NA_character_
}

as_bool <- function(x) toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
numify <- function(x) suppressWarnings(as.numeric(as.character(x)))
normalize_mir <- function(x) tolower(trimws(as.character(x)))

col_chr <- function(df, col, default = NA_character_) {
  if (col %in% names(df)) as.character(df[[col]]) else rep(default, nrow(df))
}
col_num <- function(df, col) {
  if (col %in% names(df)) numify(df[[col]]) else rep(NA_real_, nrow(df))
}

axis_from_directions <- function(nc_dir, mir_dir, mrna_dir) {
  case_when(
    nc_dir == "DOWN" & mir_dir == "UP"   & mrna_dir == "DOWN" ~ "DOWN–UP–DOWN",
    nc_dir == "UP"   & mir_dir == "DOWN" & mrna_dir == "UP"   ~ "UP–DOWN–UP",
    TRUE ~ NA_character_
  )
}

save_dual <- function(plot_obj, filename, width = 12, height = 7.5, dpi = 600) {
  ggsave(file.path(out_dir, paste0(filename, ".pdf")), plot_obj,
         width = width, height = height, units = "in", device = "pdf", bg = "white")
  ggsave(file.path(out_dir, paste0(filename, ".png")), plot_obj,
         width = width, height = height, units = "in", dpi = dpi, bg = "white")
}

axis_colors <- c(
  "DOWN–UP–DOWN" = "#1B9E77",
  "UP–DOWN–UP" = "#D95F02"
)

base_theme <- theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 19, hjust = 0.5),
    plot.subtitle = element_text(size = 11.5, hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11.5),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# -----------------------------
# Inputs
# -----------------------------
lnc_triplet_file <- find_input(
  c("COHERENT_lncRNA_miRNA_mRNA_triplets.csv",
    
    "Fig8_priority_ceRNA_triplets.txt"),
  regex_pattern = "^(COHERENT_lncRNA_miRNA_mRNA_triplets|Fig8_priority_ceRNA_triplets).*(csv|txt)$"
)

lnc6_file <- find_input(
  c("NE_vs_Ctrl_6h_lncRNA_DE.txt"),
  regex_pattern = "^NE_vs_Ctrl_6h_lncRNA_DE.*\\.txt$"
)
lnc24_file <- find_input(
  c("NE_vs_Ctrl_24h_lncRNA_DE.txt"),
  regex_pattern = "^NE_vs_Ctrl_24h_lncRNA_DE.*\\.txt$"
)
mir6_file <- find_input(
  c("DE_miRNA_All_with_status_NE6.txt"),
  regex_pattern = "^DE_miRNA_All_with_status_NE6.*\\.txt$"
)
mir24_file <- find_input(
  c("DE_miRNA_All_with_status_NE24.txt"),
  regex_pattern = "^DE_miRNA_All_with_status_NE24.*\\.txt$"
)

circ_triplet_file <- find_input(
  c("12_STRICT_DE_coherent_circRNA_triplets.csv"),
  regex_pattern = "^12_STRICT_DE_coherent_circRNA_triplets.*\\.csv$",
  required = FALSE
)
circ6_file <- find_input(
  c("06_circRNA_sequences_de6_DE_complete_table.txt"),
  regex_pattern = "^06_circRNA_sequences_de6_DE_complete_table.*\\.txt$",
  required = FALSE
)
circ24_file <- find_input(
  c("07_circRNA_sequences_de24_DE_complete_table.txt"),
  regex_pattern = "^07_circRNA_sequences_de24_DE_complete_table.*\\.txt$",
  required = FALSE
)

# -----------------------------
# DE lookup maps
# -----------------------------
read_lnc_de <- function(path, tp) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  tibble(
    timepoint = tp,
    gene_id = as.character(x$gene_id),
    ncRNA_log2FC = numify(x$log2FoldChange)
  ) %>% distinct(timepoint, gene_id, .keep_all = TRUE)
}

read_mir_de <- function(path, tp) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  mir_col <- if ("miRNA_id" %in% names(x)) "miRNA_id" else if ("miRNA" %in% names(x)) "miRNA" else names(x)[1]
  tibble(
    timepoint = tp,
    miRNA_key = normalize_mir(x[[mir_col]]),
    miRNA_label = as.character(x[[mir_col]]),
    miRNA_log2FC = numify(x$log2FoldChange)
  ) %>% distinct(timepoint, miRNA_key, .keep_all = TRUE)
}

lnc_de <- bind_rows(
  read_lnc_de(lnc6_file, "6h"),
  read_lnc_de(lnc24_file, "24h")
)
mir_de <- bind_rows(
  read_mir_de(mir6_file, "6h"),
  read_mir_de(mir24_file, "24h")
)

# -----------------------------
# lncRNA triplets
# -----------------------------
read_lnc_triplets <- function(path) {
  if (grepl("Fig8_priority_ceRNA_triplets", basename(path), ignore.case = TRUE)) {
    x <- read.delim(path, skip = 2, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  tibble(
    timepoint = as.character(x$timepoint),
    gene_id = as.character(x$gene_id),
    ncRNA_label = if ("ncRNA_label" %in% names(x)) as.character(x$ncRNA_label) else as.character(x$gene_id),
    ncRNA_phase = col_chr(x, "phase"),
    ncRNA_direction = as.character(x$ncRNA_direction),
    miRNA = as.character(x$miRNA),
    miRNA_key = if ("miRNA_key" %in% names(x)) as.character(x$miRNA_key) else normalize_mir(x$miRNA),
    miRNA_direction = as.character(x$miRNA_direction),
    MFE = if ("E" %in% names(x)) numify(x$E) else col_num(x, "best_MFE"),
    mRNA = if ("gene_symbol" %in% names(x)) as.character(x$gene_symbol) else col_chr(x, "mRNA_gene_final"),
    mRNA_direction = if ("mRNA_direction" %in% names(x)) as.character(x$mRNA_direction) else col_chr(x, "mRNA_DE_direction"),
    mRNA_log2FC = col_num(x, "mRNA_log2FC"),
    coherent_triplet = if ("coherent_triplet" %in% names(x)) as_bool(x$coherent_triplet) else rep(TRUE, nrow(x))
  ) %>%
    mutate(
      ncRNA_phase = stringr::str_to_title(trimws(ncRNA_phase)),
      axis = axis_from_directions(ncRNA_direction, miRNA_direction, mRNA_direction)
    ) %>%
    filter(
      coherent_triplet,
      !is.na(MFE), MFE <= mfe_cutoff,
      ncRNA_phase %in% c("Early", "Sustained", "Late"),
      !is.na(axis)
    ) %>%
    distinct(timepoint, gene_id, miRNA_key, mRNA, .keep_all = TRUE) %>%
    left_join(lnc_de, by = c("timepoint", "gene_id")) %>%
    left_join(mir_de, by = c("timepoint", "miRNA_key")) %>%
    mutate(
      ncRNA_type = "lncRNA",
      panel = case_when(
        ncRNA_phase == "Early" & timepoint == "6h" ~ "Early (6 h)",
        ncRNA_phase == "Sustained" & timepoint == "6h" ~ "Sustained (6 h)",
        ncRNA_phase == "Sustained" & timepoint == "24h" ~ "Sustained (24 h)",
        ncRNA_phase == "Late" & timepoint == "24h" ~ "Late (24 h)",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(panel))
}

lnc_trip <- read_lnc_triplets(lnc_triplet_file)

# -----------------------------
# circRNA phase assignment and triplets
# -----------------------------
read_circ_ids <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- if ("circ_id" %in% names(x)) "circ_id" else if ("circRNA_ID" %in% names(x)) "circRNA_ID" else names(x)[1]
  unique(as.character(x[[id_col]]))
}

circ_phase_map <- tibble()
if (include_circRNA && !is.na(circ_triplet_file) && !is.na(circ6_file) && !is.na(circ24_file)) {
  ids6 <- read_circ_ids(circ6_file)
  ids24 <- read_circ_ids(circ24_file)
  all_ids <- union(ids6, ids24)
  circ_phase_map <- tibble(
    circ_locus = all_ids,
    ncRNA_phase = case_when(
      all_ids %in% ids6 & !(all_ids %in% ids24) ~ "Early",
      all_ids %in% ids6 & all_ids %in% ids24 ~ "Sustained",
      !(all_ids %in% ids6) & all_ids %in% ids24 ~ "Late",
      TRUE ~ NA_character_
    )
  )
}

read_circ_triplets <- function(path) {
  if (is.na(path) || !file.exists(path) || nrow(circ_phase_map) == 0) return(tibble())
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  tibble(
    timepoint = as.character(x$timepoint),
    circ_locus = as.character(x$circ_locus),
    ncRNA_label = as.character(x$circ_locus),
    ncRNA_direction = as.character(x$circ_direction),
    ncRNA_log2FC = numify(x$circ_log2FC),
    miRNA = as.character(x$miRNA),
    miRNA_key = as.character(x$miRNA_key),
    miRNA_direction = as.character(x$miRNA_DE_direction),
    miRNA_log2FC = numify(x$miRNA_log2FC),
    MFE = numify(x$best_MFE),
    mRNA = as.character(x$mRNA_gene_final),
    mRNA_direction = as.character(x$mRNA_DE_direction),
    mRNA_log2FC = numify(x$mRNA_log2FC),
    coherent_triplet = as_bool(x$strict_DE_coherent_circRNA_triplet)
  ) %>%
    left_join(circ_phase_map, by = "circ_locus") %>%
    mutate(
      axis = axis_from_directions(ncRNA_direction, miRNA_direction, mRNA_direction),
      ncRNA_type = "circRNA",
      panel = case_when(
        ncRNA_phase == "Early" & timepoint == "6h" ~ "Early (6 h)",
        ncRNA_phase == "Sustained" & timepoint == "6h" ~ "Sustained (6 h)",
        ncRNA_phase == "Sustained" & timepoint == "24h" ~ "Sustained (24 h)",
        ncRNA_phase == "Late" & timepoint == "24h" ~ "Late (24 h)",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      coherent_triplet,
      !is.na(MFE), MFE <= mfe_cutoff,
      !is.na(axis), !is.na(panel)
    ) %>%
    distinct(timepoint, circ_locus, miRNA_key, mRNA, .keep_all = TRUE)
}

circ_trip <- read_circ_triplets(circ_triplet_file)

# -----------------------------
# Unique-node values and summaries
# -----------------------------
lnc_nodes <- bind_rows(
  lnc_trip %>% transmute(timepoint, panel, ncRNA_phase, axis, ncRNA_type, node_class = "ncRNA", feature_id = gene_id, log2FC = ncRNA_log2FC),
  lnc_trip %>% transmute(timepoint, panel, ncRNA_phase, axis, ncRNA_type, node_class = "miRNA", feature_id = miRNA, log2FC = miRNA_log2FC),
  lnc_trip %>% transmute(timepoint, panel, ncRNA_phase, axis, ncRNA_type, node_class = "mRNA", feature_id = mRNA, log2FC = mRNA_log2FC)
)

circ_nodes <- tibble()
if (nrow(circ_trip) > 0) {
  circ_nodes <- bind_rows(
    circ_trip %>% transmute(timepoint, panel, ncRNA_phase, axis, ncRNA_type, node_class = "ncRNA", feature_id = circ_locus, log2FC = ncRNA_log2FC),
    circ_trip %>% transmute(timepoint, panel, ncRNA_phase, axis, ncRNA_type, node_class = "miRNA", feature_id = miRNA, log2FC = miRNA_log2FC),
    circ_trip %>% transmute(timepoint, panel, ncRNA_phase, axis, ncRNA_type, node_class = "mRNA", feature_id = mRNA, log2FC = mRNA_log2FC)
  )
}

all_nodes <- bind_rows(lnc_nodes, circ_nodes) %>%
  filter(!is.na(log2FC), !is.na(feature_id), feature_id != "") %>%
  distinct(timepoint, panel, ncRNA_phase, axis, ncRNA_type, node_class, feature_id, .keep_all = TRUE) %>%
  mutate(
    panel = factor(panel, levels = c("Early (6 h)", "Sustained (6 h)", "Sustained (24 h)", "Late (24 h)")),
    ncRNA_phase = factor(ncRNA_phase, levels = c("Early", "Sustained", "Late")),
    axis = factor(axis, levels = c("DOWN–UP–DOWN", "UP–DOWN–UP")),
    ncRNA_type = factor(ncRNA_type, levels = c("lncRNA", "circRNA")),
    node_class = factor(node_class, levels = c("ncRNA", "miRNA", "mRNA"))
  )

summary_precise <- all_nodes %>%
  group_by(panel, ncRNA_phase, timepoint, axis, ncRNA_type, node_class) %>%
  summarise(
    median_log2FC = median(log2FC, na.rm = TRUE),
    q1 = quantile(log2FC, 0.25, na.rm = TRUE),
    q3 = quantile(log2FC, 0.75, na.rm = TRUE),
    n_unique_nodes = n_distinct(feature_id),
    .groups = "drop"
  )

write.csv(all_nodes, file.path(out_dir, "01_phase_aware_unique_node_values.csv"), row.names = FALSE)
write.csv(summary_precise, file.path(out_dir, "02_phase_timepoint_trajectory_summary_median_IQR.csv"), row.names = FALSE)

triplet_counts <- bind_rows(
  lnc_trip %>% transmute(timepoint, panel, ncRNA_phase, axis, ncRNA_type, triplet_id = paste(gene_id, miRNA_key, mRNA, sep = "__")),
  circ_trip %>% transmute(timepoint, panel, ncRNA_phase, axis, ncRNA_type, triplet_id = paste(circ_locus, miRNA_key, mRNA, sep = "__"))
) %>%
  group_by(panel, ncRNA_phase, timepoint, axis, ncRNA_type) %>%
  summarise(n_unique_triplets = n_distinct(triplet_id), .groups = "drop")
write.csv(triplet_counts, file.path(out_dir, "03_phase_timepoint_triplet_counts.csv"), row.names = FALSE)

# -----------------------------
# Plot 1: recommended precise phase-timepoint view
# -----------------------------
p_precise <- ggplot(
  summary_precise,
  aes(x = node_class, y = median_log2FC, color = axis,
      group = interaction(axis, ncRNA_type),
      linetype = ncRNA_type, shape = ncRNA_type)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.08, linewidth = 0.8) +
  geom_line(linewidth = 1.15) +
  geom_point(size = 3.4) +
  facet_wrap(~panel, ncol = 2, scales = "free_y", drop = FALSE) +
  scale_color_manual(values = axis_colors, drop = FALSE) +
  scale_linetype_manual(values = c("lncRNA" = "solid", "circRNA" = "dashed"), drop = FALSE) +
  scale_shape_manual(values = c("lncRNA" = 16, "circRNA" = 17), drop = FALSE) +
  labs(
    title = "Expression coherence of filtered ceRNA triplets across temporal programs",
    subtitle = "Early and Late are time-specific; Sustained is shown separately at 6 h and 24 h; medians and IQRs use unique nodes",
    x = "Position in the ceRNA hypothesis",
    y = "Median log2FC (NE vs control)",
    color = "Expression axis",
    linetype = "ncRNA type",
    shape = "ncRNA type"
  ) +
  base_theme

save_dual(p_precise, "Figure6A_phase_timepoint_ceRNA_coherence_RECOMMENDED", 12, 8)

# -----------------------------
# Plot 2: phase-only overview (Sustained pools node-timepoint observations)
# -----------------------------
phase_nodes <- all_nodes %>%
  mutate(feature_time_id = paste(timepoint, feature_id, sep = "__")) %>%
  distinct(ncRNA_phase, axis, ncRNA_type, node_class, feature_time_id, .keep_all = TRUE)

summary_phase_only <- phase_nodes %>%
  group_by(ncRNA_phase, axis, ncRNA_type, node_class) %>%
  summarise(
    median_log2FC = median(log2FC, na.rm = TRUE),
    q1 = quantile(log2FC, 0.25, na.rm = TRUE),
    q3 = quantile(log2FC, 0.75, na.rm = TRUE),
    n_unique_node_timepoints = n_distinct(feature_time_id),
    .groups = "drop"
  )
write.csv(summary_phase_only, file.path(out_dir, "04_phase_only_pooled_summary_median_IQR.csv"), row.names = FALSE)

p_phase_only <- ggplot(
  summary_phase_only,
  aes(x = node_class, y = median_log2FC, color = axis,
      group = interaction(axis, ncRNA_type),
      linetype = ncRNA_type, shape = ncRNA_type)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.08, linewidth = 0.8) +
  geom_line(linewidth = 1.15) +
  geom_point(size = 3.4) +
  facet_wrap(~ncRNA_phase, nrow = 1, scales = "free_y", drop = FALSE) +
  scale_color_manual(values = axis_colors, drop = FALSE) +
  scale_linetype_manual(values = c("lncRNA" = "solid", "circRNA" = "dashed"), drop = FALSE) +
  scale_shape_manual(values = c("lncRNA" = 16, "circRNA" = 17), drop = FALSE) +
  labs(
    title = "Phase-based coherence of filtered ceRNA triplets",
    subtitle = "Sustained summarizes both 6 h and 24 h node-timepoint observations; circRNA trajectories are exploratory",
    x = "Position in the ceRNA hypothesis",
    y = "Median log2FC (NE vs control)",
    color = "Expression axis",
    linetype = "ncRNA type",
    shape = "ncRNA type"
  ) +
  base_theme

save_dual(p_phase_only, "SUPP_phase_only_ceRNA_coherence_pooled", 13, 5.8)

# -----------------------------
# Run summary
# -----------------------------
summary_lines <- c(
  "Phase-aware ceRNA coherence analysis",
  "",
  paste0("lncRNA triplet file: ", lnc_triplet_file),
  paste0("circRNA triplet file: ", ifelse(is.na(circ_triplet_file), "not used", circ_triplet_file)),
  paste0("MFE cutoff: <= ", mfe_cutoff),
  "",
  "Recommended main plot:",
  "Figure6A_phase_timepoint_ceRNA_coherence_RECOMMENDED.pdf",
  "",
  "Interpretation:",
  "- Early (6 h): triplets whose ncRNA is Early.",
  "- Sustained (6 h) and Sustained (24 h): the same temporal ncRNA class evaluated separately at each sampling state.",
  "- Late (24 h): triplets whose ncRNA is Late.",
  "- Missing lines indicate that no triplets survived the strict filters for that phase/timepoint/axis/type combination.",
  ""
)
writeLines(summary_lines, file.path(out_dir, "README_phase_aware_plot.txt"))

message("Done. Phase-aware outputs written to: ", out_dir)
