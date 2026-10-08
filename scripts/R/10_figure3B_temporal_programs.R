#!/usr/bin/env Rscript
# =============================================================================
# R10 - Figure 3B: Early, Sustained and Late programmes as stacked bars.
#
# Usage:   Rscript scripts/R/10_figure3B_temporal_programs.R [--in-dir DIR] [--out-dir DIR] [--dpi N]
#          Rscript scripts/R/10_figure3B_temporal_programs.R --help
#
# Adapted from archive/final_package_v3_original/01_FINAL_FIGURE_AND_ANALYSIS_SCRIPTS/36b_Figure3_panelB_DOWN_bottom_UP_top.R
# Changes from the original are listed in docs/CHANGES_FROM_ORIGINALS.md
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
    make_option("--out-dir", type = "character", default = "results/R10_figure3B_temporal_programs",
                help = "output directory [default %default]"),
    make_option("--dpi", type = "integer", default = 600,
                help = "raster figure resolution [default %default]")),
  description = "\nR08: Figure 3B: Early/Sustained/Late stacked bars, DOWN bottom and UP top\nPorted from 36b_Figure3_panelB_DOWN_bottom_UP_top.R."))

base_dir <- normalizePath(.opt$`in-dir`, winslash = "/", mustWork = FALSE)
out_dir  <- ensure_out_dir(.opt$`out-dir`)
.suz_out_root <- out_dir
png_dpi  <- .opt$dpi
dpi_png  <- .opt$dpi
message("[R10] in-dir : ", base_dir)
message("[R10] out-dir: ", out_dir)


suppressPackageStartupMessages({
  library(tidyverse)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

# ============================================================
# Figure 3 / Panel B updater
# - Adds white count labels inside UP/DOWN stacked bars
# - Builds temporal program plots for mRNA, lncRNA, miRNA, circRNA
# - Builds detailed lncRNA annotated vs novel plot with count labels
# - Saves PDF and high-resolution PNG outputs
# ============================================================

out_dir <- ensure_out_dir(file.path(.suz_out_root, "Figure3_panelB_with_counts"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

phase_levels <- c("Early", "Sustained", "Late")
direction_levels <- c("DOWN", "UP")
direction_colors <- c("DOWN" = "#4C78A8", "UP" = "#E07B39")
source_colors <- c("annotated" = "#F07167", "novel" = "#17BEBB")

# -----------------------------
# File finder
# -----------------------------
find_input <- function(exact_names, regex_pattern = NULL, required = TRUE) {
  search_roots <- c(".")

  for (root in search_roots) {
    for (nm in exact_names) {
      candidate <- file.path(root, nm)
      if (file.exists(candidate)) return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  if (!is.null(regex_pattern)) {
    hits <- c()
    for (root in search_roots) {
      if (dir.exists(root)) {
        hits <- c(hits, list.files(root, pattern = regex_pattern, recursive = TRUE,
                                   full.names = TRUE, ignore.case = TRUE))
      }
    }
    hits <- unique(hits[file.exists(hits)])
    if (length(hits) > 0) {
      info <- file.info(hits)
      return(normalizePath(hits[order(info$mtime, decreasing = TRUE)][1], winslash = "/", mustWork = TRUE))
    }
  }

  if (required) stop("Required input not found. Tried: ", paste(exact_names, collapse = ", "))
  NA_character_
}

# -----------------------------
# Inputs
# -----------------------------
mrna6_file <- find_input(c("NE_vs_Ctrl_6h_mRNA_DE.txt"),
                         "^NE_vs_Ctrl_6h_mRNA_DE.*\\.txt$")
mrna24_file <- find_input(c("NE_vs_Ctrl_24h_mRNA_DE.txt"),
                          "^NE_vs_Ctrl_24h_mRNA_DE.*\\.txt$")
lnc_file <- find_input(c("SUPP_Fig4B2_lncDE_with_phase.csv"),
                       "^SUPP_Fig4B2_lncDE_with_phase.*\\.csv$")
mir_file <- find_input(c("miRNA_temporal_summary_table.csv"), "^miRNA_temporal_summary_table.*\\.csv$")
circ6_file <- find_input(c("06_circRNA_sequences_de6_DE_complete_table.txt"),
                         "^06_circRNA_sequences_de6_DE_complete_table.*\\.txt$")
circ24_file <- find_input(c("07_circRNA_sequences_de24_DE_complete_table.txt"),
                          "^07_circRNA_sequences_de24_DE_complete_table.*\\.txt$")

# -----------------------------
# Helpers
# -----------------------------
save_dual <- function(plot_obj, filename_base, width = 4.5, height = 4.6, dpi = 600) {
  ggsave(file.path(out_dir, paste0(filename_base, ".pdf")), plot_obj,
         width = width, height = height, units = "in", device = "pdf", bg = "white")
  ggsave(file.path(out_dir, paste0(filename_base, ".png")), plot_obj,
         width = width, height = height, units = "in", dpi = dpi, bg = "white")
}

base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(5.5, 6, 5.5, 5.5)
  )

normalize_direction <- function(x, logfc = NULL) {
  y <- toupper(trimws(as.character(x)))
  fc <- suppressWarnings(as.numeric(logfc))
  out <- case_when(
    y %in% c("UP", "UPREGULATED") ~ "UP",
    y %in% c("DOWN", "DOWNREGULATED") ~ "DOWN",
    !is.na(fc) & fc > 0 ~ "UP",
    !is.na(fc) & fc < 0 ~ "DOWN",
    TRUE ~ NA_character_
  )
  out
}

as_de_flag <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
}

panel_barplot <- function(df, title, ylab) {
  df2 <- expand_grid(
    phase = factor(phase_levels, levels = phase_levels),
    direction = factor(direction_levels, levels = direction_levels)
  ) %>%
    left_join(df %>% mutate(phase = factor(phase, levels = phase_levels),
                            direction = factor(direction, levels = direction_levels)),
              by = c("phase", "direction")) %>%
    mutate(
      n = replace_na(n, 0L),
      # ggplot stacks factor levels in reverse order by default.
      # Internal order UP, DOWN therefore draws DOWN at the bottom and UP at the top.
      direction_stack = factor(as.character(direction), levels = c("UP", "DOWN"))
    )

  ggplot(df2, aes(x = phase, y = n, fill = direction_stack)) +
    geom_col(width = 0.72, color = NA, position = "stack") +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 3.4) +
    scale_fill_manual(
      values = direction_colors,
      breaks = c("DOWN", "UP"),
      labels = c("DOWN", "UP"),
      drop = FALSE,
      name = "Direction"
    ) +
    labs(title = title, x = NULL, y = ylab) +
    base_theme +
    theme(legend.position = "bottom")
}

lnc_source_plot <- function(df) {
  df2 <- expand_grid(
    phase = factor(phase_levels, levels = phase_levels),
    lnc_source = factor(c("annotated", "novel"), levels = c("annotated", "novel"))
  ) %>%
    left_join(df %>% mutate(phase = factor(phase, levels = phase_levels),
                            lnc_source = factor(lnc_source, levels = c("annotated", "novel"))),
              by = c("phase", "lnc_source")) %>%
    mutate(n = replace_na(n, 0L))

  ggplot(df2, aes(x = phase, y = n, fill = lnc_source)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72, color = NA) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_dodge(width = 0.78),
              vjust = 0.5,
              color = "white", fontface = "bold", size = 3.4) +
    scale_fill_manual(values = source_colors, drop = FALSE) +
    labs(title = "lncRNA temporal programs: annotated vs novel", x = NULL, y = "# lncRNA genes", fill = "source") +
    base_theme +
    theme(legend.position = "right")
}

# -----------------------------
# mRNA counts
# -----------------------------
read_mrna <- function(file, suffix) {
  df <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  df %>%
    mutate(
      gene_id = as.character(gene_id),
      padj = suppressWarnings(as.numeric(padj)),
      log2FoldChange = suppressWarnings(as.numeric(log2FoldChange)),
      sig = !is.na(padj) & padj < 0.05,
      direction = normalize_direction(if ("direction" %in% names(.)) direction else NA_character_, log2FoldChange)
    ) %>%
    filter(sig, !is.na(direction), !is.na(gene_id), gene_id != "") %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    transmute(gene_id, !!paste0("dir_", suffix) := direction)
}

mrna6 <- read_mrna(mrna6_file, "6h")
mrna24 <- read_mrna(mrna24_file, "24h")

mrna_temporal <- full_join(mrna6, mrna24, by = "gene_id") %>%
  mutate(
    phase = case_when(
      !is.na(dir_6h) &  is.na(dir_24h) ~ "Early",
       is.na(dir_6h) & !is.na(dir_24h) ~ "Late",
      !is.na(dir_6h) & !is.na(dir_24h) & dir_6h == dir_24h ~ "Sustained",
      !is.na(dir_6h) & !is.na(dir_24h) & dir_6h != dir_24h ~ "Reversal",
      TRUE ~ NA_character_
    ),
    direction = case_when(
      phase %in% c("Early", "Sustained") ~ dir_6h,
      phase == "Late" ~ dir_24h,
      TRUE ~ NA_character_
    )
  )

mrna_counts <- mrna_temporal %>%
  filter(phase %in% phase_levels, direction %in% direction_levels) %>%
  count(phase, direction, name = "n")

write.csv(mrna_counts, file.path(out_dir, "Figure3B_mRNA_temporal_counts_with_labels.csv"), row.names = FALSE)

# -----------------------------
# lncRNA counts (all lncRNAs combined for panel B)
# -----------------------------
lnc_raw <- read.csv(lnc_file, stringsAsFactors = FALSE, check.names = FALSE)

lnc_temporal <- lnc_raw %>%
  mutate(
    gene_key = as.character(gene_key),
    phase = str_to_title(trimws(as.character(phase))),
    direction = normalize_direction(direction, if ("log2FC" %in% names(.)) log2FC else NA_real_),
    de_keep = if ("de_sig" %in% names(.)) as_de_flag(de_sig) else TRUE,
    lnc_source = case_when(
      str_detect(toupper(gene_key), "MSTRG") ~ "novel",
      TRUE ~ "annotated"
    )
  ) %>%
  filter(de_keep, phase %in% phase_levels, direction %in% direction_levels,
         !is.na(gene_key), gene_key != "") %>%
  distinct(gene_key, phase, direction, .keep_all = TRUE)

lnc_counts <- lnc_temporal %>%
  count(phase, direction, name = "n")

lnc_source_counts <- lnc_temporal %>%
  distinct(gene_key, phase, lnc_source) %>%
  count(phase, lnc_source, name = "n")

write.csv(lnc_counts, file.path(out_dir, "Figure3B_lncRNA_temporal_counts_with_labels.csv"), row.names = FALSE)
write.csv(lnc_source_counts, file.path(out_dir, "Figure3B_lncRNA_annotated_vs_novel_counts_with_labels.csv"), row.names = FALSE)

# -----------------------------
# miRNA counts
# -----------------------------
mir_raw <- read.csv(mir_file, stringsAsFactors = FALSE, check.names = FALSE)

mir_temporal <- mir_raw %>%
  mutate(
    miRNA_id = as.character(miRNA_id),
    temporal_class = as.character(temporal_class),
    temporal_pattern = as.character(temporal_pattern),
    phase = case_when(
      temporal_class == "Early_6h_only" ~ "Early",
      temporal_class == "Sustained_both" ~ "Sustained",
      temporal_class == "Late_24h_only" ~ "Late",
      TRUE ~ NA_character_
    ),
    direction = case_when(
      str_detect(temporal_pattern, "UP") ~ "UP",
      str_detect(temporal_pattern, "DOWN") ~ "DOWN",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(phase %in% phase_levels, direction %in% direction_levels)

mir_counts <- mir_temporal %>%
  count(phase, direction, name = "n")

write.csv(mir_counts, file.path(out_dir, "Figure3B_miRNA_temporal_counts_with_labels.csv"), row.names = FALSE)

# -----------------------------
# circRNA counts
# -----------------------------
read_circ <- function(file, suffix) {
  df <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- if ("circ_id" %in% names(df)) "circ_id" else if ("circRNA_ID" %in% names(df)) "circRNA_ID" else stop("No circRNA ID column found in ", file)

  df %>%
    mutate(
      circ_id = as.character(.data[[id_col]]),
      pvalue = suppressWarnings(as.numeric(pvalue)),
      log2FC = suppressWarnings(as.numeric(log2FC)),
      length = if ("length" %in% names(.)) suppressWarnings(as.numeric(length)) else NA_real_,
      direction = normalize_direction(if ("regulations" %in% names(.)) regulations else NA_character_, log2FC),
      status = case_when(
        !is.na(pvalue) & pvalue <= 0.05 & log2FC >= 0.5 ~ "UP",
        !is.na(pvalue) & pvalue <= 0.05 & log2FC <= -0.5 ~ "DOWN",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(circ_id, pvalue, desc(length)) %>%
    distinct(circ_id, .keep_all = TRUE) %>%
    filter(!is.na(status)) %>%
    transmute(circ_id, !!paste0("dir_", suffix) := status)
}

circ6 <- read_circ(circ6_file, "6h")
circ24 <- read_circ(circ24_file, "24h")

circ_temporal <- full_join(circ6, circ24, by = "circ_id") %>%
  mutate(
    phase = case_when(
      !is.na(dir_6h) &  is.na(dir_24h) ~ "Early",
       is.na(dir_6h) & !is.na(dir_24h) ~ "Late",
      !is.na(dir_6h) & !is.na(dir_24h) & dir_6h == dir_24h ~ "Sustained",
      !is.na(dir_6h) & !is.na(dir_24h) & dir_6h != dir_24h ~ "Reversal",
      TRUE ~ NA_character_
    ),
    direction = case_when(
      phase %in% c("Early", "Sustained") ~ dir_6h,
      phase == "Late" ~ dir_24h,
      TRUE ~ NA_character_
    )
  )

circ_counts <- circ_temporal %>%
  filter(phase %in% phase_levels, direction %in% direction_levels) %>%
  count(phase, direction, name = "n")

write.csv(circ_counts, file.path(out_dir, "Figure3B_circRNA_temporal_counts_with_labels.csv"), row.names = FALSE)

# -----------------------------
# Build plots
# -----------------------------
p_mrna <- panel_barplot(mrna_counts, "mRNA temporal programs", "# mRNA genes")
p_lnc  <- panel_barplot(lnc_counts,  "lncRNA temporal programs", "# lncRNA genes")
p_mir  <- panel_barplot(mir_counts,  "miRNA temporal programs", "# miRNAs")
p_circ <- panel_barplot(circ_counts, "circRNA temporal programs", "# circRNA loci")
p_lnc_src <- lnc_source_plot(lnc_source_counts)

save_dual(p_mrna, "Figure3B_mRNA_temporal_programs_with_labels")
save_dual(p_lnc,  "Figure3B_lncRNA_temporal_programs_with_labels")
save_dual(p_mir,  "Figure3B_miRNA_temporal_programs_with_labels")
save_dual(p_circ, "Figure3B_circRNA_temporal_programs_with_labels")
save_dual(p_lnc_src, "Figure3B_lncRNA_annotated_vs_novel_with_labels", width = 6.6, height = 5.2)

if (has_patchwork) {
  p_combined <- (p_mrna + p_lnc + p_mir + p_circ) +
    plot_layout(ncol = 4, guides = "collect") & theme(legend.position = "bottom")

  ggsave(file.path(out_dir, "Figure3B_temporal_programs_all_RNAs_with_labels.pdf"), p_combined,
         width = 13.8, height = 4.8, units = "in", device = "pdf", bg = "white")
  ggsave(file.path(out_dir, "Figure3B_temporal_programs_all_RNAs_with_labels.png"), p_combined,
         width = 13.8, height = 4.8, units = "in", dpi = 600, bg = "white")
}

# -----------------------------
# Session / manifest
# -----------------------------
manifest <- tibble(
  item = c("mrna6_file", "mrna24_file", "lnc_file", "mir_file", "circ6_file", "circ24_file"),
  path = c(mrna6_file, mrna24_file, lnc_file, mir_file, circ6_file, circ24_file)
)
write.csv(manifest, file.path(out_dir, "input_manifest.csv"), row.names = FALSE)

message("Done. Outputs written to: ", out_dir)
