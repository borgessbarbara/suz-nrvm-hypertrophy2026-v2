#!/usr/bin/env Rscript
# =============================================================================
# R09 - Figures 3A, 3C miRNA, 5A and S4A.
#
# Up/down DE counts across all four RNA classes, the miRNA temporal heatmap,
# the lncRNA guilt-by-association map, and the curated KEGG representation
# among significant lncRNA-mRNA correlations.
#
# Usage:   Rscript scripts/R/09_de_counts_heatmaps_gba.R [--in-dir DIR] [--out-dir DIR] [--dpi N]
#          Rscript scripts/R/09_de_counts_heatmaps_gba.R --help
#
# Adapted from archive/final_package_v3_original/01_FINAL_FIGURE_AND_ANALYSIS_SCRIPTS/35_FINAL_missing_figures_miRNA_heatmap_DE_GBA.R
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
    make_option("--out-dir", type = "character", default = "results/R09_de_counts_heatmaps_gba",
                help = "output directory [default %default]"),
    make_option("--dpi", type = "integer", default = 600,
                help = "raster figure resolution [default %default]")),
  description = "\nR07: Figure 3A, Figure 3C miRNA, Figure 5A and Figure S4B: DE counts, miRNA temporal heatmap and lncRNA guilt-by-association\nPorted from 35_FINAL_missing_figures_miRNA_heatmap_DE_GBA.R."))

base_dir <- normalizePath(.opt$`in-dir`, winslash = "/", mustWork = FALSE)
out_dir  <- ensure_out_dir(.opt$`out-dir`)
.suz_out_root <- out_dir
png_dpi  <- .opt$dpi
dpi_png  <- .opt$dpi
message("[R09] in-dir : ", base_dir)
message("[R09] out-dir: ", out_dir)


suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

out_dir <- ensure_out_dir(file.path(.suz_out_root, "Final_missing_figures_FINAL"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# User-adjustable settings
# -----------------------------
top_miRNAs_per_phase <- 12
top_lncRNAs_per_state <- 12
rho_cutoff <- 0.80
correlation_p_cutoff <- 0.05
gba_temporal_mode <- "time_consistent"  # "time_consistent" or "all_selected"
gba_require_mRNA_DE <- TRUE

phase_levels <- c("Early", "Sustained", "Late")
direction_levels <- c("UP", "DOWN")
rna_levels <- c("mRNA", "lncRNA", "miRNA", "circRNA")

direction_colors <- c("UP" = "#D95F02", "DOWN" = "#2C7FB8")
phase_colors <- c("Early" = "#E6AB02", "Sustained" = "#4D4D4D", "Late" = "#7B3294")
source_shapes <- c("annotated" = 16, "novel" = 17)

selected_kegg <- tribble(
  ~pathway_id, ~pathway, ~kegg_program,
  "rno04060", "Cytokine-cytokine receptor interaction", "6 h state-only",
  "rno04510", "Focal adhesion", "6 h state-only",
  "rno04350", "TGF-beta signaling pathway", "6 h state-only",
  "rno04390", "Hippo signaling pathway", "Shared",
  "rno04670", "Leukocyte transendothelial migration", "Shared",
  "rno04141", "Protein processing in endoplasmic reticulum", "Shared",
  "rno04530", "Tight junction", "Shared",
  "rno04110", "Cell cycle", "24 h state-only",
  "rno04820", "Cytoskeleton in muscle cells", "24 h state-only",
  "rno00565", "Ether lipid metabolism", "24 h state-only",
  "rno03320", "PPAR signaling pathway", "24 h state-only"
)

# -----------------------------
# Helpers
# -----------------------------
find_input <- function(exact_names, regex_pattern = NULL, required = TRUE) {
  for (nm in exact_names) if (!is.na(nm) && file.exists(nm)) return(nm)
  if (!is.null(regex_pattern)) {
    hits <- list.files(".", recursive = TRUE, full.names = TRUE,
                       pattern = regex_pattern, ignore.case = TRUE)
    if (length(hits) > 0) {
      info <- file.info(hits)
      return(hits[order(info$mtime, decreasing = TRUE)][1])
    }
  }
  if (required) stop("Required file not found: ", paste(exact_names, collapse = ", "))
  NA_character_
}

read_flex <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

save_dual <- function(plot_obj, filename, width = 10, height = 7, dpi = 600) {
  ggsave(file.path(out_dir, paste0(filename, ".pdf")), plot_obj,
         width = width, height = height, units = "in", device = "pdf", bg = "white")
  ggsave(file.path(out_dir, paste0(filename, ".png")), plot_obj,
         width = width, height = height, units = "in", dpi = dpi, bg = "white")
}

normalize_direction <- function(x, logfc = NULL) {
  y <- toupper(trimws(as.character(x)))
  fc <- if (is.null(logfc)) rep(NA_real_, length(y)) else suppressWarnings(as.numeric(logfc))
  case_when(
    y %in% c("UP", "UPREGULATED", "SIGUP") ~ "UP",
    y %in% c("DOWN", "DOWNREGULATED", "SIGDOWN") ~ "DOWN",
    !is.na(fc) & fc > 0 ~ "UP",
    !is.na(fc) & fc < 0 ~ "DOWN",
    TRUE ~ NA_character_
  )
}

z_score <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

classify_two_states <- function(tbl6, tbl24, rna_class) {
  full_join(
    tbl6 %>% select(feature_id, dir_6h = direction),
    tbl24 %>% select(feature_id, dir_24h = direction),
    by = "feature_id"
  ) %>%
    mutate(
      phase = case_when(
        !is.na(dir_6h) & is.na(dir_24h) ~ "Early",
        is.na(dir_6h) & !is.na(dir_24h) ~ "Late",
        !is.na(dir_6h) & !is.na(dir_24h) & dir_6h == dir_24h ~ "Sustained",
        !is.na(dir_6h) & !is.na(dir_24h) & dir_6h != dir_24h ~ "Direction_switch",
        TRUE ~ NA_character_
      ),
      direction = case_when(
        phase == "Early" ~ dir_6h,
        phase == "Late" ~ dir_24h,
        phase == "Sustained" ~ dir_6h,
        TRUE ~ NA_character_
      ),
      direction_switch = phase == "Direction_switch",
      rna_class = rna_class
    )
}

base_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# -----------------------------
# Locate inputs
# -----------------------------
mir_temporal_file <- find_input(c("miRNA_temporal_summary_table.csv"),
                                "^miRNA_temporal_summary_table.*\\.csv$")
mir6_file <- find_input(c("DE_miRNA_All_with_status_NE6.txt"),
                        "^DE_miRNA_All_with_status_NE6.*\\.txt$")
mir24_file <- find_input(c("DE_miRNA_All_with_status_NE24.txt"),
                         "^DE_miRNA_All_with_status_NE24.*\\.txt$")

mrna6_file <- find_input(c("NE_vs_Ctrl_6h_mRNA_DE.txt"),
                         "^NE_vs_Ctrl_6h_mRNA_DE.*\\.txt$")
mrna24_file <- find_input(c("NE_vs_Ctrl_24h_mRNA_DE.txt"),
                          "^NE_vs_Ctrl_24h_mRNA_DE.*\\.txt$")
lnc6_file <- find_input(c("NE_vs_Ctrl_6h_lncRNA_DE.txt"),
                        "^NE_vs_Ctrl_6h_lncRNA_DE.*\\.txt$")
lnc24_file <- find_input(c("NE_vs_Ctrl_24h_lncRNA_DE.txt"),
                         "^NE_vs_Ctrl_24h_lncRNA_DE.*\\.txt$")
lnc_phase_file <- find_input(
  c("SUPP_Fig4B2_lncDE_with_phase.csv"),
  "^SUPP_Fig4B2_lncDE_with_phase.*\\.csv$"
)

circ6_file <- find_input(c("06_circRNA_sequences_de6_DE_complete_table.txt"),
                         "^06_circRNA_sequences_de6_DE_complete_table.*\\.txt$")
circ24_file <- find_input(c("07_circRNA_sequences_de24_DE_complete_table.txt"),
                          "^07_circRNA_sequences_de24_DE_complete_table.*\\.txt$")

cor6_file <- find_input(c("cor_early_annotated.csv", "correlations_lost_after_6h.csv"),
                        "^(cor_early_annotated|correlations_lost_after_6h).*\\.csv$")
cor24_file <- find_input(c("cor_late_annotated.csv", "correlations_gained_at_24h.csv", "correlations_late_only.csv"),
                         "^(cor_late_annotated|correlations_gained_at_24h|correlations_late_only).*\\.csv$")

kegg6_file <- find_input(c("KEGG_6h_vs_24h_Only6h_long_format_SYMBOLS.txt"),
                         "^KEGG_6h_vs_24h_Only6h_long_format_SYMBOLS.*\\.txt$")
kegg_shared_file <- find_input(c("KEGG_6h_vs_24h_Shared_long_format_SYMBOLS.txt"),
                               "^KEGG_6h_vs_24h_Shared_long_format_SYMBOLS.*\\.txt$")
kegg24_file <- find_input(c("KEGG_6h_vs_24h_Only24h_long_format_SYMBOLS.txt"),
                          "^KEGG_6h_vs_24h_Only24h_long_format_SYMBOLS.*\\.txt$")

# =====================================================================
# FIGURE 1: miRNA temporal heatmap in the same format as the existing
# mRNA/lncRNA/circRNA temporal heatmaps
# =====================================================================
mir_temporal <- read.csv(mir_temporal_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  transmute(
    miRNA_id = as.character(miRNA_id),
    log2FC_6h = suppressWarnings(as.numeric(log2FC_6h)),
    log2FC_24h = suppressWarnings(as.numeric(log2FC_24h)),
    phase = case_when(
      temporal_class == "Early_6h_only" ~ "Early",
      temporal_class == "Sustained_both" ~ "Sustained",
      temporal_class == "Late_24h_only" ~ "Late",
      TRUE ~ NA_character_
    ),
    direction = case_when(
      str_detect(toupper(as.character(temporal_pattern)), "UP") ~ "UP",
      str_detect(toupper(as.character(temporal_pattern)), "DOWN") ~ "DOWN",
      TRUE ~ NA_character_
    ),
    rank_score = pmax(abs(log2FC_6h), abs(log2FC_24h), na.rm = TRUE)
  ) %>%
  filter(phase %in% phase_levels, direction %in% direction_levels) %>%
  mutate(
    # Grey cells reproduce the temporal-class logic used in the other heatmaps:
    # Early = only the 6 h value is displayed; Late = only the 24 h value.
    log2FC_6h_plot = if_else(phase == "Late", NA_real_, log2FC_6h),
    log2FC_24h_plot = if_else(phase == "Early", NA_real_, log2FC_24h),
    phase = factor(phase, levels = phase_levels),
    direction = factor(direction, levels = direction_levels)
  )

mir_top <- mir_temporal %>%
  group_by(phase) %>%
  slice_max(rank_score, n = top_miRNAs_per_phase, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(phase, desc(rank_score), miRNA_id)

write.csv(mir_top,
          file.path(out_dir, "Figure1_miRNA_temporal_heatmap_TOP_data.csv"),
          row.names = FALSE)

make_miRNA_temporal_heatmap <- function(df, filename) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install it with install.packages('pheatmap').")
  }

  df <- as.data.frame(df, stringsAsFactors = FALSE)
  row_ids <- make.unique(as.character(df$miRNA_id))

  mat <- cbind(
    log2FC_6h = suppressWarnings(as.numeric(df$log2FC_6h_plot)),
    log2FC_24h = suppressWarnings(as.numeric(df$log2FC_24h_plot))
  )
  rownames(mat) <- row_ids
  storage.mode(mat) <- "double"

  finite_values <- mat[is.finite(mat)]
  if (length(finite_values) == 0) {
    stop("The miRNA heatmap matrix contains no finite numeric log2FC values.")
  }

  heat_limit <- max(4, ceiling(max(abs(finite_values))))
  heat_breaks <- seq(-heat_limit, heat_limit, length.out = 101)
  heat_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

  annotation_row <- data.frame(
    phase = factor(as.character(df$phase), levels = phase_levels),
    direction = factor(as.character(df$direction), levels = direction_levels),
    stringsAsFactors = FALSE,
    row.names = row_ids
  )

  annotation_colors <- list(
    phase = c("Early" = "#7CAE00", "Sustained" = "#F8766D", "Late" = "#00BFC4"),
    direction = c("UP" = "#D95F02", "DOWN" = "#2C7FB8")
  )

  phase_counts <- table(factor(df$phase, levels = phase_levels))
  gaps_row <- cumsum(as.integer(phase_counts))
  gaps_row <- gaps_row[gaps_row < nrow(mat)]

  draw_heatmap <- function() {
    pheatmap::pheatmap(
      mat,
      color = heat_colors,
      breaks = heat_breaks,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      annotation_row = annotation_row,
      annotation_colors = annotation_colors,
      gaps_row = gaps_row,
      fontsize_row = 8,
      fontsize_col = 11,
      main = "Top miRNAs by response class",
      border_color = "grey70",
      na_col = "grey80",
      angle_col = 90,
      silent = TRUE
    )
  }

  hm <- draw_heatmap()
  figure_height <- max(8, 0.28 * nrow(mat) + 2.4)

  grDevices::pdf(file.path(out_dir, paste0(filename, ".pdf")),
                 width = 8.5, height = figure_height, onefile = FALSE)
  grid::grid.newpage()
  grid::grid.draw(hm$gtable)
  grDevices::dev.off()

  grDevices::png(file.path(out_dir, paste0(filename, ".png")),
                 width = 8.5, height = figure_height,
                 units = "in", res = 600, bg = "white")
  grid::grid.newpage()
  grid::grid.draw(hm$gtable)
  grDevices::dev.off()
}

make_miRNA_temporal_heatmap(
  mir_top,
  "Figure1_miRNA_temporal_heatmap_TOP"
)

# =====================================================================
# FIGURE 2: global DE RNA summary and temporal composition
# =====================================================================
read_gene_de <- function(path, rna_class, timepoint) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  x %>%
    transmute(
      feature_id = as.character(gene_id),
      direction = normalize_direction(direction, log2FoldChange),
      rna_class = rna_class,
      timepoint = timepoint
    ) %>%
    filter(!is.na(feature_id), feature_id != "", direction %in% direction_levels) %>%
    distinct(feature_id, .keep_all = TRUE)
}

read_mir_de <- function(path, timepoint) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  x %>%
    transmute(
      feature_id = as.character(miRNA_id),
      direction = case_when(
        toupper(trimws(as.character(status))) == "UPREGULATED" ~ "UP",
        toupper(trimws(as.character(status))) == "DOWNREGULATED" ~ "DOWN",
        TRUE ~ NA_character_
      ),
      rna_class = "miRNA",
      timepoint = timepoint
    ) %>%
    filter(direction %in% direction_levels) %>%
    distinct(feature_id, .keep_all = TRUE)
}

read_circ_de <- function(path, timepoint) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- if ("circ_id" %in% names(x)) "circ_id" else if ("circRNA_ID" %in% names(x)) "circRNA_ID" else names(x)[1]
  x %>%
    transmute(
      feature_id = as.character(.data[[id_col]]),
      direction = normalize_direction(regulations, log2FC),
      rna_class = "circRNA",
      timepoint = timepoint
    ) %>%
    filter(!is.na(feature_id), feature_id != "", direction %in% direction_levels) %>%
    group_by(feature_id) %>%
    summarise(direction = first(direction), rna_class = first(rna_class), timepoint = first(timepoint), .groups = "drop")
}

mrna6 <- read_gene_de(mrna6_file, "mRNA", "6h")
mrna24 <- read_gene_de(mrna24_file, "mRNA", "24h")
lnc6 <- read_gene_de(lnc6_file, "lncRNA", "6h")
lnc24 <- read_gene_de(lnc24_file, "lncRNA", "24h")
mir6 <- read_mir_de(mir6_file, "6h")
mir24 <- read_mir_de(mir24_file, "24h")
circ6 <- read_circ_de(circ6_file, "6h")
circ24 <- read_circ_de(circ24_file, "24h")

mrna_temporal <- classify_two_states(mrna6, mrna24, "mRNA")
circ_temporal <- classify_two_states(circ6, circ24, "circRNA")

lnc_phase_raw <- read.csv(lnc_phase_file, stringsAsFactors = FALSE, check.names = FALSE)
lnc_temporal <- lnc_phase_raw %>%
  transmute(
    feature_id = as.character(gene_key),
    phase = str_to_title(as.character(phase)),
    direction = normalize_direction(direction, log2FC),
    rna_class = "lncRNA"
  ) %>%
  filter(phase %in% phase_levels, direction %in% direction_levels) %>%
  group_by(feature_id, phase, rna_class) %>%
  summarise(
    direction_switch = n_distinct(direction) > 1,
    direction = if_else(direction_switch, NA_character_, first(direction)),
    .groups = "drop"
  ) %>%
  mutate(phase = if_else(direction_switch, "Direction_switch", phase))

mir_temporal_class <- mir_temporal %>%
  transmute(
    feature_id = miRNA_id,
    phase,
    direction,
    direction_switch = FALSE,
    rna_class = "miRNA"
  )

temporal_all <- bind_rows(
  mrna_temporal %>% select(feature_id, phase, direction, direction_switch, rna_class),
  lnc_temporal %>% select(feature_id, phase, direction, direction_switch, rna_class),
  mir_temporal_class,
  circ_temporal %>% select(feature_id, phase, direction, direction_switch, rna_class)
) %>%
  mutate(
    rna_class = factor(rna_class, levels = rna_levels),
    phase = factor(phase, levels = c(phase_levels, "Direction_switch")),
    direction = factor(direction, levels = direction_levels)
  )

write.csv(temporal_all, file.path(out_dir, "Figure2_all_RNA_temporal_classification.csv"), row.names = FALSE)
write.csv(filter(temporal_all, direction_switch),
          file.path(out_dir, "Figure2_direction_switch_features_excluded_from_directional_bars.csv"),
          row.names = FALSE)

# Panel A: unique DE features across the experiment, counted once
union_counts <- temporal_all %>%
  filter(!direction_switch, phase %in% phase_levels, !is.na(direction)) %>%
  count(rna_class, direction, name = "n")
write.csv(union_counts, file.path(out_dir, "Figure2A_unique_DE_RNA_counts_by_class_direction.csv"), row.names = FALSE)

p_union <- ggplot(union_counts, aes(rna_class, n, fill = direction)) +
  geom_col(width = 0.72) +
  facet_wrap(~rna_class, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = direction_colors, drop = FALSE) +
  labs(
    title = "Differentially expressed RNA features across the experiment",
    x = NULL, y = "Number of unique features", fill = "Direction"
  ) +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_dual(p_union, "Figure2A_overall_unique_DE_RNA_counts", 11, 5.3)

# Panel B: Early/Sustained/Late counts
phase_counts <- temporal_all %>%
  filter(!direction_switch, phase %in% phase_levels, !is.na(direction)) %>%
  count(rna_class, phase, direction, name = "n", .drop = FALSE)
write.csv(phase_counts, file.path(out_dir, "Figure2B_temporal_program_counts_by_RNA_class_direction.csv"), row.names = FALSE)

p_phase <- ggplot(phase_counts, aes(phase, n, fill = direction)) +
  geom_col(width = 0.72) +
  facet_wrap(~rna_class, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = direction_colors, drop = FALSE) +
  labs(
    title = "Temporal programs of differentially expressed RNA classes",
    x = "Temporal program", y = "Number of features", fill = "Direction"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_dual(p_phase, "Figure2B_DE_RNA_temporal_programs", 12, 5.8)

# Panel C: timepoint-specific counts
state_counts <- bind_rows(mrna6, mrna24, lnc6, lnc24, mir6, mir24, circ6, circ24) %>%
  count(rna_class, timepoint, direction, name = "n") %>%
  mutate(
    rna_class = factor(rna_class, levels = rna_levels),
    timepoint = factor(timepoint, levels = c("6h", "24h")),
    direction = factor(direction, levels = direction_levels)
  )
write.csv(state_counts, file.path(out_dir, "Figure2C_timepoint_DE_RNA_counts_by_class_direction.csv"), row.names = FALSE)

p_state <- ggplot(state_counts, aes(timepoint, n, fill = direction)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.8) +
  facet_wrap(~rna_class, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = direction_colors, drop = FALSE) +
  labs(
    title = "Differentially expressed RNA classes at the 6 h and 24 h states",
    x = NULL, y = "Number of features", fill = "Direction"
  ) +
  base_theme
save_dual(p_state, "Figure2C_DE_RNA_counts_6h_24h", 11, 5.3)

if (requireNamespace("patchwork", quietly = TRUE)) {
  combined_de <- patchwork::wrap_plots(list(p_union, p_phase, p_state), ncol = 1) +
    patchwork::plot_annotation(tag_levels = "A")
  save_dual(combined_de, "Figure2_COMBINED_DE_RNA_overview", 12, 15.5)
}

# =====================================================================
# FIGURE 3: lncRNA guilt-by-association with exact curated KEGG terms
# =====================================================================
read_kegg_reference <- function(path, program) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  tibble(
    pathway_id = as.character(x$ID),
    pathway = as.character(x$Description),
    kegg_program = program,
    gene_symbol = as.character(x$gene_symbol),
    gene_symbol_key = toupper(trimws(as.character(x$gene_symbol))),
    NES_6h = if ("NES_6h" %in% names(x)) as.numeric(x$NES_6h) else if (program == "6 h state-only") as.numeric(x$NES) else NA_real_,
    NES_24h = if ("NES_24h" %in% names(x)) as.numeric(x$NES_24h) else if (program == "24 h state-only") as.numeric(x$NES) else NA_real_
  ) %>% distinct()
}

kegg_reference <- bind_rows(
  read_kegg_reference(kegg6_file, "6 h state-only"),
  read_kegg_reference(kegg_shared_file, "Shared"),
  read_kegg_reference(kegg24_file, "24 h state-only")
) %>%
  semi_join(selected_kegg, by = c("pathway_id", "pathway", "kegg_program"))

read_correlations <- function(path, state_label) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  x %>%
    transmute(
      correlation_state = state_label,
      lncRNA = as.character(lncRNA),
      lnc_gene_name = as.character(lnc_gene_name),
      lnc_source = tolower(as.character(lnc_source)),
      mRNA_id = as.character(mRNA),
      mRNA_symbol = as.character(mRNA_gene_name),
      mRNA_symbol_key = toupper(trimws(as.character(mRNA_gene_name))),
      rho = as.numeric(rho),
      pvalue = as.numeric(pvalue)
    ) %>%
    filter(
      !is.na(lncRNA), lncRNA != "",
      !is.na(mRNA_symbol_key), mRNA_symbol_key != "",
      !is.na(rho), abs(rho) >= rho_cutoff,
      !is.na(pvalue), pvalue < correlation_p_cutoff
    ) %>%
    distinct(correlation_state, lncRNA, mRNA_id, .keep_all = TRUE)
}

cor_all <- bind_rows(
  read_correlations(cor6_file, "6 h"),
  read_correlations(cor24_file, "24 h")
)

lnc_phase_map <- lnc_phase_raw %>%
  transmute(
    lncRNA = as.character(gene_key),
    lnc_phase = str_to_title(as.character(phase)),
    lnc_source_phase = if_else(str_detect(toupper(as.character(gene_key)), "^MSTRG"), "novel", tolower(as.character(lnc_source)))
  ) %>%
  filter(lnc_phase %in% phase_levels) %>%
  distinct(lncRNA, .keep_all = TRUE)

if (gba_require_mRNA_DE) {
  de_mrna_keys <- bind_rows(
    read.delim(mrna6_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
      transmute(correlation_state = "6 h", mRNA_symbol_key = toupper(trimws(as.character(gene_name)))),
    read.delim(mrna24_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
      transmute(correlation_state = "24 h", mRNA_symbol_key = toupper(trimws(as.character(gene_name))))
  ) %>% distinct()
  cor_all <- semi_join(cor_all, de_mrna_keys, by = c("correlation_state", "mRNA_symbol_key"))
}

# Base merge is used deliberately because one mRNA can belong to multiple exact KEGG pathways.
gba_map <- merge(
  cor_all,
  as.data.frame(kegg_reference),
  by.x = "mRNA_symbol_key",
  by.y = "gene_symbol_key",
  all = FALSE,
  sort = FALSE
) %>%
  as_tibble() %>%
  left_join(lnc_phase_map, by = "lncRNA") %>%
  mutate(
    lnc_source = case_when(
      str_detect(toupper(lncRNA), "^MSTRG") ~ "novel",
      lnc_source %in% c("novel", "annotated") ~ lnc_source,
      TRUE ~ coalesce(lnc_source_phase, "annotated")
    ),
    correlation_sign = if_else(rho > 0, "Positive", "Negative"),
    lnc_label = case_when(
      str_detect(toupper(lncRNA), "^MSTRG") ~ lncRNA,
      !is.na(lnc_gene_name) & lnc_gene_name != "" & lnc_gene_name != "NA" ~ lnc_gene_name,
      TRUE ~ lncRNA
    )
  )

if (gba_temporal_mode == "time_consistent") {
  gba_map <- gba_map %>%
    filter(
      (correlation_state == "6 h" & kegg_program %in% c("6 h state-only", "Shared")) |
      (correlation_state == "24 h" & kegg_program %in% c("Shared", "24 h state-only"))
    )
} else if (gba_temporal_mode != "all_selected") {
  stop("gba_temporal_mode must be 'time_consistent' or 'all_selected'.")
}

gba_map <- gba_map %>%
  distinct(correlation_state, lncRNA, mRNA_id, pathway_id, .keep_all = TRUE)

write.csv(gba_map,
          file.path(out_dir, "Figure3_GBA_exact_lncRNA_mRNA_KEGG_mapping.csv"),
          row.names = FALSE)

# Audit: which selected pathways are represented in each correlation state?
gba_pathway_audit <- tidyr::crossing(
  tibble(correlation_state = c("6 h", "24 h")),
  selected_kegg
) %>%
  left_join(
    gba_map %>%
      group_by(correlation_state, pathway_id, pathway, kegg_program) %>%
      summarise(
        n_lncRNAs = n_distinct(lncRNA),
        n_mRNAs = n_distinct(mRNA_symbol_key),
        n_pairs = n_distinct(paste(lncRNA, mRNA_id, sep = "__")),
        terminal_mRNAs = paste(sort(unique(mRNA_symbol)), collapse = ";"),
        .groups = "drop"
      ),
    by = c("correlation_state", "pathway_id", "pathway", "kegg_program")
  ) %>%
  mutate(
    across(c(n_lncRNAs, n_mRNAs, n_pairs), ~replace_na(.x, 0L)),
    represented = n_pairs > 0
  )
write.csv(gba_pathway_audit,
          file.path(out_dir, "Figure3_GBA_selected_KEGG_pathway_representation_audit.csv"),
          row.names = FALSE)

# Rank lncRNAs by curated-pathway coverage, separately for 6 h and 24 h.
lnc_rank <- gba_map %>%
  group_by(correlation_state, lncRNA, lnc_label, lnc_source, lnc_phase) %>%
  summarise(
    n_curated_mRNAs = n_distinct(mRNA_symbol_key),
    n_curated_pathways = n_distinct(pathway_id),
    n_correlated_pairs = n_distinct(paste(lncRNA, mRNA_id, sep = "__")),
    .groups = "drop"
  ) %>%
  group_by(correlation_state) %>%
  arrange(desc(n_curated_pathways), desc(n_curated_mRNAs), desc(n_correlated_pairs), .by_group = TRUE) %>%
  slice_head(n = top_lncRNAs_per_state) %>%
  ungroup()
write.csv(lnc_rank, file.path(out_dir, "Figure3_GBA_top_lncRNAs_selected.csv"), row.names = FALSE)

bubble_summary <- gba_map %>%
  semi_join(lnc_rank, by = c("correlation_state", "lncRNA")) %>%
  group_by(correlation_state, lncRNA, lnc_label, lnc_source, lnc_phase,
           pathway_id, pathway, kegg_program, correlation_sign) %>%
  summarise(
    n_mRNAs = n_distinct(mRNA_symbol_key),
    median_rho = median(rho),
    mean_abs_rho = mean(abs(rho)),
    mRNAs = paste(sort(unique(mRNA_symbol)), collapse = ";"),
    .groups = "drop"
  )
write.csv(bubble_summary, file.path(out_dir, "Figure3_GBA_bubble_summary.csv"), row.names = FALSE)

lnc_order <- lnc_rank %>%
  arrange(factor(correlation_state, levels = c("6 h", "24 h")),
          n_curated_pathways, n_curated_mRNAs) %>%
  pull(lnc_label) %>% unique()

bubble_summary <- bubble_summary %>%
  mutate(
    correlation_state = factor(correlation_state, levels = c("6 h", "24 h")),
    correlation_sign = factor(correlation_sign, levels = c("Positive", "Negative")),
    lnc_label = factor(lnc_label, levels = lnc_order),
    pathway = factor(pathway, levels = selected_kegg$pathway)
  )

p_gba <- ggplot(bubble_summary,
                aes(pathway, lnc_label, size = n_mRNAs, color = median_rho, shape = lnc_source)) +
  geom_point(alpha = 0.90, stroke = 0.8) +
  facet_grid(correlation_state ~ correlation_sign, scales = "free_y", space = "free_y") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), oob = squish,
                         name = "Median rho") +
  scale_size_continuous(range = c(2.5, 9), breaks = pretty_breaks(4), name = "Correlated mRNAs") +
  scale_shape_manual(values = source_shapes, name = "lncRNA source") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 22)) +
  labs(
    title = "Guilt-by-association of temporally regulated lncRNAs with curated KEGG programs",
    x = NULL,
    y = "Selected lncRNAs"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 8.5),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92")
  )
save_dual(p_gba, "Figure3_GBA_lncRNA_curated_KEGG_bubble_RECOMMENDED", 15, 11)

# Process-level summary: unique correlated mRNAs by sign
pathway_sign_summary <- gba_map %>%
  group_by(correlation_state, pathway_id, pathway, kegg_program, correlation_sign) %>%
  summarise(
    n_unique_lncRNAs = n_distinct(lncRNA),
    n_unique_mRNAs = n_distinct(mRNA_symbol_key),
    n_unique_pairs = n_distinct(paste(lncRNA, mRNA_id, sep = "__")),
    median_rho = median(rho),
    .groups = "drop"
  ) %>%
  mutate(
    correlation_state = factor(correlation_state, levels = c("6 h", "24 h")),
    pathway = factor(pathway, levels = rev(selected_kegg$pathway)),
    correlation_sign = factor(correlation_sign, levels = c("Positive", "Negative"))
  )
write.csv(pathway_sign_summary, file.path(out_dir, "Figure3_GBA_pathway_sign_summary.csv"), row.names = FALSE)

p_gba_pathway <- ggplot(pathway_sign_summary,
                        aes(pathway, n_unique_pairs, fill = correlation_sign)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_text(aes(label = n_unique_pairs, group = correlation_sign),
            position = position_dodge(width = 0.78),
            vjust = -0.25, size = 3.2) +
  facet_wrap(~correlation_state, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Positive" = "#B2182B", "Negative" = "#2166AC")) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 18)) +
  labs(
    title = "Curated KEGG representation among significant lncRNA–mRNA correlations",
    x = "Curated KEGG pathway",
    y = "Unique significant correlation pairs",
    fill = "Correlation"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
save_dual(p_gba_pathway, "Figure3_GBA_curated_KEGG_pathway_pair_counts", 15, 7.8)

sink(file.path(out_dir, "README_run_summary.txt"))
cat("Final missing-figure script summary\n")
cat("===================================\n\n")
cat("FIGURE 1: miRNA temporal heatmaps\n")
cat("- TOP version: up to ", top_miRNAs_per_phase, " miRNAs per temporal program\n", sep = "")
cat("- Values are temporal log2FC values at 6 h and 24 h.\n")
cat("- Early and Late non-applicable cells are displayed in grey.\n\n")
cat("FIGURE 2: global DE RNA overview\n")
cat("- Panel A: unique DE features counted once across the experiment\n")
cat("- Panel B: Early/Sustained/Late composition\n")
cat("- Panel C: 6 h and 24 h state-specific counts\n")
cat("- circRNAs are nominally altered candidate loci and remain exploratory.\n")
cat("- Direction-switching mRNAs/lncRNAs are exported separately and omitted from UP/DOWN bars.\n\n")
cat("FIGURE 3: lncRNA guilt-by-association\n")
cat("- Correlation threshold: |rho| >= ", rho_cutoff, ", p < ", correlation_p_cutoff, "\n", sep = "")
cat("- KEGG terms: exact 11 curated pathways used in the TF/pathway figure\n")
cat("- Temporal mode: ", gba_temporal_mode, "\n", sep = "")
cat("- Requires correlated mRNA to be DE at the matching state: ", gba_require_mRNA_DE, "\n", sep = "")
cat("- The available correlation analysis contains Early and Late networks; sustained correlations = 0 in the supplied summary.\n")
cat("- Correlation/pathway overlap is hypothesis-generating and does not establish lncRNA function or causality.\n\n")
cat("Input files used:\n")
cat("miRNA temporal: ", mir_temporal_file, "\n", sep = "")
cat("miRNA 6 h: ", mir6_file, "\n", sep = "")
cat("miRNA 24 h: ", mir24_file, "\n", sep = "")
cat("mRNA 6 h: ", mrna6_file, "\n", sep = "")
cat("mRNA 24 h: ", mrna24_file, "\n", sep = "")
cat("lncRNA phase: ", lnc_phase_file, "\n", sep = "")
cat("correlations 6 h: ", cor6_file, "\n", sep = "")
cat("correlations 24 h: ", cor24_file, "\n", sep = "")
sink()

message("Done. Outputs written to: ", normalizePath(out_dir))
