#!/usr/bin/env Rscript
# =============================================================================
# R18 — Figure S3C: directional overlap matrices and annotated vs novel stacked barplot
#
# PORTED, NOT REWRITTEN.
# Source: archive/.../39e_lncRNA_directional_overlap_barplot_UP_top_DOWN_bottom_FINAL.R (382 lines)
# Produced by scripts/dev/port_from_archive.sh, which applies exactly five
# classes of change and passes everything else through unmodified:
#   1. this preamble (optparse CLI, lib/ sourcing, base_dir/out_dir/dpi)
#   2. hard-coded C:/Users/sebau/... base_dir and out_dir commented out
#   3. base = "C:/Users/sebau/..." default arguments -> base = base_dir
#   4. browser-download duplicates "name(1).csv" removed from candidate lists
#   5. the /mnt/data container search path removed
#
# Diff this file against the archived original to audit the port:
#   diff archive/.../39e_lncRNA_directional_overlap_barplot_UP_top_DOWN_bottom_FINAL.R scripts/R/17_directional_overlap.R
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
    make_option("--out-dir", type = "character", default = "results/R18_directional_overlap",
                help = "output directory [default %default]"),
    make_option("--dpi", type = "integer", default = 600,
                help = "raster figure resolution [default %default]")),
  description = "\nR17: Figure S3C: directional overlap matrices and annotated vs novel stacked barplot\nPorted from 39e_lncRNA_directional_overlap_barplot_UP_top_DOWN_bottom_FINAL.R."))

base_dir <- normalizePath(.opt$`in-dir`, winslash = "/", mustWork = FALSE)
out_dir  <- ensure_out_dir(.opt$`out-dir`)
.suz_out_root <- out_dir   ## [ported] preserved: the body below rebinds out_dir
png_dpi  <- .opt$dpi
dpi_png  <- .opt$dpi
message("[R18] in-dir : ", base_dir)
message("[R18] out-dir: ", out_dir)

# ----------------------------- ported body ----------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

# ============================================================
# lncRNA supplementary figures — fixed identifier handling
# 1) Directional overlap heatmap (6h vs 24h)
#    - all lncRNAs
#    - annotated only
#    - novel only
#    - optional faceted figure
# 2) Annotated vs novel stacked barplot with counts
#    - Early / Sustained / Late
#    - DOWN at bottom, UP at top
#    - white numbers inside bars
# ============================================================

# -----------------------------
# Helpers
# -----------------------------
find_input <- function(preferred = character(), pattern = NULL, search_dir = ".") {
  candidates <- unique(preferred[nzchar(preferred)])
  for (f in candidates) {
    if (file.exists(file.path(search_dir, f))) return(file.path(search_dir, f))
    if (file.exists(f)) return(f)
  }
  if (!is.null(pattern)) {
    hits <- list.files(search_dir, pattern = pattern, full.names = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  stop("Could not find input file. Preferred: ", paste(preferred, collapse = ", "),
       if (!is.null(pattern)) paste0(" | pattern: ", pattern) else "")
}

save_both <- function(plot_obj, out_stub, width = 8, height = 6, dpi = 400) {
  target_dir <- dirname(out_stub)
  if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(target_dir)) stop("Could not create figure directory: ", target_dir)

  pdf_file <- paste0(out_stub, ".pdf")
  png_file <- paste0(out_stub, ".png")

  ggsave(filename = pdf_file, plot = plot_obj, width = width, height = height,
         units = "in", device = "pdf", bg = "white")
  ggsave(filename = png_file, plot = plot_obj, width = width, height = height,
         units = "in", dpi = dpi, bg = "white")
}

mode_or_first <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

# -----------------------------
# Load input
# -----------------------------
lnc_phase_file <- find_input(
  preferred = c(
    ## [ported] the "(1)" and "(2)" browser-download duplicates were removed here.
    "SUPP_Fig4B2_lncDE_with_phase.csv"
  ),
  pattern = "^SUPP_Fig4B2_lncDE_with_phase.*\\.csv$"
)

out_dir <- ensure_out_dir(file.path(.suz_out_root, "lnc_overlap_out"))  ## [ported] was a CWD-relative path
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create output directory: ", normalizePath(getwd(), winslash = "/", mustWork = FALSE), "/", out_dir)
if (file.access(out_dir, 2) != 0) stop("Output directory is not writable: ", out_dir)

message("Using file: ", lnc_phase_file)

# Read first, then resolve optional identifier columns safely.
# Referencing .data$gene_id inside mutate() fails when that column does not exist,
# even when gene_key is already present.
lnc_raw <- read_csv(lnc_phase_file, show_col_types = FALSE) %>%
  rename_with(~gsub("\\.+", "_", .x))

required_cols <- c("time", "direction", "phase", "log2FC")
missing_required <- setdiff(required_cols, names(lnc_raw))
if (length(missing_required) > 0) {
  stop("The lncRNA input is missing required columns: ",
       paste(missing_required, collapse = ", "))
}

# Pick the first identifier column that actually exists.
id_candidates <- c("gene_key", "gene_id", "gene_name")
id_col <- id_candidates[id_candidates %in% names(lnc_raw)][1]
if (is.na(id_col)) {
  stop("No lncRNA identifier column was found. Expected one of: ",
       paste(id_candidates, collapse = ", "))
}

gene_key_vec <- as.character(lnc_raw[[id_col]])
gene_name_vec <- if ("gene_name" %in% names(lnc_raw)) {
  as.character(lnc_raw$gene_name)
} else if ("gene_name_raw" %in% names(lnc_raw)) {
  as.character(lnc_raw$gene_name_raw)
} else {
  gene_key_vec
}
gene_name_vec[is.na(gene_name_vec) | gene_name_vec == ""] <-
  gene_key_vec[is.na(gene_name_vec) | gene_name_vec == ""]

source_vec <- if ("lnc_source" %in% names(lnc_raw)) {
  as.character(lnc_raw$lnc_source)
} else {
  rep(NA_character_, nrow(lnc_raw))
}

de_sig_vec <- if ("de_sig" %in% names(lnc_raw)) {
  x <- lnc_raw$de_sig
  if (is.logical(x)) {
    replace_na(x, FALSE)
  } else {
    toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
  }
} else {
  rep(TRUE, nrow(lnc_raw))
}

padj_vec <- if ("padj" %in% names(lnc_raw)) {
  suppressWarnings(as.numeric(lnc_raw$padj))
} else {
  rep(NA_real_, nrow(lnc_raw))
}

lnc_long <- lnc_raw %>%
  mutate(
    gene_key = gene_key_vec,
    gene_name = gene_name_vec,
    time = as.character(time),
    time = case_when(
      time %in% c("6h", "6", "6 h") ~ "6h",
      time %in% c("24h", "24", "24 h") ~ "24h",
      TRUE ~ time
    ),
    direction = toupper(as.character(direction)),
    direction = case_when(
      direction %in% c("UP", "DOWN") ~ direction,
      TRUE ~ NA_character_
    ),
    de_sig = de_sig_vec,
    padj = padj_vec,
    log2FC = suppressWarnings(as.numeric(log2FC)),
    lnc_source_original = source_vec,
    lnc_source = case_when(
      str_detect(toupper(gene_key), "^MSTRG\\.?") |
        str_detect(toupper(gene_name), "^MSTRG\\.?") ~ "novel",
      !is.na(lnc_source_original) &
        str_to_lower(lnc_source_original) %in% c("novel", "annotated") ~
        str_to_lower(lnc_source_original),
      TRUE ~ "annotated"
    ),
    phase = case_when(
      str_to_lower(as.character(phase)) == "early" ~ "Early",
      str_to_lower(as.character(phase)) == "sustained" ~ "Sustained",
      str_to_lower(as.character(phase)) == "late" ~ "Late",
      TRUE ~ as.character(phase)
    )
  ) %>%
  filter(time %in% c("6h", "24h"), !is.na(gene_key), gene_key != "")

write_csv(
  lnc_long %>% count(lnc_source_original, lnc_source, name = "n_rows"),
  file.path(out_dir, "00A_lncRNA_source_reclassification_audit.csv")
)

# -----------------------------
# Build gene-level table
# -----------------------------
# Direction per timepoint; missing timepoint becomes Not DE
state_wide <- lnc_long %>%
  mutate(state = ifelse(!is.na(direction) & de_sig, direction, "Not DE")) %>%
  select(gene_key, gene_name, lnc_source, phase, time, log2FC, state) %>%
  group_by(gene_key, time) %>%
  summarise(
    gene_name = mode_or_first(gene_name),
    lnc_source = mode_or_first(lnc_source),
    phase = mode_or_first(phase),
    log2FC = if (all(is.na(log2FC))) NA_real_ else log2FC[which.max(abs(replace_na(log2FC, 0)))],
    state = mode_or_first(state),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = gene_key,
    names_from = time,
    values_from = c(gene_name, lnc_source, phase, log2FC, state),
    names_sep = "_"
  ) %>%
  mutate(
    gene_name = coalesce(gene_name_6h, gene_name_24h, gene_key),
    lnc_source = coalesce(lnc_source_6h, lnc_source_24h,
                          ifelse(str_detect(gene_key, "^MSTRG") | str_detect(gene_name, "^MSTRG"), "novel", "annotated")),
    phase = coalesce(phase_6h, phase_24h),
    log2FC_6h = log2FC_6h,
    log2FC_24h = log2FC_24h,
    state_6h = coalesce(state_6h, "Not DE"),
    state_24h = coalesce(state_24h, "Not DE")
  ) %>%
  select(gene_key, gene_name, lnc_source, phase, log2FC_6h, log2FC_24h, state_6h, state_24h) %>%
  distinct()

# Override source one more time from MSTRG if needed
state_wide <- state_wide %>%
  mutate(lnc_source = ifelse(str_detect(gene_key, "^MSTRG") | str_detect(gene_name, "^MSTRG"), "novel", lnc_source))

write_csv(state_wide, file.path(out_dir, "00_lncRNA_gene_level_states.csv"))

# -----------------------------
# 1) Directional overlap heatmaps
# -----------------------------
make_overlap_counts <- function(df) {
  x_levels <- c("UP", "DOWN", "Not DE")
  y_levels <- c("Not DE", "DOWN", "UP")
  df %>%
    mutate(
      state_6h = factor(state_6h, levels = y_levels),
      state_24h = factor(state_24h, levels = x_levels)
    ) %>%
    count(state_6h, state_24h, name = "n") %>%
    complete(state_6h = factor(y_levels, levels = y_levels),
             state_24h = factor(x_levels, levels = x_levels),
             fill = list(n = 0)) %>%
    mutate(fill_n = ifelse(n == 0, NA_real_, n))
}

plot_overlap_heatmap <- function(count_df, title_text) {
  ggplot(count_df, aes(x = state_24h, y = state_6h, fill = fill_n)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(data = subset(count_df, n > 0), aes(label = n), size = 8, fontface = "bold", color = "black") +
    scale_fill_gradient(low = "#f3d6df", high = "#c51b8a", na.value = "grey85") +
    labs(
      title = title_text,
      x = "24 h state",
      y = "6 h state",
      fill = "lncRNAs"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "#4d4d4d"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 28),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(1.2, "cm")
    )
}

all_counts <- make_overlap_counts(state_wide)
ann_counts <- make_overlap_counts(filter(state_wide, lnc_source == "annotated"))
nov_counts <- make_overlap_counts(filter(state_wide, lnc_source == "novel"))

write_csv(all_counts, file.path(out_dir, "01_overlap_counts_all_lncRNAs.csv"))
write_csv(ann_counts, file.path(out_dir, "01_overlap_counts_annotated_lncRNAs.csv"))
write_csv(nov_counts, file.path(out_dir, "01_overlap_counts_novel_lncRNAs.csv"))

p_all <- plot_overlap_heatmap(all_counts, "Directional overlap of lncRNA DE between 6 h and 24 h")
p_ann <- plot_overlap_heatmap(ann_counts, "Directional overlap of annotated lncRNA DE between 6 h and 24 h")
p_nov <- plot_overlap_heatmap(nov_counts, "Directional overlap of novel lncRNA DE between 6 h and 24 h")

save_both(p_all, file.path(out_dir, "overlap_all"), width = 12, height = 10)
save_both(p_ann, file.path(out_dir, "overlap_annotated"), width = 12, height = 10)
save_both(p_nov, file.path(out_dir, "overlap_novel"), width = 12, height = 10)

# Optional faceted version
facet_counts <- bind_rows(
  mutate(ann_counts, lnc_source = "annotated"),
  mutate(nov_counts, lnc_source = "novel")
)

p_facet <- ggplot(facet_counts, aes(x = state_24h, y = state_6h, fill = fill_n)) +
  geom_tile(color = "white", linewidth = 1.1) +
  geom_text(data = subset(facet_counts, n > 0), aes(label = n), size = 6.5, fontface = "bold", color = "black") +
  scale_fill_gradient(low = "#f3d6df", high = "#c51b8a", na.value = "grey85") +
  facet_wrap(~lnc_source, nrow = 1) +
  labs(
    title = "Directional overlap of lncRNA DE between 6 h and 24 h",
    x = "24 h state",
    y = "6 h state",
    fill = "lncRNAs"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 22),
    axis.text = element_text(color = "#4d4d4d"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 26)
  )

save_both(p_facet, file.path(out_dir, "overlap_annotated_vs_NOVEL_faceted"), width = 18, height = 9)

# -----------------------------
# 2) Annotated vs novel stacked barplot (DOWN bottom, UP top)
# -----------------------------
bar_df <- state_wide %>%
  filter(phase %in% c("Early", "Sustained", "Late")) %>%
  mutate(direction_group = case_when(
    phase == "Early" & state_6h %in% c("UP", "DOWN") ~ state_6h,
    phase == "Late" & state_24h %in% c("UP", "DOWN") ~ state_24h,
    phase == "Sustained" & state_6h == state_24h & state_6h %in% c("UP", "DOWN") ~ state_6h,
    TRUE ~ "MIXED"
  ))

write_csv(bar_df, file.path(out_dir, "02_lncRNA_barplot_gene_assignments.csv"))
write_csv(filter(bar_df, direction_group == "MIXED"), file.path(out_dir, "02B_lncRNA_direction_switching_or_unassigned.csv"))

bar_counts <- bar_df %>%
  filter(direction_group %in% c("UP", "DOWN")) %>%
  mutate(
    phase = factor(phase, levels = c("Early", "Sustained", "Late")),
    lnc_source = factor(lnc_source, levels = c("annotated", "novel")),
    # For ggplot stacked bars, using levels = c("UP", "DOWN") places
    # DOWN at the bottom and UP at the top in the final stack.
    direction_stack = factor(direction_group, levels = c("UP", "DOWN"))
  ) %>%
  count(lnc_source, phase, direction_stack, name = "n")

write_csv(bar_counts, file.path(out_dir, "03_lncRNA_annotated_novel_barplot_counts.csv"))

p_bar <- ggplot(bar_counts, aes(x = phase, y = n, fill = direction_stack)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            color = "white", size = 6, fontface = "bold") +
  facet_wrap(~lnc_source, nrow = 1) +
  scale_fill_manual(values = c("DOWN" = "#4C78A8", "UP" = "#E07B39"),
                    breaks = c("UP", "DOWN"),
                    labels = c("UP", "DOWN"),
                    drop = FALSE) +
  labs(
    title = "lncRNA temporal programs: annotated vs novel",
    x = "Response class",
    y = "# lncRNA genes",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 22),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 28),
    legend.title = element_text(face = "bold")
  )

save_both(p_bar, file.path(out_dir, "bar_annotated_novel"), width = 14, height = 8)

# -----------------------------
# Session info / manifest
# -----------------------------
manifest <- tibble::tibble(
  file = c(
    basename(lnc_phase_file),
    "00_lncRNA_gene_level_states.csv",
    "01_overlap_counts_all_lncRNAs.csv",
    "01_overlap_counts_annotated_lncRNAs.csv",
    "01_overlap_counts_novel_lncRNAs.csv",
    "02_lncRNA_barplot_gene_assignments.csv",
    "02B_lncRNA_direction_switching_or_unassigned.csv",
    "03_lncRNA_annotated_novel_barplot_counts.csv"
  ),
  description = c(
    "Input lncRNA long-format DE/phase table",
    "Gene-level wide table used for all outputs",
    "Directional overlap counts for all lncRNAs",
    "Directional overlap counts for annotated lncRNAs",
    "Directional overlap counts for novel lncRNAs",
    "Per-gene assignments used for the annotated-vs-novel barplot",
    "Genes excluded from the barplot because sustained direction was mixed or missing",
    "Final counts shown in the stacked barplot"
  )
)
write_csv(manifest, file.path(out_dir, "00_manifest.csv"))

message("Done. Outputs written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
