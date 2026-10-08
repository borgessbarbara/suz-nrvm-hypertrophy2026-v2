#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

# -----------------------------
# Helper functions
# -----------------------------
find_input <- function(candidates = NULL, pattern = NULL) {
  all_files <- list.files('.', recursive = TRUE, full.names = TRUE)
  all_files <- unique(c(all_files, candidates[file.exists(candidates)]))

  if (!is.null(candidates)) {
    for (cand in candidates) {
      direct_hit <- candidates[file.exists(candidates) & basename(candidates) == basename(cand)]
      if (length(direct_hit) > 0) return(direct_hit[1])

      hit <- all_files[basename(all_files) == basename(cand)]
      if (length(hit) > 0) return(hit[1])
    }
  }

  if (!is.null(pattern)) {
    hit <- all_files[grepl(pattern, basename(all_files), ignore.case = TRUE)]
    if (length(hit) > 0) {
      info <- file.info(hit)
      return(hit[order(info$mtime, decreasing = TRUE)][1])
    }
  }

  stop('Could not find the required input file in the working directory or its subfolders.')
}

safe_label <- function(gene_name, gene_key) {
  gene_name <- ifelse(is.na(gene_name) | gene_name == '' | gene_name == 'NA', gene_key, gene_name)
  gene_name
}

# MSTRG/StringTie identifiers always override a mislabeled input lnc_source field.
classify_lnc_source <- function(gene_key, gene_name = NULL, gene_name_raw = NULL, source_original = NULL) {
  key <- toupper(trimws(as.character(gene_key)))
  gname <- if (is.null(gene_name)) rep('', length(key)) else toupper(trimws(as.character(gene_name)))
  graw <- if (is.null(gene_name_raw)) rep('', length(key)) else toupper(trimws(as.character(gene_name_raw)))
  src <- if (is.null(source_original)) rep('', length(key)) else tolower(trimws(as.character(source_original)))

  is_mstrg <- str_detect(key, '(^|[|:_-])MSTRG\\.?[0-9]') |
              str_detect(gname, '(^|[|:_-])MSTRG\\.?[0-9]') |
              str_detect(graw, '(^|[|:_-])MSTRG\\.?[0-9]')

  case_when(
    is_mstrg ~ 'novel',
    src %in% c('novel', 'stringtie', 'reconstructed') ~ 'novel',
    TRUE ~ 'annotated'
  )
}

calc_row_gaps <- function(df, gap_cols) {
  grp <- apply(df[, gap_cols, drop = FALSE], 1, paste, collapse = ' | ')
  idx <- which(grp[-1] != grp[-length(grp)])
  idx
}

make_unique_rownames <- function(x) make.unique(as.character(x))

build_lnc_wide <- function(df) {
  req <- c('gene_key','gene_name','log2FC','time','phase')
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop('Missing columns: ', paste(miss, collapse = ', '))

  if ('de_sig' %in% colnames(df)) {
    df <- df %>% filter(is.na(de_sig) | de_sig)
  }

  source_original <- if ('lnc_source' %in% colnames(df)) df$lnc_source else rep(NA_character_, nrow(df))
  gene_name_raw <- if ('gene_name_raw' %in% colnames(df)) df$gene_name_raw else rep(NA_character_, nrow(df))

  df <- df %>%
    mutate(
      time = as.character(time),
      phase = factor(as.character(phase), levels = c('Early', 'Sustained', 'Late')),
      lnc_source_original = source_original,
      lnc_source = classify_lnc_source(gene_key, gene_name, gene_name_raw, source_original),
      lnc_source = factor(lnc_source, levels = c('annotated', 'novel')),
      display_id = safe_label(gene_name, gene_key)
    )

  wide <- df %>%
    select(gene_key, display_id, phase, lnc_source, time, log2FC) %>%
    distinct() %>%
    pivot_wider(names_from = time, values_from = log2FC, names_prefix = 'log2FC_') %>%
    mutate(
      abs_max = pmax(abs(log2FC_6h), abs(log2FC_24h), na.rm = TRUE)
    ) %>%
    arrange(phase, lnc_source, desc(abs_max), display_id)

  wide
}

plot_heatmap <- function(plot_df, out_prefix, title_text, show_names = TRUE, width = 8.5, height = 12) {
  if (nrow(plot_df) == 0) return(invisible(NULL))

  mat <- plot_df %>%
    select(log2FC_6h, log2FC_24h) %>%
    as.data.frame()
  rownames(mat) <- make_unique_rownames(plot_df$display_id)
  mat <- as.matrix(mat)

  ann_row <- plot_df %>%
    select(lnc_source, phase) %>%
    as.data.frame()
  rownames(ann_row) <- rownames(mat)

  ann_colors <- list(
    lnc_source = c(annotated = '#C77CDE', novel = '#00BFC4'),
    phase      = c(Early = '#7CAE00', Sustained = '#F8766D', Late = '#00BFC4')
  )

  gaps <- calc_row_gaps(
    plot_df %>% mutate(phase = as.character(phase), lnc_source = as.character(lnc_source)),
    c('phase', 'lnc_source')
  )

  color_fun <- colorRampPalette(c('#2166AC', '#F7F7F7', '#B2182B'))
  hm_colors <- color_fun(100)

  max_abs <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  breaks <- seq(-max_abs, max_abs, length.out = 101)

  p <- pheatmap(
    mat,
    color = hm_colors,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = show_names,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 11,
    border_color = 'grey75',
    main = title_text,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    gaps_row = gaps,
    na_col = 'grey80',
    cellwidth = 90,
    cellheight = ifelse(nrow(plot_df) > 40, 12, 18)
  )

  pdf(paste0(out_prefix, '.pdf'), width = width, height = height, useDingbats = FALSE)
  grid.newpage(); grid.draw(p$gtable)
  dev.off()

  png(paste0(out_prefix, '.png'), width = width, height = height, units = 'in', res = 400)
  grid.newpage(); grid.draw(p$gtable)
  dev.off()
}

select_top_by_phase_source <- function(df, top_n = 10) {
  df %>%
    group_by(phase, lnc_source) %>%
    slice_max(order_by = abs_max, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(phase, lnc_source, desc(abs_max), display_id)
}

# -----------------------------
# Inputs / outputs
# -----------------------------
lnc_file <- find_input(
  candidates = c('SUPP_Fig4B2_lncDE_with_phase.csv',
                 'SUPP_Fig4B2_lncDE_with_phase(1).csv',
                 'SUPP_Fig4B2_lncDE_with_phase(2).csv'),
  pattern = '^SUPP_Fig4B2_lncDE_with_phase.*\\.csv$'
)

out_dir <- file.path(getwd(), 'Figure_lncRNA_temporal_heatmaps_MSTRG_FIXED')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message('Using input: ', lnc_file)
message('Output dir: ', out_dir)

# -----------------------------
# Build plotting table
# -----------------------------
lnc_raw  <- read_csv(lnc_file, show_col_types = FALSE)
lnc_wide <- build_lnc_wide(lnc_raw)
write_csv(lnc_wide, file.path(out_dir, '00_lncRNA_temporal_heatmap_input_table.csv'))

source_diagnostic <- lnc_raw %>%
  mutate(
    source_original = if ('lnc_source' %in% colnames(lnc_raw)) as.character(lnc_source) else NA_character_,
    source_corrected = classify_lnc_source(
      gene_key,
      if ('gene_name' %in% colnames(lnc_raw)) gene_name else NA_character_,
      if ('gene_name_raw' %in% colnames(lnc_raw)) gene_name_raw else NA_character_,
      source_original
    ),
    is_MSTRG = str_detect(toupper(as.character(gene_key)), 'MSTRG')
  ) %>%
  count(source_original, source_corrected, is_MSTRG, name = 'n_rows')

write_csv(source_diagnostic, file.path(out_dir, '00A_lncRNA_source_reclassification_diagnostic.csv'))

source_by_phase <- lnc_wide %>%
  count(phase, lnc_source, name = 'n_unique_lncRNAs') %>%
  arrange(phase, lnc_source)
write_csv(source_by_phase, file.path(out_dir, '00B_lncRNA_source_counts_by_phase_FIXED.csv'))

# -----------------------------
# Plot outputs
# -----------------------------
# 1) All lncRNAs, combined, with annotated vs novel shown in row annotation
plot_heatmap(
  lnc_wide,
  out_prefix = file.path(out_dir, 'Figure_lncRNA_temporal_heatmap_ALL_annotated_vs_novel'),
  title_text = 'lncRNAs by response class (annotated vs novel)',
  show_names = TRUE,
  width = 8.5,
  height = max(9, min(18, 0.25 * nrow(lnc_wide) + 3))
)

# 2) Top lncRNAs per phase and source for a cleaner main-figure option
lnc_top <- select_top_by_phase_source(lnc_wide, top_n = 10)
write_csv(lnc_top, file.path(out_dir, '00B_lncRNA_temporal_heatmap_top_table.csv'))
plot_heatmap(
  lnc_top,
  out_prefix = file.path(out_dir, 'Figure_lncRNA_temporal_heatmap_TOP_annotated_vs_novel'),
  title_text = 'Top lncRNAs by response class (annotated vs novel)',
  show_names = TRUE,
  width = 8.5,
  height = max(8, min(14, 0.30 * nrow(lnc_top) + 3))
)

# 3) Annotated-only panel
lnc_ann <- lnc_wide %>% filter(lnc_source == 'annotated')
plot_heatmap(
  lnc_ann,
  out_prefix = file.path(out_dir, 'Figure_lncRNA_temporal_heatmap_ANNOTATED_only'),
  title_text = 'Annotated lncRNAs by response class',
  show_names = TRUE,
  width = 8,
  height = max(8, min(16, 0.22 * nrow(lnc_ann) + 3))
)

# 4) Novel-only panel
lnc_nov <- lnc_wide %>% filter(lnc_source == 'novel')
plot_heatmap(
  lnc_nov,
  out_prefix = file.path(out_dir, 'Figure_lncRNA_temporal_heatmap_NOVEL_only'),
  title_text = 'Novel lncRNAs by response class',
  show_names = TRUE,
  width = 8,
  height = max(6, min(14, 0.25 * nrow(lnc_nov) + 3))
)

message('Corrected source counts by phase:')
print(source_by_phase)
message('Done. Generated heatmaps in: ', out_dir)
