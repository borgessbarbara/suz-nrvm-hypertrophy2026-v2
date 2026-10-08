#!/usr/bin/env Rscript
# =============================================================================
# R20 - Figure S3A-B circRNA: exploratory volcano plots on a nominal p threshold.
#
# Usage:   Rscript scripts/R/20_circRNA_volcanoes.R [--in-dir DIR] [--out-dir DIR] [--dpi N]
#          Rscript scripts/R/20_circRNA_volcanoes.R --help
#
# Adapted from archive/final_package_v3_original/01_FINAL_FIGURE_AND_ANALYSIS_SCRIPTS/25_recreate_circRNA_volcanoes_with_all_NS.R
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
    make_option("--out-dir", type = "character", default = "results/R20_circRNA_volcanoes",
                help = "output directory [default %default]"),
    make_option("--dpi", type = "integer", default = 600,
                help = "raster figure resolution [default %default]")),
  description = "\nR19: Figure S3A-B circRNA: exploratory nominal-p volcanoes\nPorted from 25_recreate_circRNA_volcanoes_with_all_NS.R."))

base_dir <- normalizePath(.opt$`in-dir`, winslash = "/", mustWork = FALSE)
out_dir  <- ensure_out_dir(.opt$`out-dir`)
png_dpi  <- .opt$dpi
dpi_png  <- .opt$dpi
message("[R20] in-dir : ", base_dir)
message("[R20] out-dir: ", out_dir)


suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})


# base_dir <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis"

# out_dir <- file.path(base_dir, "CircRNA_full_volcanoes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw6_name  <- "raw_count_matrix_all_6.txt"
raw24_name <- "raw_count_matrix_all_24.txt"
annotation6_name  <- "06_circRNA_sequences_de6_DE_complete_table.txt"
annotation24_name <- "07_circRNA_sequences_de24_DE_complete_table.txt"

lfc_cut <- 0.5
alpha_nominal <- 0.05
alpha_fdr <- 0.05
n_labels_per_direction <- 10
png_dpi <- 600

# The supplied count matrices already represent the circRNA inventory used downstream.
# Keep all rows by default so the volcano background includes non-significant circRNAs.
# Change this to TRUE only if an additional count filter is explicitly required.
apply_extra_count_filter <- FALSE
minimum_count <- 2
minimum_samples <- 1

find_exact_file <- function(filename, base = ".", required = TRUE) {
  hits <- list.files(base, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == filename]
  if (length(hits) == 0) {
    if (required) stop("File not found: ", filename)
    return(NA_character_)
  }
  hits[1]
}

read_table_auto <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  if (length(first) == 0) stop("Empty file: ", path)
  if (str_detect(first, "\\t")) return(read_tsv(path, show_col_types = FALSE))
  if (str_detect(first, ";") && !str_detect(first, ",")) {
    return(read_delim(path, delim = ";", show_col_types = FALSE, trim_ws = TRUE))
  }
  read_csv(path, show_col_types = FALSE)
}

first_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

parse_host <- function(x) {
  x <- as.character(x)
  out <- if_else(str_detect(x, "\\|"), str_trim(str_replace(x, "^.*\\|\\s*", "")), NA_character_)
  na_if(out, "")
}

read_optional_annotation <- function(path, timepoint) {
  if (is.na(path) || !file.exists(path)) {
    return(tibble(circRNA_ID = character(), host_gene = character(), timepoint = character()))
  }
  x <- read_table_auto(path)
  id_col <- first_col(x, c("circ_id", "circRNA_ID", "circRNA_id", "ID"))
  id2_col <- first_col(x, c("circRNA_ID2", "circRNA_ID_2", "circRNA_id2"))
  if (is.na(id_col)) {
    return(tibble(circRNA_ID = character(), host_gene = character(), timepoint = character()))
  }
  tibble(
    circRNA_ID = as.character(x[[id_col]]),
    host_gene = if (!is.na(id2_col)) parse_host(x[[id2_col]]) else NA_character_,
    timepoint = timepoint
  ) %>%
    group_by(circRNA_ID, timepoint) %>%
    summarise(
      host_gene = {
        z <- unique(na.omit(host_gene))
        if (length(z) == 0) NA_character_ else paste(z, collapse = ";")
      },
      .groups = "drop"
    )
}

prepare_counts <- function(path, id_col, sample_cols) {
  x <- read_table_auto(path)
  missing_cols <- setdiff(c(id_col, sample_cols), names(x))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", basename(path), ": ", paste(missing_cols, collapse = ", "))
  }

  x <- x %>%
    transmute(
      circRNA_ID = as.character(.data[[id_col]]),
      across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))
    )

  if (anyDuplicated(x$circRNA_ID) > 0) {
    message("Duplicated circRNA IDs detected in ", basename(path), "; counts will be summed by locus.")
    x <- x %>%
      group_by(circRNA_ID) %>%
      summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  }

  counts <- as.matrix(x[, sample_cols, drop = FALSE])
  rownames(counts) <- x$circRNA_ID
  storage.mode(counts) <- "integer"

  if (any(is.na(counts))) stop("NA counts detected in ", basename(path))
  if (any(counts < 0)) stop("Negative counts detected in ", basename(path))

  if (apply_extra_count_filter) {
    keep <- rowSums(counts >= minimum_count) >= minimum_samples
    message("Additional count filter retained ", sum(keep), " of ", nrow(counts), " circRNAs in ", basename(path), ".")
    counts <- counts[keep, , drop = FALSE]
  }

  counts
}

run_pairwise_deseq <- function(counts, control_samples, ne_samples, timepoint) {
  sample_order <- c(control_samples, ne_samples)
  counts <- counts[, sample_order, drop = FALSE]

  coldata <- data.frame(
    row.names = sample_order,
    condition = factor(c(rep("Ctrl", length(control_samples)), rep("NE", length(ne_samples))), levels = c("Ctrl", "NE"))
  )

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = coldata,
    design = ~ condition
  )

  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "NE", "Ctrl"), independentFiltering = TRUE)

  as.data.frame(res) %>%
    rownames_to_column("circRNA_ID") %>%
    as_tibble() %>%
    mutate(
      timepoint = timepoint,
      baseMean = as.numeric(baseMean),
      log2FoldChange = as.numeric(log2FoldChange),
      lfcSE = as.numeric(lfcSE),
      stat = as.numeric(stat),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj)
    )
}

add_status <- function(df) {
  df %>%
    mutate(
      status_nominal = case_when(
        !is.na(pvalue) & pvalue < alpha_nominal & log2FoldChange >=  lfc_cut ~ "Up",
        !is.na(pvalue) & pvalue < alpha_nominal & log2FoldChange <= -lfc_cut ~ "Down",
        TRUE ~ "NotSig"
      ),
      status_fdr = case_when(
        !is.na(padj) & padj < alpha_fdr & log2FoldChange >=  lfc_cut ~ "Up",
        !is.na(padj) & padj < alpha_fdr & log2FoldChange <= -lfc_cut ~ "Down",
        TRUE ~ "NotSig"
      ),
      label = if_else(!is.na(host_gene), paste0(host_gene, " (", circRNA_ID, ")"), circRNA_ID)
    )
}

choose_labels <- function(df, status_col, metric_col, n_each = 10) {
  metric_sym <- rlang::sym(metric_col)
  status_sym <- rlang::sym(status_col)

  bind_rows(
    df %>%
      filter(!!status_sym == "Up", !is.na(!!metric_sym)) %>%
      arrange(!!metric_sym, desc(abs(log2FoldChange))) %>%
      slice_head(n = n_each),
    df %>%
      filter(!!status_sym == "Down", !is.na(!!metric_sym)) %>%
      arrange(!!metric_sym, desc(abs(log2FoldChange))) %>%
      slice_head(n = n_each)
  ) %>%
    distinct(circRNA_ID, .keep_all = TRUE)
}

make_volcano <- function(df, timepoint, metric = c("pvalue", "padj"), x_limit) {
  metric <- match.arg(metric)
  status_col <- if (metric == "pvalue") "status_nominal" else "status_fdr"
  alpha <- if (metric == "pvalue") alpha_nominal else alpha_fdr
  metric_title <- if (metric == "pvalue") "nominal P-value" else "adjusted P-value"
  status_sym <- rlang::sym(status_col)

  plot_df <- df %>%
    mutate(
      metric_value = .data[[metric]],
      metric_for_plot = case_when(
        is.na(metric_value) ~ NA_real_,
        metric_value <= 0 ~ .Machine$double.xmin,
        TRUE ~ metric_value
      ),
      neglog10 = -log10(metric_for_plot)
    ) %>%
    filter(!is.na(log2FoldChange), !is.na(neglog10))

  label_df <- choose_labels(plot_df, status_col, "metric_value", n_labels_per_direction)

  n_tested <- nrow(plot_df)
  n_up <- sum(plot_df[[status_col]] == "Up", na.rm = TRUE)
  n_down <- sum(plot_df[[status_col]] == "Down", na.rm = TRUE)

  ggplot(plot_df, aes(log2FoldChange, neglog10, color = !!status_sym)) +
    geom_point(size = 1.7, alpha = 0.72) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 3.0,
      max.overlaps = Inf,
      seed = 2026,
      box.padding = 0.35,
      point.padding = 0.20,
      min.segment.length = 0.05,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("Down" = "#2F6BDE", "NotSig" = "grey72", "Up" = "#C83737"),
      breaks = c("Down", "NotSig", "Up"),
      drop = FALSE
    ) +
    coord_cartesian(xlim = c(-x_limit, x_limit)) +
    labs(
      title = paste0("circRNAs at ", timepoint, ": NE vs Ctrl"),
      subtitle = paste0(
        "All circRNAs from the count matrix; ", metric_title,
        " < ", alpha, " and |log2FC| >= ", lfc_cut,
        "  |  plotted = ", n_tested,
        ", Up = ", n_up, ", Down = ", n_down
      ),
      x = expression(log[2]~"Fold Change (NE vs Ctrl)"),
      y = if (metric == "pvalue") expression(-log[10](P~value)) else expression(-log[10](adjusted~P~value)),
      color = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle = element_text(hjust = 0.5, size = 9.5),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

save_plot <- function(p, stub, width = 10, height = 8) {
  ggsave(file.path(out_dir, paste0(stub, ".pdf")), p, width = width, height = height, bg = "white")
  ggsave(file.path(out_dir, paste0(stub, ".png")), p, width = width, height = height, dpi = png_dpi, bg = "white")
}

raw6_path <- find_exact_file(raw6_name, base_dir)
raw24_path <- find_exact_file(raw24_name, base_dir)
ann6_path <- find_exact_file(annotation6_name, base_dir, required = FALSE)
ann24_path <- find_exact_file(annotation24_name, base_dir, required = FALSE)

counts6 <- prepare_counts(
  raw6_path,
  id_col = "circRNA_ID",
  sample_cols = c("c6_1", "c6_2", "c6_4", "n6_1", "n6_2", "n6_4")
)
counts24 <- prepare_counts(
  raw24_path,
  id_col = "circRNA_ID",
  sample_cols = c("c24_1", "c24_2", "c24_4", "n24_1", "n24_2", "n24_4")
)

res6 <- run_pairwise_deseq(
  counts6,
  control_samples = c("c6_1", "c6_2", "c6_4"),
  ne_samples = c("n6_1", "n6_2", "n6_4"),
  timepoint = "6 h"
)
res24 <- run_pairwise_deseq(
  counts24,
  control_samples = c("c24_1", "c24_2", "c24_4"),
  ne_samples = c("n24_1", "n24_2", "n24_4"),
  timepoint = "24 h"
)

annotation <- bind_rows(
  read_optional_annotation(ann6_path, "6 h"),
  read_optional_annotation(ann24_path, "24 h")
)

res6 <- res6 %>%
  left_join(annotation %>% filter(timepoint == "6 h") %>% select(-timepoint), by = "circRNA_ID") %>%
  add_status()
res24 <- res24 %>%
  left_join(annotation %>% filter(timepoint == "24 h") %>% select(-timepoint), by = "circRNA_ID") %>%
  add_status()

write_csv(res6, file.path(out_dir, "circRNA_full_DESeq2_results_6h_all_rows.csv"), na = "")
write_csv(res24, file.path(out_dir, "circRNA_full_DESeq2_results_24h_all_rows.csv"), na = "")
write_tsv(res6, file.path(out_dir, "circRNA_full_DESeq2_results_6h_all_rows.txt"), na = "")
write_tsv(res24, file.path(out_dir, "circRNA_full_DESeq2_results_24h_all_rows.txt"), na = "")

all_lfc <- c(res6$log2FoldChange, res24$log2FoldChange)
all_lfc <- all_lfc[is.finite(all_lfc)]
x_limit <- max(5, ceiling(max(abs(all_lfc), na.rm = TRUE)))

p6_nominal <- make_volcano(res6, "6 h", metric = "pvalue", x_limit = x_limit)
p24_nominal <- make_volcano(res24, "24 h", metric = "pvalue", x_limit = x_limit)
p6_fdr <- make_volcano(res6, "6 h", metric = "padj", x_limit = x_limit)
p24_fdr <- make_volcano(res24, "24 h", metric = "padj", x_limit = x_limit)

save_plot(p6_nominal, "SUPP_FigCircA_volcano_6h_ALL_with_NS_nominal")
save_plot(p24_nominal, "SUPP_FigCircB_volcano_24h_ALL_with_NS_nominal")
save_plot(p6_fdr, "AUDIT_volcano_6h_ALL_with_NS_FDR")
save_plot(p24_fdr, "AUDIT_volcano_24h_ALL_with_NS_FDR")

summary_table <- bind_rows(res6, res24) %>%
  group_by(timepoint) %>%
  summarise(
    rows_in_full_result = n(),
    rows_with_pvalue = sum(!is.na(pvalue)),
    rows_with_padj = sum(!is.na(padj)),
    nominal_up = sum(status_nominal == "Up", na.rm = TRUE),
    nominal_down = sum(status_nominal == "Down", na.rm = TRUE),
    fdr_up = sum(status_fdr == "Up", na.rm = TRUE),
    fdr_down = sum(status_fdr == "Down", na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_table, file.path(out_dir, "SUMMARY_circRNA_full_volcano_counts.csv"))

sink(file.path(out_dir, "README_circRNA_full_volcanoes.txt"))
cat("Full circRNA volcano reconstruction\n")
cat("===================================\n\n")
cat("The previous 06/07 circRNA sequence tables contain only nominally selected circRNAs and alternative sequence records.\n")
cat("They are not suitable as the full background for a volcano plot.\n\n")
cat("This script reruns pairwise DESeq2 from the complete 6 h and 24 h circRNA count matrices and retains all result rows.\n")
cat("Primary exploratory circRNA plots use nominal P < 0.05 and |log2FC| >= 0.5, matching the previous circRNA figures.\n")
cat("Separate FDR audit plots use padj < 0.05 and |log2FC| >= 0.5.\n\n")
print(summary_table)
cat("\nPrimary figure files:\n")
cat(" - SUPP_FigCircA_volcano_6h_ALL_with_NS_nominal.pdf/png\n")
cat(" - SUPP_FigCircB_volcano_24h_ALL_with_NS_nominal.pdf/png\n")
cat("\nFull DESeq2 tables:\n")
cat(" - circRNA_full_DESeq2_results_6h_all_rows.csv/txt\n")
cat(" - circRNA_full_DESeq2_results_24h_all_rows.csv/txt\n")
sink()

message("Done. Outputs written to: ", normalizePath(out_dir))
