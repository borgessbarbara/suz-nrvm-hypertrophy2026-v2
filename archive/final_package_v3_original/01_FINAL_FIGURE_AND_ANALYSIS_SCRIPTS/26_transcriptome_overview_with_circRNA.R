#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# =========================
# User-adjustable paths
# =========================
annot_file   <- "hisat2_6_24_meta_ann.csv"
novel_file   <- "hisat2_6_24_meta_novel_lncRNA.csv"
circ6_file   <- "raw_count_matrix_all_6.txt"
circ24_file  <- "raw_count_matrix_all_24.txt"
out_dir      <- "Transcriptome_overview_with_circRNA"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Helper functions
# =========================
normalize_biotype <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    tolower()
}

# Format numbers with commas for text labels
fmt_n <- function(x) scales::comma(x)

save_dual <- function(plot_obj, filename_base, out_dir, w = 16, h = 7.5, dpi = 600) {
  ggsave(file.path(out_dir, paste0(filename_base, ".pdf")), plot_obj, width = w, height = h,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(filename_base, ".png")), plot_obj, width = w, height = h,
         units = "in", dpi = dpi, bg = "white")
}

# =========================
# 1) Annotated gene biotypes
# =========================
ann_all <- read.csv(annot_file, stringsAsFactors = FALSE)

ann_counts <- ann_all %>%
  transmute(
    gene_id = ensembl_gene_id,
    gene_biotype = normalize_biotype(gene_biotype)
  ) %>%
  distinct(gene_id, gene_biotype) %>%
  filter(!is.na(gene_id), gene_id != "", !is.na(gene_biotype), gene_biotype != "") %>%
  count(gene_biotype, name = "count")

# Separate annotated lncRNA explicitly
annot_lnc_count <- ann_counts %>%
  filter(gene_biotype == "lncrna") %>%
  pull(count)
if (length(annot_lnc_count) == 0) annot_lnc_count <- 0

ann_counts_main <- ann_counts %>%
  filter(gene_biotype != "lncrna")

# =========================
# 2) Novel lncRNA genes
# =========================
ann_novel_raw <- read.csv(novel_file, stringsAsFactors = FALSE)

novel_lnc_count <- ann_novel_raw %>%
  filter(type == "transcript") %>%
  transmute(gene_id = gene_id, pred_biotype = normalize_biotype(pred_biotype)) %>%
  filter(pred_biotype == "lncrna") %>%
  distinct(gene_id) %>%
  nrow()

# =========================
# 3) circRNA count from union across 6 h and 24 h
# =========================
read_circ_ids <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  first_col <- names(df)[1]
  unique(as.character(df[[first_col]]))
}

circ6_ids <- read_circ_ids(circ6_file)
circ24_ids <- read_circ_ids(circ24_file)
all_circ_ids <- union(circ6_ids, circ24_ids)
circ_count <- length(all_circ_ids)

# =========================
# 4) Build summary table
# =========================
summary_tbl <- bind_rows(
  tibble(gene_biotype = "protein_coding", count = ann_counts_main %>% filter(gene_biotype == "protein_coding") %>% pull(count) %>% {if(length(.)==0) 0 else .}),
  tibble(gene_biotype = "annotated_lncrna", count = annot_lnc_count),
  tibble(gene_biotype = "novel_lncrna", count = novel_lnc_count),
  tibble(gene_biotype = "circrna", count = circ_count),
  ann_counts_main %>% filter(!gene_biotype %in% "protein_coding")
) %>%
  group_by(gene_biotype) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(
    label = case_when(
      gene_biotype == "annotated_lncrna" ~ "annotated lncRNA",
      gene_biotype == "novel_lncrna" ~ "novel lncRNA",
      gene_biotype == "circrna" ~ "circRNA",
      TRUE ~ gene_biotype
    )
  ) %>%
  arrange(desc(count)) %>%
  mutate(label = factor(label, levels = label))

# Save source table for manuscript / supplement
write.csv(summary_tbl %>% select(label, count),
          file.path(out_dir, "transcriptome_overview_with_circRNA_counts.csv"),
          row.names = FALSE)

# =========================
# 5) Plot
# =========================
p <- ggplot(summary_tbl, aes(x = label, y = count, fill = label)) +
  geom_col(width = 0.92, color = NA) +
  geom_text(aes(label = fmt_n(count)), vjust = -0.35, size = 4.5) +
  scale_y_log10(
    labels = label_comma(),
    breaks = c(1, 10, 100, 1000, 10000, 100000),
    expand = expansion(mult = c(0.02, 0.14))
  ) +
  scale_fill_manual(values = grDevices::hcl.colors(nrow(summary_tbl), palette = "Spectral", rev = TRUE)) +
  labs(
    title = "Transcriptome Overview 6 and 24 h",
    x = "Gene Biotypes",
    y = "Count (log10 scale)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

save_dual(p,
          filename_base = "Transcriptome_Overview_6h_24h_with_circRNA",
          out_dir = out_dir,
          w = 16, h = 7.5, dpi = 600)

message("Done. Files written to: ", out_dir)
message("circRNA unique loci across 6 h and 24 h = ", circ_count)

