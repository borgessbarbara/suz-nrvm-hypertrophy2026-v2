#!/usr/bin/env Rscript
# =============================================================================
# R06 — Figure 2B: RNA-class inventory including circRNA loci.
#
# Ported from 26_transcriptome_overview_with_circRNA.R (151 lines), which was
# the most portable script in the whole collection: it was the only one in the
# FINAL package with CWD-relative rather than absolute paths.
#
# WHAT CHANGED
#  * optparse CLI; the four "User-adjustable paths" at the top of the original
#    become arguments.
#  * The two hisat2_6_24_meta_* inputs now resolve through read_frozen(), so
#    the Ensembl-114 tables under data/frozen/ are used by default. Pass
#    --meta-ann / --meta-novel to point at freshly generated R01 output
#    instead - but see docs/OPEN_ITEMS.md before doing so for the manuscript.
#  * save_dual() -> lib/plotting.R save_figure(), with --dpi.
#  * The analysis itself (biotype normalisation, the explicit annotated vs
#    novel lncRNA split, the circRNA union across 6 h and 24 h, the log10
#    axis) is unchanged.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(stringr)
  library(tibble); library(ggplot2); library(scales); library(readr)
})

here <- dirname(normalizePath(sub("^--file=", "",
         commandArgs(trailingOnly = FALSE)[grepl("^--file=",
         commandArgs(trailingOnly = FALSE))][1]), winslash = "/"))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "plotting.R"))

opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--meta-ann", type = "character", default = NULL,
                help = "hisat2_6_24_meta_ann.csv [default: data/frozen/]"),
    make_option("--meta-novel", type = "character", default = NULL,
                help = "hisat2_6_24_meta_novel_lncRNA.csv [default: data/frozen/]"),
    make_option("--circ-6", type = "character", default = NULL,
                help = "raw_count_matrix_all_6.txt [R06]"),
    make_option("--circ-24", type = "character", default = NULL,
                help = "raw_count_matrix_all_24.txt [R06]"),
    make_option("--in-dir", type = "character", default = "results",
                help = "root to search for inputs [default %default]"),
    make_option("--out-dir", type = "character",
                default = "results/R06_transcriptome_overview",
                help = "[default %default]"),
    make_option("--dpi", type = "integer", default = 600, help = "[default %default]"),
    make_option("--skip-circrna", action = "store_true", default = FALSE,
                help = "Omit the circRNA bar if the circRNA branch has not been run.")),
  description = "\nR04: Figure 2B, RNA-class inventory across 6 h and 24 h."))

out_dir <- ensure_out_dir(opt$`out-dir`)
repo_root <- normalizePath(file.path(here, "..", ".."), winslash = "/", mustWork = FALSE)

normalize_biotype <- function(x) tolower(str_trim(as.character(x)))

# --- 1) annotated gene biotypes ----------------------------------------------
ann_all <- if (is.null(opt$`meta-ann`)) {
  message("[R06] using the frozen Ensembl-114 meta_ann table (data/frozen/)")
  read_frozen("hisat2_6_24_meta_ann.csv", repo_root)
} else read_auto(opt$`meta-ann`)

ann_counts <- ann_all %>%
  transmute(gene_id = ensembl_gene_id,
            gene_biotype = normalize_biotype(gene_biotype)) %>%
  distinct(gene_id, gene_biotype) %>%
  filter(!is.na(gene_id), gene_id != "",
         !is.na(gene_biotype), gene_biotype != "") %>%
  count(gene_biotype, name = "count")

annot_lnc_count <- ann_counts %>% filter(gene_biotype == "lncrna") %>% pull(count)
if (length(annot_lnc_count) == 0) annot_lnc_count <- 0
ann_counts_main <- ann_counts %>% filter(gene_biotype != "lncrna")

# --- 2) novel lncRNA genes ----------------------------------------------------
ann_novel_raw <- if (is.null(opt$`meta-novel`)) {
  read_frozen("hisat2_6_24_meta_novel_lncRNA.csv", repo_root)
} else read_auto(opt$`meta-novel`)

novel_lnc_count <- ann_novel_raw %>%
  filter(type == "transcript") %>%
  transmute(gene_id = gene_id, pred_biotype = normalize_biotype(pred_biotype)) %>%
  filter(pred_biotype == "lncrna") %>%
  distinct(gene_id) %>% nrow()

# --- 3) circRNA loci, union of 6 h and 24 h ----------------------------------
circ_count <- 0L
inputs <- list(meta_ann = opt$`meta-ann` %||% "data/frozen/hisat2_6_24_meta_ann.csv.gz",
               meta_novel = opt$`meta-novel` %||% "data/frozen/hisat2_6_24_meta_novel_lncRNA.csv.gz")

if (!isTRUE(opt$`skip-circrna`)) {
  c6  <- opt$`circ-6`  %||% find_input("raw_count_matrix_all_6.txt",
                                       in_dir = opt$`in-dir`, optional = TRUE,
                                       what = "circRNA raw counts 6h")
  c24 <- opt$`circ-24` %||% find_input("raw_count_matrix_all_24.txt",
                                       in_dir = opt$`in-dir`, optional = TRUE,
                                       what = "circRNA raw counts 24h")
  if (is.na(c6) || is.na(c24)) {
    stop("The circRNA count matrices were not found.\n",
         "They are produced by scripts/R/04b_circRNA_de_temporal.R, which is a\n",
         "RECONSTRUCTION of a step that no received package contained.\n",
         "Run the circRNA branch first, or pass --skip-circrna to build Figure 2B\n",
         "without the circRNA bar. See metadata/circRNA_input_manifest.tsv.",
         call. = FALSE)
  }
  read_circ_ids <- function(path) {
    df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    unique(as.character(df[[names(df)[1]]]))
  }
  circ_count <- length(union(read_circ_ids(c6), read_circ_ids(c24)))
  inputs$circ_6 <- c6; inputs$circ_24 <- c24
  message("[R06] circRNA unique loci across 6 h and 24 h = ", circ_count)
} else {
  message("[R06] --skip-circrna: Figure 2B will omit the circRNA class.")
}
write_input_manifest(inputs, out_dir)

# --- 4) summary table ---------------------------------------------------------
pc <- ann_counts_main %>% filter(gene_biotype == "protein_coding") %>% pull(count)
if (length(pc) == 0) pc <- 0

rows <- list(
  tibble(gene_biotype = "protein_coding",   count = pc),
  tibble(gene_biotype = "annotated_lncrna", count = annot_lnc_count),
  tibble(gene_biotype = "novel_lncrna",     count = novel_lnc_count))
if (!isTRUE(opt$`skip-circrna`)) {
  rows <- c(rows, list(tibble(gene_biotype = "circrna", count = circ_count)))
}
rows <- c(rows, list(ann_counts_main %>% filter(gene_biotype != "protein_coding")))

summary_tbl <- bind_rows(rows) %>%
  group_by(gene_biotype) %>% summarise(count = sum(count), .groups = "drop") %>%
  mutate(label = case_when(
    gene_biotype == "annotated_lncrna" ~ "annotated lncRNA",
    gene_biotype == "novel_lncrna"     ~ "novel lncRNA",
    gene_biotype == "circrna"          ~ "circRNA",
    TRUE ~ gene_biotype)) %>%
  arrange(desc(count)) %>%
  mutate(label = factor(label, levels = label))

write_out(summary_tbl %>% select(label, count), out_dir,
          "transcriptome_overview_with_circRNA_counts.csv")

# --- 5) plot ------------------------------------------------------------------
p <- ggplot(summary_tbl, aes(x = label, y = count, fill = label)) +
  geom_col(width = 0.92, colour = NA) +
  geom_text(aes(label = comma(count)), vjust = -0.35, size = 4.5) +
  scale_y_log10(labels = label_comma(),
                breaks = c(1, 10, 100, 1000, 10000, 100000),
                expand = expansion(mult = c(0.02, 0.14))) +
  scale_fill_manual(values = grDevices::hcl.colors(nrow(summary_tbl),
                                                   palette = "Spectral", rev = TRUE)) +
  labs(title = "Transcriptome overview, 6 h and 24 h",
       x = "Gene biotype", y = "Count (log10 scale)") +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

save_figure(p, out_dir, "Figure2B_transcriptome_overview_6h_24h_with_circRNA",
            width = 16, height = 7.5, dpi = opt$dpi)

message("[R06] complete -> ", out_dir)
