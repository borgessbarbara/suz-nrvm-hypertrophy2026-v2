#!/usr/bin/env Rscript
# =============================================================================
# R17 — Figure S3A-B mRNA and miRNA: volcano plots.
#
#   *** NEW SCRIPT. The final plotting code was embedded in the large
#   *** historical scripts and was never isolated.
#
# 04_METADATA/final_figure_script_map.tsv marks this panel OPEN:
#     "Figure S3 | A-B mRNA/miRNA | Volcano plots |
#      Final standalone scripts not isolated | OPEN |
#      Code exists in larger historical scripts; collaborator should
#      extract/parameterize only final plotting logic."
# and remaining_open_items.tsv repeats it as category S3_MRNA_MIRNA.
#
# The equivalent lncRNA panel is 07_Figure4_clean_volcanoes.R, ported as R15.
# This script deliberately mirrors R15's conventions so the S3A-B panels are
# consistent once assembled:
#   * padj threshold for mRNA and miRNA (FDR), as in the DESeq2 stages;
#   * the same UP/DOWN palette from lib/plotting.R;
#   * top-N labelling per direction.
#
# NOTE ON THE lncRNA PANEL'S THRESHOLD
# ------------------------------------
# 07_Figure4_clean_volcanoes.R sets `absfc_cutoff <- 0`, i.e. the lncRNA
# volcanoes call significance on padj ALONE with no fold-change requirement.
# This script defaults to --lfc-cut 0 for the same reason, so the three S3
# panels use one rule. Override if the manuscript says otherwise.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(stringr)
  library(tibble); library(readr); library(ggplot2); library(ggrepel); library(scales)
})

.here <- dirname(normalizePath(sub("^--file=", "",
          commandArgs(trailingOnly = FALSE)[grepl("^--file=",
          commandArgs(trailingOnly = FALSE))][1]), winslash = "/"))
source(file.path(.here, "lib", "io.R"))
source(file.path(.here, "lib", "plotting.R"))

opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--in-dir", type = "character", default = "results",
                help = "root searched for the DE tables [default %default]"),
    make_option("--out-dir", type = "character",
                default = "results/R17_volcanoes_mRNA_miRNA",
                help = "[default %default]"),
    make_option("--padj-cut", type = "double", default = 0.05, help = "[default %default]"),
    make_option("--lfc-cut", type = "double", default = 0,
                help = "0 matches the lncRNA panel (absfc_cutoff <- 0) [default %default]"),
    make_option("--labels-per-direction", type = "integer", default = 6,
                help = "matches R15 [default %default]"),
    make_option("--dpi", type = "integer", default = 600, help = "[default %default]")),
  description = paste("\nR16: Figure S3A-B mRNA and miRNA volcano plots.",
                      "\nNEW - extracted from the historical scripts (audit item S3_MRNA_MIRNA).")))

out_dir <- ensure_out_dir(opt$`out-dir`)

# --- locate the four DE tables ------------------------------------------------
spec <- tribble(
  ~rna,    ~timepoint, ~candidates,
  "mRNA",  "6h",       "NE_vs_Ctrl_6h_mRNA_DE.txt",
  "mRNA",  "24h",      "NE_vs_Ctrl_24h_mRNA_DE.txt",
  "miRNA", "6h",       "DE_miRNA_All_with_status_NE6.txt",
  "miRNA", "24h",      "DE_miRNA_All_with_status_NE24.txt")

paths <- vapply(spec$candidates, function(f)
  find_input(f, in_dir = opt$`in-dir`, what = f), character(1))
spec$path <- unname(paths)
write_input_manifest(setNames(as.list(spec$path),
                              paste0(spec$rna, "_", spec$timepoint)), out_dir)

read_one <- function(path, rna, tp) {
  x <- read_auto(path)
  id  <- choose_col(x,
          c("gene_id", "ensembl_gene_id", "gene", "miRNA", "miRNA_id", "mirna", "Row.names"),
          c("^gene.?id$", "^mir", "^id$"), what = paste(rna, "identifier"))
  lfc <- choose_col(x, c("log2FoldChange", "log2FC", "logFC"), "log2.?f", what = "log2FC")
  padj <- choose_col(x, c("padj", "FDR", "adj.P.Val", "p_adj"), "adj", what = "adjusted p")
  nm  <- choose_col(x, c("gene_name", "external_gene_name", "symbol", "Symbol"),
                    "name|symbol", required = FALSE)
  tibble(feature = as.character(x[[id]]),
         label_name = if (is.na(nm)) as.character(x[[id]]) else
                      ifelse(is.na(x[[nm]]) | x[[nm]] == "", as.character(x[[id]]),
                             as.character(x[[nm]])),
         log2FC = suppressWarnings(as.numeric(x[[lfc]])),
         padj   = suppressWarnings(as.numeric(x[[padj]])),
         rna = rna, timepoint = tp)
}

de <- bind_rows(lapply(seq_len(nrow(spec)), function(i)
  read_one(spec$path[i], spec$rna[i], spec$timepoint[i])))

de <- de %>%
  filter(!is.na(log2FC)) %>%
  mutate(direction = case_when(
           !is.na(padj) & padj < opt$`padj-cut` & log2FC >   opt$`lfc-cut` ~ "UP",
           !is.na(padj) & padj < opt$`padj-cut` & log2FC < - opt$`lfc-cut` ~ "DOWN",
           TRUE ~ "NS"),
         neglog10padj = -log10(pmax(padj, .Machine$double.xmin)),
         timepoint = factor(timepoint, levels = c("6h", "24h")))

write_out(de, out_dir, "SUPP_FigS3AB_mRNA_miRNA_volcano_data.csv")
write_out(de %>% count(rna, timepoint, direction, name = "n"), out_dir,
          "SUPP_FigS3AB_counts_by_rna_time_direction.csv")

lab <- de %>%
  filter(direction != "NS") %>%
  group_by(rna, timepoint, direction) %>%
  slice_max(neglog10padj, n = opt$`labels-per-direction`, with_ties = FALSE) %>%
  ungroup()

for (r in unique(de$rna)) {
  d <- de %>% filter(rna == r)
  l <- lab %>% filter(rna == r)
  p <- ggplot(d, aes(x = log2FC, y = neglog10padj, colour = direction)) +
    geom_point(data = d %>% filter(direction == "NS"),
               colour = suz_direction_colours[["NS"]], alpha = 0.45, size = 0.8) +
    geom_point(data = d %>% filter(direction != "NS"), alpha = 0.8, size = 1.1) +
    geom_hline(yintercept = -log10(opt$`padj-cut`), linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    { if (opt$`lfc-cut` > 0)
        geom_vline(xintercept = c(-opt$`lfc-cut`, opt$`lfc-cut`), linetype = "dashed",
                   colour = "grey40", linewidth = 0.3) else NULL } +
    ggrepel::geom_text_repel(data = l, aes(label = label_name), size = 2.6,
                             max.overlaps = 25, min.segment.length = 0,
                             segment.size = 0.2, show.legend = FALSE) +
    scale_colour_manual(values = suz_direction_colours, name = NULL) +
    facet_wrap(~ timepoint, nrow = 1, scales = "free") +
    labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~p),
         title = paste0("Figure S3A-B - ", r, " differential expression"),
         subtitle = sprintf("padj < %.2g%s", opt$`padj-cut`,
                            if (opt$`lfc-cut` > 0)
                              sprintf(" and |log2FC| > %.2g", opt$`lfc-cut`)
                            else "  (no fold-change cutoff, matching the lncRNA panel)")) +
    suz_theme(base_size = 11) +
    theme(legend.position = "bottom")

  save_figure(p, out_dir, paste0("FigureS3AB_", r, "_volcano"),
              width = 9, height = 4.8, dpi = opt$dpi)
}

message("[R17] complete -> ", out_dir)
