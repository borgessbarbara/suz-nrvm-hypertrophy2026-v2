#!/usr/bin/env Rscript
# =============================================================================
# R13 — Figure 3C circRNA: exploratory temporal heatmap.
#
#   *** NEW SCRIPT. No standalone original existed. ***
#
# 04_METADATA/final_figure_script_map.tsv marks this panel OPEN:
#     "Figure 3 | C circRNA | Exploratory temporal heatmap |
#      Exact standalone script not isolated | OPEN |
#      Repository may need a small dedicated script from the final circRNA
#      DE/temporal table."
# and remaining_open_items.tsv repeats it as category FIG3C_CIRC. This is that
# script.
#
# It deliberately mirrors the structure of the lncRNA and mRNA temporal
# heatmaps (R09, R10) so the three Figure 3C panels are directly comparable
# once assembled: same phase definitions, same z-scoring, same ordering.
#
# EXPLORATORY STATUS
# ------------------
# 31_phase_aware_ceRNA_coherence_plots.R states that "circRNA trajectories
# remain exploratory", and the circRNA DE uses a NOMINAL p threshold rather
# than FDR. That is carried through here: the panel is labelled exploratory
# and the caption records the thresholds actually applied.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(stringr)
  library(tibble); library(readr); library(ggplot2); library(scales)
})

.here <- dirname(normalizePath(sub("^--file=", "",
          commandArgs(trailingOnly = FALSE)[grepl("^--file=",
          commandArgs(trailingOnly = FALSE))][1]), winslash = "/"))
source(file.path(.here, "lib", "io.R"))
source(file.path(.here, "lib", "plotting.R"))

opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--de-6", type = "character", default = NULL,
                help = "06_circRNA_sequences_de6_DE_complete_table.txt [R13]"),
    make_option("--de-24", type = "character", default = NULL,
                help = "07_circRNA_sequences_de24_DE_complete_table.txt [R13]"),
    make_option("--in-dir", type = "character", default = "results",
                help = "[default %default]"),
    make_option("--out-dir", type = "character",
                default = "results/R13_figure3C_circRNA_heatmap",
                help = "[default %default]"),
    make_option("--top-per-phase", type = "integer", default = 12,
                help = "circRNAs shown per temporal phase [default %default]"),
    make_option("--alpha", type = "double", default = 0.05,
                help = "nominal p threshold, as in R19 [default %default]"),
    make_option("--lfc-cut", type = "double", default = 0.5,
                help = "[default %default]"),
    make_option("--dpi", type = "integer", default = 600, help = "[default %default]")),
  description = paste("\nR11: Figure 3C circRNA exploratory temporal heatmap.",
                      "\nNEW - no standalone original existed (audit item FIG3C_CIRC).")))

out_dir <- ensure_out_dir(opt$`out-dir`)

de6_path  <- opt$`de-6`  %||% find_input("06_circRNA_sequences_de6_DE_complete_table.txt",
                                         in_dir = opt$`in-dir`, what = "circRNA DE table 6h")
de24_path <- opt$`de-24` %||% find_input("07_circRNA_sequences_de24_DE_complete_table.txt",
                                         in_dir = opt$`in-dir`, what = "circRNA DE table 24h")
write_input_manifest(list(de_6h = de6_path, de_24h = de24_path), out_dir)

read_de <- function(path, tp) {
  x <- read_auto(path)
  id  <- choose_col(x, c("circRNA_ID", "circRNA", "id"), "^circ.*id$|^id$", what = "circRNA id")
  lfc <- choose_col(x, c("log2FoldChange", "log2FC"), "log2.?f", what = "log2FC")
  pv  <- choose_col(x, c("pvalue", "p_value", "P.Value"), "^p.?val", what = "p value")
  hg  <- choose_col(x, c("host_gene", "geneName", "gene_name"), "gene", required = FALSE)
  tibble(circRNA_ID = as.character(x[[id]]),
         host_gene  = if (is.na(hg)) NA_character_ else as.character(x[[hg]]),
         log2FC     = suppressWarnings(as.numeric(x[[lfc]])),
         pvalue     = suppressWarnings(as.numeric(x[[pv]])),
         timepoint  = tp)
}

de <- bind_rows(read_de(de6_path, "6h"), read_de(de24_path, "24h")) %>%
  mutate(direction = case_when(
    !is.na(pvalue) & pvalue < opt$alpha & log2FC >=  opt$`lfc-cut` ~ "UP",
    !is.na(pvalue) & pvalue < opt$alpha & log2FC <= -opt$`lfc-cut` ~ "DOWN",
    TRUE ~ "NS"))

# --- temporal phase assignment, matching the mRNA/lncRNA convention ----------
#   Early     significant at 6 h only
#   Sustained significant at both, same direction
#   Late      significant at 24 h only
#   Reversal  significant at both, opposite directions
wide <- de %>%
  select(circRNA_ID, host_gene, timepoint, log2FC, direction) %>%
  pivot_wider(id_cols = c(circRNA_ID, host_gene),
              names_from = timepoint,
              values_from = c(log2FC, direction)) %>%
  mutate(
    d6  = replace_na(direction_6h,  "NS"),
    d24 = replace_na(direction_24h, "NS"),
    phase = case_when(
      d6 != "NS" & d24 == "NS"              ~ "Early",
      d6 != "NS" & d24 != "NS" & d6 == d24  ~ "Sustained",
      d6 == "NS" & d24 != "NS"              ~ "Late",
      d6 != "NS" & d24 != "NS" & d6 != d24  ~ "Reversal",
      TRUE ~ "NS"))

write_out(wide, out_dir, "00_circRNA_temporal_classification.csv")
write_out(wide %>% count(phase, d6, d24, name = "n_circRNA"), out_dir,
          "01_circRNA_temporal_phase_counts.csv")

# The per-phase circRNA lists. 20_ceRNA_coherent_pairs_and_triplets reads
# early/sustained/late/reversal_circRNA.csv alongside the equivalent mRNA and
# lncRNA phase files, but nothing in the material received produced the
# circRNA ones -- the temporal classification existed only inside whichever
# script drew the heatmap. Emitting them here closes that gap.
for (ph in c("Early", "Sustained", "Late", "Reversal")) {
  write_out(wide %>% filter(phase == ph) %>%
              select(circRNA_ID, host_gene, log2FC_6h, log2FC_24h,
                     direction_6h = d6, direction_24h = d24),
            out_dir, paste0(tolower(ph), "_circRNA.csv"))
}

plot_tbl <- wide %>% filter(phase %in% c("Early", "Sustained", "Late", "Reversal"))
if (!nrow(plot_tbl)) {
  stop("No circRNA passed the exploratory DE thresholds (nominal p < ", opt$alpha,
       ", |log2FC| >= ", opt$`lfc-cut`, "), so Figure 3C circRNA cannot be drawn.",
       call. = FALSE)
}

top <- plot_tbl %>%
  mutate(max_abs = pmax(abs(replace_na(log2FC_6h, 0)), abs(replace_na(log2FC_24h, 0)))) %>%
  group_by(phase) %>% slice_max(max_abs, n = opt$`top-per-phase`, with_ties = FALSE) %>%
  ungroup()

long <- top %>%
  select(circRNA_ID, host_gene, phase, log2FC_6h, log2FC_24h) %>%
  pivot_longer(c(log2FC_6h, log2FC_24h), names_to = "timepoint", values_to = "log2FC") %>%
  mutate(timepoint = factor(ifelse(timepoint == "log2FC_6h", "6 h", "24 h"),
                            levels = c("6 h", "24 h")),
         phase = factor(phase, levels = c("Early", "Sustained", "Late", "Reversal")),
         label = ifelse(is.na(host_gene) | host_gene == "", circRNA_ID,
                        paste0(circRNA_ID, "  (", host_gene, ")")))

ord <- long %>% group_by(label, phase) %>%
  summarise(k = max(log2FC, na.rm = TRUE), .groups = "drop") %>%
  arrange(phase, k)
long$label <- factor(long$label, levels = ord$label)

write_out(long, out_dir, "02_circRNA_temporal_heatmap_table.csv")

lim <- max(abs(long$log2FC), na.rm = TRUE)
p <- ggplot(long, aes(x = timepoint, y = label, fill = log2FC)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(phase ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = suz_direction_colours[["DOWN"]], mid = "white",
                       high = suz_direction_colours[["UP"]],
                       midpoint = 0, limits = c(-lim, lim),
                       name = expression(log[2]~FC)) +
  labs(x = NULL, y = NULL,
       title = "Figure 3C - circRNA temporal programmes (exploratory)",
       subtitle = sprintf("nominal p < %.2g and |log2FC| >= %.2g; top %d per phase",
                          opt$alpha, opt$`lfc-cut`, opt$`top-per-phase`)) +
  suz_theme(base_size = 10) +
  theme(axis.text.y = element_text(size = 7),
        strip.text.y = element_text(angle = 0, face = "bold"),
        panel.grid = element_blank())

save_figure(p, out_dir, "Figure3C_circRNA_temporal_heatmap",
            width = 7.5, height = max(5, 0.19 * nlevels(long$label) + 2),
            dpi = opt$dpi)

writeLines(c(
  "Figure 3C circRNA - exploratory temporal heatmap",
  "",
  "This panel had no standalone script in any package received; it is new,",
  "written to the specification in 04_METADATA/final_figure_script_map.tsv",
  "(audit item FIG3C_CIRC).",
  "",
  sprintf("Thresholds: nominal p < %.3g, |log2FC| >= %.3g", opt$alpha, opt$`lfc-cut`),
  sprintf("Top %d circRNAs per phase by maximum |log2FC|.", opt$`top-per-phase`),
  "",
  "circRNA differential expression is treated as EXPLORATORY throughout the",
  "manuscript: the threshold is nominal p, not FDR."),
  file.path(out_dir, "README_Figure3C_circRNA.txt"))

message("[R13] complete -> ", out_dir)
