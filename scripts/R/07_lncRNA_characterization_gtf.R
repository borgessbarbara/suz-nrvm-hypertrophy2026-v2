#!/usr/bin/env Rscript
# =============================================================================
# R07 - Figure 2C and Figure S2A-C: lncRNA structural characterisation.
#
# Parses the reconstructed transcriptome GTF and summarises, for annotated
# (ENSRNOG) versus novel (MSTRG) loci: transcript length, exon count, isoform
# complexity per locus, and the gene/transcript inventory.
#
# Usage:   Rscript scripts/R/07_lncRNA_characterization_gtf.R [--gtf FILE] [--out-dir DIR] [--dpi N]
#          Rscript scripts/R/07_lncRNA_characterization_gtf.R --help
#
# Adapted from archive/final_package_v3_original/01_FINAL_FIGURE_AND_ANALYSIS_SCRIPTS/15_lncRNA_characterization_from_GTF.R
# Changes from the original are listed in docs/CHANGES_FROM_ORIGINALS.md
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(tidyr)
  library(stringr); library(ggplot2); library(scales)
})

here <- dirname(normalizePath(sub("^--file=", "",
         commandArgs(trailingOnly = FALSE)[grepl("^--file=",
         commandArgs(trailingOnly = FALSE))][1]), winslash = "/"))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "plotting.R"))

opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--gtf", type = "character", default = NULL,
                help = "merged annotated GTF [stage 07]"),
    make_option("--in-dir", type = "character", default = "results",
                help = "root to search for the GTF [default %default]"),
    make_option("--out-dir", type = "character",
                default = "results/R07_lncRNA_characterization",
                help = "[default %default]"),
    make_option("--dpi", type = "integer", default = 600,
                help = "figure DPI; the original used 400 [default %default]")),
  description = paste("\nR05: lncRNA structural characterisation from the",
                      "reconstructed transcriptome GTF.\nServes Figure 2C and",
                      "Supplementary Figure S2A-C.")))

out_dir <- ensure_out_dir(opt$`out-dir`)
gtf_file <- opt$gtf %||% find_input(
  c("merged.annotated.gtf", "merged_CNE_6_24_HISAT2.annotated.gtf"),
  in_dir = opt$`in-dir`, what = "merged annotated GTF")
stopifnot(file.exists(gtf_file))
write_input_manifest(list(merged_annotated_gtf = gtf_file), out_dir)
message("[R07] GTF: ", gtf_file)

read_gtf <- function(path) {
  cols <- c("seqname","source_col","feature","start","end","score","strand","frame","attribute")
  readr::read_tsv(path, comment = "#", col_names = cols,
    col_types = readr::cols(
      seqname = col_character(), source_col = col_character(),
      feature = col_character(), start = col_integer(), end = col_integer(),
      score = col_character(), strand = col_character(),
      frame = col_character(), attribute = col_character()),
    progress = FALSE)
}
extract_attr <- function(attr, key) str_match(attr, paste0(key, ' "([^"]+)"'))[, 2]
infer_source <- function(gene_id) case_when(
  str_detect(gene_id, "^MSTRG")   ~ "novel",
  str_detect(gene_id, "^ENSRNOG") ~ "annotated",
  TRUE ~ "other")

gtf <- read_gtf(gtf_file) %>%
  mutate(transcript_id = extract_attr(attribute, "transcript_id"),
         gene_id       = extract_attr(attribute, "gene_id"),
         gene_name     = extract_attr(attribute, "gene_name"),
         ref_gene_id   = extract_attr(attribute, "ref_gene_id"),
         cmp_ref       = extract_attr(attribute, "cmp_ref"),
         class_code    = extract_attr(attribute, "class_code"),
         xloc          = extract_attr(attribute, "xloc"),
         exon_number   = suppressWarnings(as.integer(extract_attr(attribute, "exon_number"))),
         source        = infer_source(gene_id))

transcripts <- gtf %>%
  filter(feature == "transcript", !is.na(transcript_id), !is.na(gene_id)) %>%
  transmute(seqname, strand, transcript_id, gene_id,
            gene_name = coalesce(gene_name, gene_id),
            ref_gene_id, cmp_ref, class_code, xloc, source,
            transcript_start = start, transcript_end = end,
            transcript_length = end - start + 1L) %>%
  distinct()

exons <- gtf %>%
  filter(feature == "exon", !is.na(transcript_id), !is.na(gene_id)) %>%
  transmute(transcript_id, gene_id, exon_start = start, exon_end = end,
            exon_length = end - start + 1L) %>%
  distinct()

exon_summary <- exons %>%
  group_by(transcript_id, gene_id) %>%
  summarise(exon_count = n(),
            summed_exon_length = sum(exon_length, na.rm = TRUE), .groups = "drop")

tx_tbl <- transcripts %>%
  left_join(exon_summary, by = c("transcript_id", "gene_id")) %>%
  mutate(exon_count = replace_na(exon_count, 1L),
         summed_exon_length = coalesce(summed_exon_length, transcript_length),
         transcript_class = case_when(
           source == "annotated" & !is.na(cmp_ref) & !is.na(class_code) & class_code == "=" ~ "known_reference_transcript",
           source == "annotated" & !is.na(cmp_ref) & !is.na(class_code) & class_code != "=" ~ "known_gene_alt_isoform",
           source == "annotated" & is.na(cmp_ref) ~ "annotated_gene_transcript_without_cmp_ref",
           source == "novel" ~ "novel_reconstructed_transcript",
           TRUE ~ "other"))

gene_tbl <- tx_tbl %>%
  group_by(gene_id, source) %>%
  summarise(gene_name = first(gene_name),
            n_transcripts = n_distinct(transcript_id),
            min_tx_length = min(transcript_length, na.rm = TRUE),
            median_tx_length = median(transcript_length, na.rm = TRUE),
            max_tx_length = max(transcript_length, na.rm = TRUE),
            median_exon_count = median(exon_count, na.rm = TRUE),
            class_code_values = paste(sort(unique(na.omit(class_code))), collapse = ","),
            .groups = "drop")

inventory_summary <- bind_rows(
  tx_tbl   %>% count(source, name = "n_transcripts") %>%
               mutate(metric = "total_transcripts") %>% rename(value = n_transcripts),
  gene_tbl %>% count(source, name = "n_genes") %>%
               mutate(metric = "total_genes") %>% rename(value = n_genes),
  tx_tbl   %>% count(source, transcript_class, name = "n") %>%
               mutate(metric = paste0("transcript_class__", transcript_class)) %>%
               select(source, metric, value = n))

isoform_summary <- gene_tbl %>%
  mutate(isoform_bin = case_when(
    n_transcripts == 1 ~ "1 transcript",
    n_transcripts == 2 ~ "2 transcripts",
    n_transcripts == 3 ~ "3 transcripts",
    n_transcripts >= 4 ~ "4+ transcripts",
    TRUE ~ "other")) %>%
  count(source, isoform_bin, name = "n_genes")

write_out_tsv(tx_tbl,            out_dir, "lncRNA_transcript_characterization_from_GTF.tsv")
write_out_tsv(gene_tbl,          out_dir, "lncRNA_gene_characterization_from_GTF.tsv")
write_out_tsv(inventory_summary, out_dir, "lncRNA_inventory_summary.tsv")
write_out_tsv(isoform_summary,   out_dir, "lncRNA_isoform_complexity_summary.tsv")

plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_text(face = "bold"))

p_len <- tx_tbl %>%
  filter(source %in% c("annotated", "novel")) %>%
  mutate(source = factor(source, levels = c("annotated", "novel"))) %>%
  ggplot(aes(x = source, y = transcript_length, fill = source)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.9, colour = "grey30") +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "Transcript length (nt)",
       title = "Annotated vs novel lncRNA transcript length") +
  plot_theme + theme(legend.position = "none")
save_figure(p_len, out_dir, "FigureS2C_lncRNA_transcript_length_violin",
            width = 6.4, height = 5.1, dpi = opt$dpi)

p_exon <- tx_tbl %>%
  filter(source %in% c("annotated", "novel")) %>%
  count(source, exon_count, name = "n") %>%
  mutate(source = factor(source, levels = c("annotated", "novel"))) %>%
  ggplot(aes(x = exon_count, y = n, fill = source)) +
  geom_col(position = "dodge", width = 0.8) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(labels = comma) +
  labs(x = "Exons per transcript", y = "# transcripts",
       title = "Exon count distribution of annotated vs novel lncRNA transcripts") +
  plot_theme
save_figure(p_exon, out_dir, "FigureS2A_lncRNA_exon_count_distribution",
            width = 6.8, height = 5.1, dpi = opt$dpi)

p_inv <- inventory_summary %>%
  filter(metric %in% c("total_transcripts", "total_genes")) %>%
  mutate(metric = recode(metric, total_transcripts = "Transcripts", total_genes = "Genes"),
         source = factor(source, levels = c("annotated", "novel"))) %>%
  ggplot(aes(x = metric, y = value, fill = source)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(aes(label = comma(value)), position = position_dodge(width = 0.72),
            vjust = -0.25, size = 3.5) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Count",
       title = "Annotated and novel lncRNA inventory from transcriptome reconstruction") +
  plot_theme
save_figure(p_inv, out_dir, "Figure2C_lncRNA_gene_transcript_inventory",
            width = 6.8, height = 5.1, dpi = opt$dpi)

p_iso <- isoform_summary %>%
  mutate(source = factor(source, levels = c("annotated", "novel")),
         isoform_bin = factor(isoform_bin,
           levels = c("1 transcript", "2 transcripts", "3 transcripts", "4+ transcripts"))) %>%
  ggplot(aes(x = source, y = n_genes, fill = isoform_bin)) +
  geom_col(width = 0.72) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "# lncRNA genes", fill = "Isoform complexity",
       title = "Isoform complexity per lncRNA gene") +
  plot_theme
save_figure(p_iso, out_dir, "FigureS2B_lncRNA_isoform_complexity",
            width = 6.8, height = 5.1, dpi = opt$dpi)

writeLines(c(
  "lncRNA transcriptome reconstruction summary",
  "",
  paste0("Annotated lncRNA genes: ",       comma(sum(gene_tbl$source == "annotated"))),
  paste0("Novel lncRNA genes: ",           comma(sum(gene_tbl$source == "novel"))),
  paste0("Annotated lncRNA transcripts: ", comma(sum(tx_tbl$source == "annotated"))),
  paste0("Novel lncRNA transcripts: ",     comma(sum(tx_tbl$source == "novel"))),
  "",
  paste0("Annotated median transcript length (nt): ",
         round(median(tx_tbl$transcript_length[tx_tbl$source == "annotated"], na.rm = TRUE), 1)),
  paste0("Novel median transcript length (nt): ",
         round(median(tx_tbl$transcript_length[tx_tbl$source == "novel"], na.rm = TRUE), 1)),
  paste0("Annotated median exon count: ",
         round(median(tx_tbl$exon_count[tx_tbl$source == "annotated"], na.rm = TRUE), 1)),
  paste0("Novel median exon count: ",
         round(median(tx_tbl$exon_count[tx_tbl$source == "novel"], na.rm = TRUE), 1)),
  "",
  "Interpretation note:",
  paste("Because this is a transcriptome reconstruction framework, novel MSTRG",
        "entries are best described conservatively as novel reconstructed lncRNA",
        "loci/transcripts rather than definitively novel lncRNA genes unless",
        "orthogonal evidence supports locus-level novelty.")),
  file.path(out_dir, "lncRNA_reconstruction_summary.txt"))
message("[out] lncRNA_reconstruction_summary.txt")

message("[R07] complete -> ", out_dir)
