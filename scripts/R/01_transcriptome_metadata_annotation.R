#!/usr/bin/env Rscript
# =============================================================================
# R01 - Transcript metadata and annotation.
#
# Builds, for the HISAT2 6+24 h dataset:
#   <prefix>_meta_ann.csv              annotated transcripts with exon counts
#   <prefix>_meta_novel_lncRNA.csv     novel lncRNA loci (class codes u, i, x)
#   <prefix>_exon_lengths.csv          per-exon table
#   <prefix>_incompatible_lncRNA.csv   lncRNA predictions Ensembl calls otherwise
#   <prefix>.compatible.gtf            the retained transcript set
#   <prefix>_feelnc_validated.csv      FEELnc relations for validated lncRNAs
#   gtf_outputs/<prefix>_final_annotated.gtf
#
# This is the only stage that queries Ensembl. The archive host is pinned in
# scripts/R/lib/ensembl.R, driven by metadata/reference_release.tsv.
#
# Usage:   Rscript scripts/R/01_transcriptome_metadata_annotation.R [--in-dir DIR] [--out-dir DIR]
#          Rscript scripts/R/01_transcriptome_metadata_annotation.R --dry-run-biomart
#          Rscript scripts/R/01_transcriptome_metadata_annotation.R --help
#
# Requires bash stages 07 and 09 through 09e.
#
# Adapted from the transcriptome-metadata notebook (B. Borges). The shell steps
# it depends on are in archive/transcriptome_metadata_original/.
# Changes from the original are listed in docs/CHANGES_FROM_ORIGINALS.md
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
})

here <- dirname(normalizePath(sub("^--file=", "",
         commandArgs(trailingOnly = FALSE)[grepl("^--file=",
         commandArgs(trailingOnly = FALSE))][1]), winslash = "/"))
source(file.path(here, "lib", "io.R"))
source(file.path(here, "lib", "ensembl.R"))

# --- CLI ---------------------------------------------------------------------
option_list <- list(
  make_option("--tmap", type = "character", default = NULL,
              help = "gffcompare .tmap file [stage 07]"),
  make_option("--feelnc-classes", type = "character", default = NULL,
              help = "FEELnc classifier output, lncRNA_classes_*_1kb_100kb.txt [stage 09]"),
  make_option("--lncrna-gtf", type = "character", default = NULL,
              help = "<id>.annotated.lncRNA_tx.gtf [stage 09d]"),
  make_option("--other-rna-gtf", type = "character", default = NULL,
              help = "filtered_merged_<id>.annotated.gtf [stage 09e]"),
  make_option("--in-dir", type = "character", default = "results",
              help = "root to search for inputs not given explicitly [default %default]"),
  make_option("--out-dir", type = "character", default = "results/R01_transcriptome_metadata",
              help = "output directory [default %default]"),
  make_option("--prefix", type = "character", default = "hisat2_6_24",
              help = "output filename prefix [default %default]"),
  make_option("--biomart-host", type = "character", default = NULL,
              help = paste0("Ensembl archive host [default ",
                            SUZ_ENSEMBL_HOST_DEFAULT, "]")),
  make_option("--dry-run-biomart", action = "store_true", default = FALSE,
              help = "Only verify the pinned Ensembl host resolves, then exit.")
)
opt <- parse_args(OptionParser(
  option_list = option_list,
  description = paste(
    "\nR01: build the transcript metadata tables (_meta_ann.csv and",
    "\n_meta_novel_lncRNA.csv) plus exon lengths, the compatible GTF, the",
    "\nFEELnc-validated relations and the final annotated GTF.",
    "\n\nRequires bash stages 07, 09, 09b, 09c, 09d and 09e.")))

# --- smoke test ---------------------------------------------------------------
if (isTRUE(opt$`dry-run-biomart`)) {
  message("[R01] --dry-run-biomart: verifying the pinned Ensembl archive only.")
  info <- check_ensembl_connection(host = opt$`biomart-host`)
  message("[R01] OK. host=", info$host, "  dataset=", info$dataset,
          "  version=", info$version)
  quit(status = 0)
}

out_dir <- ensure_out_dir(opt$`out-dir`)
gtf_out_dir <- ensure_out_dir(file.path(out_dir, "gtf_outputs"))
prefix <- opt$prefix

inputs <- list(
  tmap          = opt$tmap          %||% find_input(sprintf("CNE_6_24_HISAT2.tmap"),
                                                    in_dir = opt$`in-dir`, what = "gffcompare .tmap"),
  feelnc        = opt$`feelnc-classes` %||% find_input(
                     sprintf("lncRNA_classes_%s_1kb_100kb.txt", "CNE_6_24_HISAT2"),
                     in_dir = opt$`in-dir`, what = "FEELnc classes"),
  lncrna_gtf    = opt$`lncrna-gtf`  %||% find_input("CNE_6_24_HISAT2.annotated.lncRNA_tx.gtf",
                                                    in_dir = opt$`in-dir`, what = "lncRNA tx GTF"),
  other_rna_gtf = opt$`other-rna-gtf` %||% find_input("filtered_merged_CNE_6_24_HISAT2.annotated.gtf",
                                                    in_dir = opt$`in-dir`, what = "lncRNA-depleted merged GTF")
)
write_input_manifest(inputs, out_dir)

summaries <- list()
sw <- function(name, df) {
  summaries[[name]] <<- df
  write_out(df, out_dir, paste0("SUMMARY_", name, ".csv"))
}

# --- Step 1: load -------------------------------------------------------------
message("[R01] Step 1/8  loading inputs")
tmap_data         <- read.table(inputs$tmap,   sep = "\t", header = TRUE)
feelnc_data       <- read.table(inputs$feelnc, sep = "\t", header = TRUE)
lncrna_gtf_raw    <- rtracklayer::import(inputs$lncrna_gtf)
other_rna_gtf_raw <- rtracklayer::import(inputs$other_rna_gtf)

sw("01_load", data.frame(
  File = c("TMAP", "FEELnc", "lncRNA GTF", "Other RNA GTF"),
  Entries = c(nrow(tmap_data), nrow(feelnc_data),
              length(lncrna_gtf_raw), length(other_rna_gtf_raw)),
  Unique_Transcripts = c(length(unique(tmap_data$qry_id)),
                         length(unique(feelnc_data$lncRNA_transcript)),
                         length(unique(lncrna_gtf_raw$transcript_id)),
                         length(unique(other_rna_gtf_raw$transcript_id)))))

# --- Step 2: strand filter + stranded exports ---------------------------------
message("[R01] Step 2/8  strand filtering")
lncrna_gtf    <- lncrna_gtf_raw[strand(lncrna_gtf_raw) != "*"]
other_rna_gtf <- other_rna_gtf_raw[strand(other_rna_gtf_raw) != "*"]

export(lncrna_gtf,
       file.path(out_dir, sub("\\.gtf$", ".stranded.gtf", basename(inputs$lncrna_gtf))),
       format = "gtf")
export(other_rna_gtf,
       file.path(out_dir, sub("\\.gtf$", ".stranded.gtf", basename(inputs$other_rna_gtf))),
       format = "gtf")

sw("02_strand_filter", data.frame(
  RNA_Type = c("lncRNA", "Other RNA"),
  Raw_Count = c(length(unique(lncrna_gtf_raw$transcript_id)),
                length(unique(other_rna_gtf_raw$transcript_id))),
  Stranded_Count = c(length(unique(lncrna_gtf$transcript_id)),
                     length(unique(other_rna_gtf$transcript_id))),
  Removed = c(length(unique(lncrna_gtf_raw$transcript_id)) - length(unique(lncrna_gtf$transcript_id)),
              length(unique(other_rna_gtf_raw$transcript_id)) - length(unique(other_rna_gtf$transcript_id)))))

lncrna_ids    <- unique(mcols(lncrna_gtf)$transcript_id)
other_rna_ids <- unique(mcols(other_rna_gtf)$transcript_id)

tmap_lncrna    <- tmap_data %>% filter(qry_id %in% lncrna_ids)    %>% mutate(pred_biotype = "lncRNA")
tmap_other_rna <- tmap_data %>% filter(qry_id %in% other_rna_ids) %>% mutate(pred_biotype = NA_character_)
tmap_combined  <- bind_rows(tmap_lncrna, tmap_other_rna)

# --- Step 2a: class-code summary ----------------------------------------------
sw("02a_class_codes", tmap_combined %>%
     group_by(pred_biotype, class_code) %>%
     summarise(count = n(), .groups = "drop") %>%
     pivot_wider(names_from = pred_biotype, values_from = count, values_fill = 0) %>%
     rename_with(~ ifelse(. == "NA", "other_RNA", .), .cols = everything()))

# --- Step 3: biomaRt ----------------------------------------------------------
message("[R01] Step 3/8  biomaRt annotation")
tx_attributes <- c("chromosome_name", "ensembl_gene_id", "ensembl_transcript_id",
                   "external_transcript_name", "external_gene_name",
                   "transcript_start", "transcript_end", "transcript_length",
                   "gene_biotype", "transcript_biotype", "strand")

ensembl_mart <- connect_ensembl(host = opt$`biomart-host`)
annotated_tx_ids <- tmap_combined %>%
  filter(class_code %in% c("=", "c", "k")) %>% pull(qry_id)

meta_annotated <- suz_getBM(attributes = tx_attributes,
                            filters = "ensembl_transcript_id",
                            values = annotated_tx_ids, mart = ensembl_mart)

if (nrow(meta_annotated) > 0) {
  meta_annotated$strand <- ifelse(meta_annotated$strand == 1, "+", "-")
} else {
  stop("biomaRt returned no records for ", length(unique(annotated_tx_ids)),
       " transcript IDs. The historical run retrieved ~95,000. Check the host ",
       "and that the transcript IDs are Ensembl IDs for this release.",
       call. = FALSE)
}

sw("03_biomart", data.frame(
  Step = c("Transcripts sent to biomaRt", "Metadata records retrieved",
           "Ensembl host"),
  Count = c(length(unique(annotated_tx_ids)), nrow(meta_annotated),
            NA_integer_),
  Detail = c(NA, NA, suz_ensembl_host(opt$`biomart-host`))))

meta_annotated <- left_join(
  meta_annotated,
  tmap_combined %>% select(qry_id, class_code, len, ref_match_len, pred_biotype),
  by = c("ensembl_transcript_id" = "qry_id"))

# --- Step 4: biotype incompatibility + compatible GTF -------------------------
message("[R01] Step 4/8  biotype compatibility")
incompatible_lncrna <- meta_annotated %>%
  filter(pred_biotype == "lncRNA" & transcript_biotype != "lncRNA")
write_out(incompatible_lncrna, out_dir, paste0(prefix, "_incompatible_lncRNA.csv"))

ids_to_remove <- incompatible_lncrna %>%
  filter(!class_code %in% c("=", "c")) %>% pull(ensembl_transcript_id)
meta_annotated_final <- meta_annotated %>%
  filter(!ensembl_transcript_id %in% ids_to_remove)

compatible_ids <- unique(meta_annotated_final$ensembl_transcript_id)
all_known_stranded_gtf <- c(lncrna_gtf, other_rna_gtf)
compatible_gtf_gr <- all_known_stranded_gtf[
  mcols(all_known_stranded_gtf)$transcript_id %in% compatible_ids]
export(compatible_gtf_gr, file.path(out_dir, paste0(prefix, ".compatible.gtf")),
       format = "gtf")

sw("04_incompatibility", data.frame(
  Category = c("Total Annotated", "Incompatible lncRNA",
               "Incompatible (class 'c','k') removed", "Final Annotated Set",
               "Transcripts in Compatible GTF"),
  Count = c(nrow(meta_annotated), nrow(incompatible_lncrna),
            length(ids_to_remove), nrow(meta_annotated_final),
            length(unique(compatible_gtf_gr$transcript_id)))))

# --- Step 5: exons from the GTF -> _meta_ann.csv ------------------------------
message("[R01] Step 5/8  exon lengths and counts (from the GTF, not biomaRt)")
exons_gr <- compatible_gtf_gr[mcols(compatible_gtf_gr)$type == "exon"]
exon_lengths_df <- NULL

if (length(exons_gr) > 0) {
  biotype_lookup <- meta_annotated_final %>%
    select(ensembl_transcript_id, gene_biotype, transcript_biotype) %>% distinct()
  exon_lengths_df <- as.data.frame(exons_gr) %>%
    dplyr::select(transcript_id, exon_id = exon_number, seqnames, start, end,
                  strand, width) %>%
    left_join(biotype_lookup, by = c("transcript_id" = "ensembl_transcript_id"))
  write_out(exon_lengths_df, out_dir, paste0(prefix, "_exon_lengths.csv"))
  exon_counts <- exon_lengths_df %>%
    group_by(transcript_id) %>% summarise(num_exons = n(), .groups = "drop")
} else {
  exon_counts <- data.frame(transcript_id = character(), num_exons = integer())
}

meta_annotated_with_exons <- meta_annotated_final %>%
  left_join(exon_counts, by = c("ensembl_transcript_id" = "transcript_id")) %>%
  mutate(num_exons = replace_na(num_exons, 0))

sw("05_exons", data.frame(
  Description = c("Known transcripts for exon search", "Total exons found in GTF",
                  "Transcripts with exon info"),
  Count = c(nrow(meta_annotated_final),
            if (!is.null(exon_lengths_df)) nrow(exon_lengths_df) else 0L,
            nrow(exon_counts))))
sw("05_biotypes", meta_annotated_with_exons %>%
     group_by(transcript_biotype) %>% summarise(Count = n(), .groups = "drop") %>%
     arrange(desc(Count)))

write_out(meta_annotated_with_exons, out_dir, paste0(prefix, "_meta_ann.csv"))

# --- Step 6: novel lncRNAs -> _meta_novel_lncRNA.csv --------------------------
message("[R01] Step 6/8  novel lncRNA characterisation")
meta_novel_lncrna <- tmap_lncrna[
  !tmap_lncrna$qry_id %in% meta_annotated_final$ensembl_transcript_id, ]
meta_novel_lncrna <- meta_novel_lncrna[
  meta_novel_lncrna$class_code %in% c("u", "i", "x"), ]

novel_lncrna_gtf_subset <- lncrna_gtf[
  mcols(lncrna_gtf)$transcript_id %in% meta_novel_lncrna$qry_id]
novel_lncrna_gtf_df <- as.data.frame(novel_lncrna_gtf_subset) %>%
  dplyr::select(seqnames, type, transcript_id, gene_id, start, end, strand)

meta_novel_lncrna_final <- merge(novel_lncrna_gtf_df, meta_novel_lncrna,
                                 by.x = "transcript_id", by.y = "qry_id",
                                 all.x = TRUE)
is_exon <- meta_novel_lncrna_final$type == "exon"
meta_novel_lncrna_final$len[is_exon] <-
  meta_novel_lncrna_final$end[is_exon] - meta_novel_lncrna_final$start[is_exon] + 1

sw("06_novel", data.frame(
  Description = "Novel lncRNA transcripts identified",
  Count = length(unique(meta_novel_lncrna_final$transcript_id))))
write_out(meta_novel_lncrna_final, out_dir,
          paste0(prefix, "_meta_novel_lncRNA.csv"))

# --- Step 7: FEELnc validation ------------------------------------------------
message("[R01] Step 7/8  FEELnc validation")
known_lncrna_ids <- meta_annotated_with_exons %>%
  filter(transcript_biotype == "lncRNA") %>% pull(ensembl_transcript_id)
novel_lncrna_ids <- unique(meta_novel_lncrna_final$transcript_id)
validated_lncrna_ids <- c(known_lncrna_ids, novel_lncrna_ids)

feelnc_validated <- feelnc_data %>%
  filter(lncRNA_transcript %in% validated_lncrna_ids)
write_out(feelnc_validated, out_dir, paste0(prefix, "_feelnc_validated.csv"))

# Recorded for the record: in the historical input every row already had
# isBest == 1 (42,040 of 42,040), so the isBest filter applied downstream in
# R06/R13 is a no-op on this file. See docs/OPEN_ISSUES.md.
n_best <- if ("isBest" %in% names(feelnc_data)) sum(feelnc_data$isBest == 1) else NA_integer_
sw("07_final", data.frame(
  Dataset = c("Final Known Transcripts", "Novel lncRNAs",
              "Validated Known lncRNAs for FEELnc", "Total Validated lncRNAs",
              "Final FEELnc Relations", "FEELnc rows with isBest == 1"),
  Count = c(nrow(meta_annotated_with_exons), length(novel_lncrna_ids),
            length(known_lncrna_ids), length(validated_lncrna_ids),
            nrow(feelnc_validated), n_best)))

# --- Step 8: final combined GTF -----------------------------------------------
message("[R01] Step 8/8  final annotated GTF")
annotated_ids <- unique(meta_annotated_with_exons$ensembl_transcript_id)
novel_ids     <- unique(meta_novel_lncrna_final$transcript_id)
final_ids     <- c(annotated_ids, novel_ids)

all_stranded_gtf <- c(lncrna_gtf, other_rna_gtf)
final_gtf_gr <- all_stranded_gtf[mcols(all_stranded_gtf)$transcript_id %in% final_ids]

ann_meta <- meta_annotated_with_exons %>%
  select(transcript_id = ensembl_transcript_id,
         gene_name = external_gene_name,
         transcript_name = external_transcript_name,
         gene_biotype, transcript_biotype) %>%
  mutate(source_new = "Ensembl",
         gene_name = ifelse(is.na(gene_name), "NA", gene_name),
         transcript_name = ifelse(is.na(transcript_name), "NA", transcript_name))

ref_name_lookup <- meta_annotated_with_exons %>%
  select(ref_id = ensembl_gene_id, ref_name = external_gene_name) %>%
  filter(!is.na(ref_id)) %>% distinct()

novel_meta <- meta_novel_lncrna_final %>%
  distinct(transcript_id, gene_id, class_code, ref_id) %>%
  left_join(ref_name_lookup, by = "ref_id") %>%
  mutate(ref_name_final = ifelse(is.na(ref_name), ref_id, ref_name),
         source_new = "StringTie",
         gene_name = case_when(
           class_code == "u" ~ paste0("u_", gene_id),
           TRUE ~ paste0(class_code, "_novel_gene_vs_", ref_name_final)),
         transcript_name = case_when(
           class_code == "u" ~ paste0("u_", transcript_id),
           TRUE ~ paste0(transcript_id, "_vs_", ref_name_final)),
         gene_biotype = paste0(class_code, "_novel_lncRNA_gene"),
         transcript_biotype = paste0(class_code, "_novel_lncRNA")) %>%
  select(transcript_id, source_new, gene_name, transcript_name,
         gene_biotype, transcript_biotype)

all_meta_to_add <- bind_rows(ann_meta, novel_meta)
mcols_df_merged <- left_join(as.data.frame(mcols(final_gtf_gr)),
                             all_meta_to_add, by = "transcript_id")

final_gtf_gr$source             <- mcols_df_merged$source_new
final_gtf_gr$gene_name          <- mcols_df_merged$gene_name
final_gtf_gr$transcript_name    <- mcols_df_merged$transcript_name
final_gtf_gr$gene_biotype       <- mcols_df_merged$gene_biotype
final_gtf_gr$transcript_biotype <- mcols_df_merged$transcript_biotype
final_gtf_gr$xloc <- NULL; final_gtf_gr$ref_gene_id <- NULL
final_gtf_gr$cmp_ref <- NULL; final_gtf_gr$tss_id <- NULL

final_gtf_path <- file.path(gtf_out_dir, paste0(prefix, "_final_annotated.gtf"))
export(final_gtf_gr, final_gtf_path, format = "gtf")
message("[out] ", basename(final_gtf_path))

sw("08_final_gtf", data.frame(
  Description = c("Annotated transcripts in final GTF", "Novel lncRNAs in final GTF",
                  "Total transcripts in final GTF"),
  Count = c(length(annotated_ids), length(novel_ids), length(final_ids))))

message("[R01] complete -> ", out_dir)
message("[R01] NOTE: the manuscript figures were produced from the Ensembl-114 ",
        "copies of these tables, shipped under data/frozen/. Compare before ",
        "substituting these outputs. See docs/OPEN_ISSUES.md.")
