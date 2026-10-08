#!/usr/bin/env Rscript
# =============================================================================
# R04 - circRNA annotation and CIRCexplorer/CIRIquant locus merge.
#
# Parses both callers per sample, reconciles their coordinate conventions,
# merges loci, and writes an UNFILTERED back-splice-junction count matrix plus
# a locus annotation table. Deterministic: no BSJ threshold, replicate
# selection or batch correction happens here. Those all live in R04b.
#
# Outputs: all_circrnas.csv, bsj_count_matrix_unfiltered_{6,24}.txt
#
# Usage:   Rscript scripts/R/04_circRNA_annotation_merge.R [--circexplorer-dir DIR] [--ciriquant-dir DIR] [--out-dir DIR]
#          Rscript scripts/R/04_circRNA_annotation_merge.R --help
#
# RECONSTRUCTED: no original script for this step was supplied. Rebuilt from
# archive/final_package_v3_original/02_CIRCRNA_UPSTREAM_RAW_PROVENANCE/03_tool_comparation.Rmd
# See docs/CHANGES_FROM_ORIGINALS.md and docs/OPEN_ISSUES.md before relying on it.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(stringr)
  library(tibble); library(readr); library(purrr)
})

.here <- dirname(normalizePath(sub("^--file=", "",
          commandArgs(trailingOnly = FALSE)[grepl("^--file=",
          commandArgs(trailingOnly = FALSE))][1]), winslash = "/"))
source(file.path(.here, "lib", "io.R"))

opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--circexplorer-dir", type = "character", default = "results/13_circexplorer",
                help = "stage 13 output [default %default]"),
    make_option("--ciriquant-dir", type = "character", default = "results/14_ciriquant",
                help = "stage 14 output [default %default]"),
    make_option("--samples", type = "character", default = "metadata/longRNA_samples.tsv",
                help = "[default %default]"),
    make_option("--out-dir", type = "character", default = "results/R04_circRNA_annotation_merge",
                help = "[default %default]")),
  description = paste(
    "\nR04: parse the CIRCexplorer and CIRIquant per-sample calls, reconcile",
    "\ncoordinates, merge loci, and emit an unfiltered BSJ count matrix plus a",
    "\nlocus annotation table.",
    "\n\nNo BSJ, replicate or batch-correction choices are made here - see R04b.",
    "\n\n*** RECONSTRUCTED - see the header and docs/OPEN_ISSUES.md. ***")))

out_dir <- ensure_out_dir(opt$`out-dir`)

message("=====================================================================")
message("[R04] RECONSTRUCTED STAGE - confirm against the original analysis")
message("[R04] Deterministic merge only. BSJ filter, replicate selection and")
message("[R04] ComBat_seq are applied downstream in R04b.")
message("=====================================================================")

samples <- read_auto(opt$samples)
stopifnot(all(c("sample_id", "condition", "time", "replicate") %in% names(samples)))
message("[R04] samples in sheet: ", nrow(samples), " (all parsed; none excluded here)")

make_id <- function(chr, start, end) paste0(chr, ":", start, "|", end)

# --- 1) per-sample calls from each branch -------------------------------------
read_circexplorer <- function(dir, sample_id) {
  f <- file.path(dir, paste0(sample_id, "_circexplorer_results"),
                 paste0(sample_id, "_CIRCRNA_circularRNA_known.txt"))
  if (!file.exists(f)) return(NULL)
  # CIRCexplorer2 annotate output has no header; 21 columns.
  cn <- c("chrom","start","end","name","score","strand","thickStart","thickEnd",
          "itemRgb","exonCount","exonSizes","exonOffsets","readNumber","circType",
          "geneName","isoformName","index","flankIntron","FPBcirc","FBPlinear","CIRCscore")
  df <- read.delim(f, header = FALSE, stringsAsFactors = FALSE)
  names(df)[seq_len(min(ncol(df), length(cn)))] <- cn[seq_len(min(ncol(df), length(cn)))]
  df %>%
    transmute(sample_id = sample_id,
              chr = as.character(chrom),
              # 03_tool_comparation.Rmd:16 - CIRCexplorer start gets +1.
              Start = as.numeric(start) + 1,
              End = as.numeric(end),
              readNumber = suppressWarnings(as.numeric(readNumber)),
              geneName = as.character(geneName),
              pip = "c3")
}

read_ciriquant <- function(dir, sample_id) {
  cand <- c(file.path(dir, sample_id, paste0(sample_id, ".gtf")),
            file.path(dir, sample_id, "circ", paste0(sample_id, ".ciri")))
  f <- cand[file.exists(cand)][1]
  if (is.na(f)) return(NULL)
  if (grepl("\\.gtf$", f)) {
    gt <- read.delim(f, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
    attr9 <- gt[[9]]
    df <- tibble(
      chr   = as.character(gt[[1]]),
      Start = as.numeric(gt[[4]]),
      End   = as.numeric(gt[[5]]),
      readNumber = suppressWarnings(as.numeric(
        str_match(attr9, 'bsj\\s+"?([0-9.]+)"?')[, 2])),
      geneName = str_match(attr9, 'gene_name\\s+"([^"]+)"')[, 2])
  } else {
    ci <- read.delim(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    df <- tibble(
      chr   = as.character(ci[[2]]),
      Start = as.numeric(ci[[3]]),
      End   = as.numeric(ci[[4]]),
      readNumber = suppressWarnings(as.numeric(ci[[5]])),
      geneName = if ("gene_id" %in% names(ci)) as.character(ci$gene_id) else NA_character_)
  }
  # NOTE: no +1 here. The asymmetry is deliberate (03_tool_comparation.Rmd:23).
  df %>% mutate(sample_id = sample_id, pip = "cq") %>%
    select(sample_id, chr, Start, End, readNumber, geneName, pip)
}

gather_branch <- function(reader, dir, label) {
  if (!dir.exists(dir)) {
    stop("Required ", label, " output directory not found: ", dir,
         "\nRun the corresponding bash stage first.", call. = FALSE)
  }
  out <- map(samples$sample_id, ~ reader(dir, .x)) %>% compact()
  if (!length(out)) {
    stop("No ", label, " per-sample results found under ", dir, ".", call. = FALSE)
  }
  message("[R04] ", label, ": ", length(out), " samples parsed")
  bind_rows(out)
}

c3_long <- gather_branch(read_circexplorer, opt$`circexplorer-dir`, "CIRCexplorer")
cq_long <- gather_branch(read_ciriquant,    opt$`ciriquant-dir`,    "CIRIquant")
write_input_manifest(list(circexplorer_dir = opt$`circexplorer-dir`,
                          ciriquant_dir = opt$`ciriquant-dir`,
                          samples = opt$samples), out_dir)

# --- 2) tool merge on chr/Start/End ------------------------------------------
c3_loci <- c3_long %>% distinct(chr, Start, End) %>% mutate(pip.x = "c3")
cq_loci <- cq_long %>% distinct(chr, Start, End) %>% mutate(pip.y = "cq")

both   <- inner_join(c3_loci, cq_loci, by = c("chr", "Start", "End")) %>%
            mutate(software = "circexp3_and_ciriquant")
only_3 <- anti_join(c3_loci, cq_loci, by = c("chr", "Start", "End")) %>%
            mutate(software = "circexp3")
only_q <- anti_join(cq_loci, c3_loci, by = c("chr", "Start", "End")) %>%
            mutate(software = "ciriquant")

host <- bind_rows(c3_long, cq_long) %>%
  mutate(circRNA_ID = make_id(chr, Start, End)) %>%
  filter(!is.na(geneName), geneName != "") %>%
  distinct(circRNA_ID, geneName) %>%
  group_by(circRNA_ID) %>%
  summarise(host_gene = paste(sort(unique(geneName)), collapse = ","), .groups = "drop")

all_circ <- bind_rows(both, only_3, only_q) %>%
  mutate(circRNA_ID = make_id(chr, Start, End)) %>%
  left_join(host, by = "circRNA_ID") %>%
  select(circRNA_ID, chr, Start, End, host_gene, software) %>%
  arrange(circRNA_ID)

recorded <- c(intersection = 46302, circexplorer_only = 3419, ciriquant_only = 23739)
observed <- c(intersection = nrow(both), circexplorer_only = nrow(only_3),
              ciriquant_only = nrow(only_q))
for (k in names(recorded)) {
  if (observed[[k]] != recorded[[k]]) {
    warning(sprintf(paste("[R04] %s = %d, but 03_tool_comparation.Rmd records %d.",
                          "The reconstruction diverges from the historical run."),
                    k, observed[[k]], recorded[[k]]), call. = FALSE)
  }
}
message("[R04] loci: intersection ", nrow(both), " | CIRCexplorer-only ",
        nrow(only_3), " | CIRIquant-only ", nrow(only_q),
        "  (recorded: 46302 / 3419 / 23739)")
write_out(all_circ, out_dir, "all_circrnas.csv")

# --- 3) UNFILTERED per-locus, per-sample BSJ counts ---------------------------
# No BSJ threshold is applied. R04b filters.
counts_long <- bind_rows(c3_long, cq_long) %>%
  mutate(circRNA_ID = make_id(chr, Start, End)) %>%
  semi_join(all_circ, by = "circRNA_ID") %>%
  group_by(circRNA_ID, sample_id) %>%
  summarise(readNumber = max(readNumber, na.rm = TRUE), .groups = "drop")

sample_label <- samples %>%
  transmute(sample_id,
            label = paste0(ifelse(treatment == "NE", "n", "c"),
                           sub("h$", "", time), "_", replicate),
            time)

for (tp in c("6h", "24h")) {
  sl <- sample_label %>% filter(time == tp)
  mat <- counts_long %>%
    filter(sample_id %in% sl$sample_id) %>%
    left_join(sl, by = "sample_id") %>%
    select(circRNA_ID, label, readNumber) %>%
    pivot_wider(names_from = label, values_from = readNumber, values_fill = 0) %>%
    arrange(circRNA_ID)
  fn <- paste0("bsj_count_matrix_unfiltered_", sub("h$", "", tp), ".txt")
  readr::write_tsv(mat, file.path(out_dir, fn))
  message("[out] ", fn, "  (", nrow(mat), " loci x ", ncol(mat) - 1, " samples)")
}

writeLines(c(
  "R04 circRNA annotation and locus merge - RECONSTRUCTED STAGE",
  "",
  "Implements stage 4 of repository_merge_plan.md.",
  "Deterministic: no BSJ threshold, replicate selection or batch correction.",
  "",
  paste0("Loci - intersection   : ", nrow(both),   "  (recorded 46302)"),
  paste0("Loci - CIRCexplorer   : ", nrow(only_3), "  (recorded 3419)"),
  paste0("Loci - CIRIquant      : ", nrow(only_q), "  (recorded 23739)"),
  paste0("Loci - total          : ", nrow(all_circ)),
  "",
  "The CIRCexplorer start coordinate carries the +1 correction; the CIRIquant",
  "one does not. That asymmetry is deliberate and comes from",
  "03_tool_comparation.Rmd:16 vs :23.",
  "",
  "Next: scripts/R/04b_circRNA_de_temporal.R"),
  file.path(out_dir, "README_RECONSTRUCTION.txt"))

message("[R04] complete -> ", out_dir)
message("[R04] Next: Rscript scripts/R/04b_circRNA_de_temporal.R")
