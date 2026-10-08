#!/usr/bin/env Rscript
# =============================================================================
# R04 — circRNA count-matrix assembly, tool merge, DE and export.
#
#   *** RECONSTRUCTED FROM PROVENANCE — NOT THE ORIGINAL SCRIPT.            ***
#   *** Requires confirmation by S. Urquiza / J. P. before the repository   ***
#   *** is frozen. See docs/OPEN_ITEMS.md.                                  ***
#
# WHY THIS EXISTS
# ---------------
# Seven downstream scripts read four families of file that NOTHING in any
# package received produces:
#
#     raw_count_matrix_all_6.txt / raw_count_matrix_all_24.txt
#     06_circRNA_sequences_de6_DE_complete_table.txt
#     07_circRNA_sequences_de24_DE_complete_table.txt
#     circRNA_{UP,DOWN}_{6,24}h_all_records.fa
#
# provenance/circRNA_upstream/03_tool_comparation.Rmd contains the merge
# logic but writes differently-named outputs (all_circrnas.txt,
# all_de6_circrnas.txt, circRNA_sequences_de6_2.fa, and a *normalised* count
# matrix where the consumers want raw counts). The step that renames and
# re-exports is genuinely absent. This stage rebuilds it.
#
# LOGIC PRESERVED VERBATIM FROM THE PROVENANCE
# --------------------------------------------
#  * The asymmetric coordinate fix: +1 is added to the CIRCexplorer start
#    coordinate and NOT to the CIRIquant one
#    (03_tool_comparation.Rmd:16 vs :23). This is a deliberate 0-based vs
#    1-based reconciliation and must survive any future refactor.
#  * The merge key is c("chr", "Start", "End"), and the software tag is one of
#    "circexp3_and_ciriquant" / "circexp3" / "ciriquant".
#  * The BSJ support rule from 01_circexp3.Rmd:91 and 02_ciriquant.Rmd:85 is
#    readNumber > 2, i.e. at least THREE back-splice reads. The methods prose
#    says "two or more". They are not the same rule; --bsj-min makes the
#    choice explicit and the run log records it.
#
# KNOWN BUG IN THE PROVENANCE, NOT REPRODUCED
# -------------------------------------------
# 03_tool_comparation.Rmd:28 does
#     c3_cq <- merge(...)
#     c3_cq <- length(c3_cq$chr)          # <- overwrites the data frame
#     c3_cq$circRNA_ID <- paste0(...)     # <- $<- on an atomic vector: error
# As committed that chunk cannot run past line 30. The intent is obvious (the
# author was noting the intersection size, 46,302) so the count is kept as an
# assertion and the data frame is not clobbered.
#
# RECORDED COUNTS used as soft assertions (warn, do not stop):
#     intersection          46302
#     CIRCexplorer-only      3419
#     CIRIquant-only        23739
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
    make_option("--ciri-as-dir", type = "character", default = "results/15_ciri_as",
                help = "stage 15 output, for reconstructed sequences [default %default]"),
    make_option("--samples", type = "character", default = "metadata/longRNA_samples.tsv",
                help = "[default %default]"),
    make_option("--out-dir", type = "character", default = "results/R04_circRNA_merge_export",
                help = "[default %default]"),
    make_option("--bsj-min", type = "integer", default = 3,
                help = paste("Minimum back-splice reads to retain a circRNA.",
                             "3 reproduces the code (readNumber > 2);",
                             "2 reproduces the methods prose. [default %default]")),
    make_option("--replicates", type = "character", default = "1,2,4",
                help = paste("Replicate indices to include per group. The default",
                             "1,2,4 matches the mRNA/lncRNA DESeq2 and the final",
                             "circRNA volcano script. [default %default]")),
    make_option("--use-combat", action = "store_true", default = FALSE,
                help = paste("Apply ComBat_seq batch correction as the historical",
                             "notebooks did. OFF by default because the final",
                             "manuscript volcano script does not use it.")),
    make_option("--lfc-cut", type = "double", default = 0.5, help = "[default %default]"),
    make_option("--alpha", type = "double", default = 0.05,
                help = "nominal p threshold [default %default]")),
  description = paste(
    "\nR18b: assemble circRNA count matrices from the CIRCexplorer and CIRIquant",
    "\nbranches, merge the two tools, run DE, and export the four file families",
    "\nthat the downstream circRNA stages require.",
    "\n\n*** RECONSTRUCTED - see the header and docs/OPEN_ITEMS.md. ***")))

out_dir <- ensure_out_dir(opt$`out-dir`)
reps <- as.integer(strsplit(opt$replicates, ",")[[1]])

message("=====================================================================")
message("[R04] RECONSTRUCTED STAGE - confirm against the original analysis")
message("[R04]   BSJ minimum reads : ", opt$`bsj-min`,
        if (opt$`bsj-min` == 3) "  (matches the notebook code, readNumber > 2)"
        else if (opt$`bsj-min` == 2) "  (matches the methods prose, >= 2)" else "")
message("[R04]   replicates kept   : ", paste(reps, collapse = ", "))
message("[R04]   ComBat_seq        : ", if (opt$`use-combat`) "ON (historical)" else "OFF (final manuscript)")
message("[R04]   nominal p < ", opt$alpha, ", |log2FC| >= ", opt$`lfc-cut`)
message("=====================================================================")

samples <- read_auto(opt$samples)
stopifnot(all(c("sample_id", "condition", "time", "replicate") %in% names(samples)))
keep <- samples %>% filter(as.integer(replicate) %in% reps)
message("[R04] samples retained: ", nrow(keep), " of ", nrow(samples))

# --- parse a circRNA locus id into chr / start / end --------------------------
# Accepts "chr:start|end" and "chr_start_end"; the canonical form written out
# is always "chr:start|end", matching canonical_locus() in R24.
parse_locus <- function(x) {
  x <- as.character(x)
  tibble(
    chr   = sub(":.*", "", x),
    Start = suppressWarnings(as.numeric(sub(".*:(\\d+)\\|.*", "\\1", x))),
    End   = suppressWarnings(as.numeric(sub(".*\\|(\\d+)", "\\1", x))))
}
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
  out <- map(keep$sample_id, ~ reader(dir, .x)) %>% compact()
  if (!length(out)) {
    stop("No ", label, " per-sample results found under ", dir,
         " for the retained samples.", call. = FALSE)
  }
  message("[R04] ", label, ": ", length(out), " samples parsed")
  bind_rows(out)
}

c3_long <- gather_branch(read_circexplorer, opt$`circexplorer-dir`, "CIRCexplorer")
cq_long <- gather_branch(read_ciriquant,    opt$`ciriquant-dir`,    "CIRIquant")

# --- 2) BSJ support filter ----------------------------------------------------
bsj_filter <- function(df) df %>% filter(!is.na(readNumber), readNumber >= opt$`bsj-min`)
n_before <- nrow(c3_long) + nrow(cq_long)
c3_long <- bsj_filter(c3_long); cq_long <- bsj_filter(cq_long)
message("[R04] BSJ filter (>= ", opt$`bsj-min`, " reads): ",
        n_before, " -> ", nrow(c3_long) + nrow(cq_long), " calls")

# --- 3) tool merge on chr/Start/End ------------------------------------------
c3_loci <- c3_long %>% distinct(chr, Start, End) %>% mutate(pip.x = "c3")
cq_loci <- cq_long %>% distinct(chr, Start, End) %>% mutate(pip.y = "cq")

both   <- inner_join(c3_loci, cq_loci, by = c("chr", "Start", "End")) %>%
            mutate(software = "circexp3_and_ciriquant")
only_3 <- anti_join(c3_loci, cq_loci, by = c("chr", "Start", "End")) %>%
            mutate(software = "circexp3")
only_q <- anti_join(cq_loci, c3_loci, by = c("chr", "Start", "End")) %>%
            mutate(software = "ciriquant")

all_circ <- bind_rows(both, only_3, only_q) %>%
  mutate(circRNA_ID = make_id(chr, Start, End)) %>%
  select(circRNA_ID, chr, Start, End, software)

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

# --- 4) raw count matrices per timepoint -------------------------------------
counts_long <- bind_rows(c3_long, cq_long) %>%
  mutate(circRNA_ID = make_id(chr, Start, End)) %>%
  semi_join(all_circ, by = "circRNA_ID") %>%
  group_by(circRNA_ID, sample_id) %>%
  summarise(readNumber = max(readNumber, na.rm = TRUE), .groups = "drop")

sample_label <- keep %>%
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
  fn <- paste0("raw_count_matrix_all_", sub("h$", "", tp), ".txt")
  readr::write_tsv(mat, file.path(out_dir, fn))
  message("[out] ", fn, "  (", nrow(mat), " loci x ", ncol(mat) - 1, " samples)")
}

# --- 5) DE and the complete annotated tables ---------------------------------
suppressPackageStartupMessages(library(DESeq2))
run_de <- function(tp, idx) {
  fn_in <- file.path(out_dir, paste0("raw_count_matrix_all_", sub("h$", "", tp), ".txt"))
  m <- read_auto(fn_in)
  rn <- m$circRNA_ID
  cm <- as.matrix(m[, -1, drop = FALSE]); rownames(cm) <- rn
  storage.mode(cm) <- "integer"
  cond <- factor(ifelse(grepl("^n", colnames(cm)), "NE", "Ctrl"), levels = c("Ctrl", "NE"))
  batch <- factor(sub(".*_", "", colnames(cm)))

  if (isTRUE(opt$`use-combat`)) {
    if (!requireNamespace("sva", quietly = TRUE)) {
      stop("--use-combat requires the 'sva' package.", call. = FALSE)
    }
    message("[R04] ", tp, ": applying ComBat_seq (historical behaviour)")
    cm <- sva::ComBat_seq(counts = cm, batch = batch, group = cond)
  }

  dds <- DESeqDataSetFromMatrix(cm, data.frame(condition = cond), design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", "NE", "Ctrl"),
                               independentFiltering = TRUE))
  res$circRNA_ID <- rownames(res)

  ann <- all_circ %>% select(circRNA_ID, chr, Start, End, software)
  host <- bind_rows(c3_long, cq_long) %>%
    mutate(circRNA_ID = make_id(chr, Start, End)) %>%
    filter(!is.na(geneName), geneName != "") %>%
    distinct(circRNA_ID, geneName) %>%
    group_by(circRNA_ID) %>%
    summarise(host_gene = paste(sort(unique(geneName)), collapse = ","), .groups = "drop")

  out <- res %>%
    left_join(ann,  by = "circRNA_ID") %>%
    left_join(host, by = "circRNA_ID") %>%
    mutate(timepoint = tp,
           status = case_when(
             !is.na(pvalue) & pvalue < opt$alpha & log2FoldChange >=  opt$`lfc-cut` ~ "Up",
             !is.na(pvalue) & pvalue < opt$alpha & log2FoldChange <= -opt$`lfc-cut` ~ "Down",
             TRUE ~ "Unchanged")) %>%
    select(circRNA_ID, chr, Start, End, host_gene, software, timepoint,
           baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, status)

  fn <- sprintf("%s_circRNA_sequences_de%s_DE_complete_table.txt",
                if (tp == "6h") "06" else "07", sub("h$", "", tp))
  readr::write_tsv(out, file.path(out_dir, fn))
  message("[out] ", fn, "  (", nrow(out), " loci; Up ",
          sum(out$status == "Up"), ", Down ", sum(out$status == "Down"), ")")
  out
}
de6  <- run_de("6h")
de24 <- run_de("24h")

# --- 6) IntaRNA target FASTAs -------------------------------------------------
# CIRI-AS reconstructs the circRNA sequences. If stage 15 output is present we
# use it; otherwise we emit the locus list and stop short of writing FASTA,
# because inventing sequence would be worse than failing.
seq_lookup <- NULL
as_dir <- opt$`ciri-as-dir`
if (dir.exists(as_dir)) {
  as_files <- list.files(as_dir, pattern = "\\.ASout(_jav\\.list)?$",
                         recursive = TRUE, full.names = TRUE)
  if (length(as_files)) {
    message("[R04] CIRI-AS files found: ", length(as_files))
    seq_lookup <- tryCatch({
      map_dfr(as_files, function(f) {
        x <- try(read.delim(f, header = TRUE, stringsAsFactors = FALSE), silent = TRUE)
        if (inherits(x, "try-error") || !ncol(x)) return(NULL)
        sc <- names(x)[grepl("seq", names(x), ignore.case = TRUE)][1]
        ic <- names(x)[grepl("circ|id", names(x), ignore.case = TRUE)][1]
        if (is.na(sc) || is.na(ic)) return(NULL)
        tibble(circRNA_ID = as.character(x[[ic]]), sequence = as.character(x[[sc]]))
      }) %>% filter(!is.na(sequence), nzchar(sequence)) %>% distinct(circRNA_ID, .keep_all = TRUE)
    }, error = function(e) NULL)
  }
}

write_fa <- function(df, path) {
  lines <- as.vector(rbind(paste0(">", df$circRNA_ID), df$sequence))
  writeLines(lines, path)
}

for (tp in c("6h", "24h")) {
  de <- if (tp == "6h") de6 else de24
  for (dir_lab in c("UP", "DOWN")) {
    sel <- de %>% filter(status == if (dir_lab == "UP") "Up" else "Down")
    fn <- sprintf("circRNA_%s_%s_all_records.fa", dir_lab, tp)
    if (!is.null(seq_lookup)) {
      withseq <- sel %>% inner_join(seq_lookup, by = "circRNA_ID")
      write_fa(withseq, file.path(out_dir, fn))
      message("[out] ", fn, "  (", nrow(withseq), " of ", nrow(sel),
              " loci had a reconstructed sequence)")
    } else {
      readr::write_tsv(sel %>% select(circRNA_ID, chr, Start, End, log2FoldChange, status),
                       file.path(out_dir, sub("\\.fa$", "_LOCI_ONLY.tsv", fn)))
      message("[warn] ", fn, " NOT written: no CIRI-AS sequences available. ",
              "Wrote the locus list instead.")
    }
  }
}
if (is.null(seq_lookup)) {
  message("=====================================================================")
  message("[R04] The four circRNA_*_all_records.fa targets were NOT produced,")
  message("[R04] because no CIRI-AS reconstructed sequences were found under")
  message("[R04]   ", as_dir)
  message("[R04] Bash stage 16 (circRNA-miRNA IntaRNA) cannot run without them.")
  message("[R04] Note that stage 15 only covers SQ13-SQ16 in the original, so")
  message("[R04] sequences for the other samples may never have existed.")
  message("[R04] See docs/OPEN_ITEMS.md.")
  message("=====================================================================")
}

writeLines(c(
  "R18b circRNA merge/export - RECONSTRUCTED STAGE",
  "",
  paste0("BSJ minimum reads     : ", opt$`bsj-min`),
  paste0("Replicates included   : ", paste(reps, collapse = ", ")),
  paste0("ComBat_seq            : ", if (opt$`use-combat`) "ON" else "OFF"),
  paste0("Nominal p threshold   : ", opt$alpha),
  paste0("|log2FC| threshold    : ", opt$`lfc-cut`),
  "",
  paste0("Loci - intersection   : ", nrow(both),   "  (recorded 46302)"),
  paste0("Loci - CIRCexplorer   : ", nrow(only_3), "  (recorded 3419)"),
  paste0("Loci - CIRIquant      : ", nrow(only_q), "  (recorded 23739)"),
  "",
  "This stage is a reconstruction of a step whose original script was never",
  "supplied. Confirm the BSJ rule, the replicate set and the ComBat_seq",
  "decision with S. Urquiza and J. P. before these outputs reach the",
  "manuscript. See docs/OPEN_ITEMS.md."),
  file.path(out_dir, "README_RECONSTRUCTION.txt"))

message("[R04] complete -> ", out_dir)
