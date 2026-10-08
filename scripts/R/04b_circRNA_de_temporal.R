#!/usr/bin/env Rscript
# =============================================================================
# R04b — final circRNA differential expression and export.
#
#   Implements stage 5 of 06_NOTES_FOR_COLLABORATOR/repository_merge_plan.md:
#     "circRNA_05_final_DE_temporal.R - final manuscript circRNA DE/temporal
#      classification using the verified sample set and verified BSJ rule."
#
#   *** RECONSTRUCTED FROM PROVENANCE - NOT THE ORIGINAL SCRIPT.            ***
#   *** Requires confirmation by S. Urquiza / J. P. See docs/OPEN_ITEMS.md. ***
#
# EVERY CONTESTED PARAMETER LIVES HERE, AND EVERY ONE IS STILL OPEN
# -----------------------------------------------------------------
#   --bsj-min      Methods say "at least two BSJ counts" (= 2).
#                  01_circexp3.Rmd:91 and 02_ciriquant.Rmd:85 use
#                  readNumber > 2 (= 3). Default 3, matching the CODE.
#                  docs/OPEN_ITEMS.md A1 - BLOCKING.
#   --replicates   Default 1,2,4. Replicate 3 is excluded from every final
#                  analysis, and the manuscript states N = 3 per condition.
#                  The historical notebooks used 4+4. docs/OPEN_ITEMS.md A3.
#   --use-combat   The notebooks applied ComBat_seq; the final volcano script
#                  does not. OFF by default. docs/OPEN_ITEMS.md A2.
#                  NOTE batch maps 1:1 to replicate index, so dropping
#                  replicate 3 removes Batch3 entirely.
#   --alpha /      Nominal p < 0.05 and |log2FC| >= 0.5, per the Methods.
#   --lfc-cut      circRNA DE is explicitly exploratory: nominal p, not FDR.
#
# Because R04 is deterministic and this stage is not, the whole point of the
# split is that you can re-run this under a different BSJ rule, replicate set
# or batch model in seconds, without re-parsing the caller output.
#
# THE TEMPORAL CLASSIFICATION lives in R13 (Figure 3C circRNA), which also
# emits early/sustained/late/reversal_circRNA.csv for R21.
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
    make_option("--merge-dir", type = "character",
                default = "results/R04_circRNA_annotation_merge",
                help = "R04 output [default %default]"),
    make_option("--ciri-as-dir", type = "character", default = "results/15_ciri_as",
                help = "stage 15 output, for reconstructed sequences [default %default]"),
    make_option("--out-dir", type = "character", default = "results/R04b_circRNA_de_temporal",
                help = "[default %default]"),
    make_option("--bsj-min", type = "integer", default = 3,
                help = paste("Minimum back-splice reads. 3 reproduces the code",
                             "(readNumber > 2); 2 reproduces the Methods prose.",
                             "[default %default]")),
    make_option("--replicates", type = "character", default = "1,2,4",
                help = "Replicate indices per group [default %default]"),
    make_option("--use-combat", action = "store_true", default = FALSE,
                help = "Apply ComBat_seq as the historical notebooks did."),
    make_option("--lfc-cut", type = "double", default = 0.5, help = "[default %default]"),
    make_option("--alpha", type = "double", default = 0.05,
                help = "nominal p threshold [default %default]")),
  description = paste(
    "\nR04b: apply the BSJ filter and replicate selection to R04's merged loci,",
    "\nrun circRNA differential expression, and export the four file families",
    "\nthe downstream circRNA stages require.",
    "\n\n*** RECONSTRUCTED - see the header and docs/OPEN_ITEMS.md. ***")))

out_dir <- ensure_out_dir(opt$`out-dir`)
reps <- as.integer(strsplit(opt$replicates, ",")[[1]])

message("=====================================================================")
message("[R04b] RECONSTRUCTED STAGE - every parameter below is an OPEN ITEM")
message("[R04b]   BSJ minimum reads : ", opt$`bsj-min`,
        if (opt$`bsj-min` == 3) "  (matches the notebook CODE, readNumber > 2)"
        else if (opt$`bsj-min` == 2) "  (matches the METHODS prose, >= 2)" else "")
message("[R04b]   replicates kept   : ", paste(reps, collapse = ", "))
message("[R04b]   ComBat_seq        : ", if (opt$`use-combat`) "ON (historical)" else "OFF (final manuscript)")
message("[R04b]   nominal p < ", opt$alpha, ", |log2FC| >= ", opt$`lfc-cut`)
message("=====================================================================")

ann_path <- find_input("all_circrnas.csv", in_dir = opt$`merge-dir`,
                       what = "R04 locus annotation")
all_circ <- read_auto(ann_path)
write_input_manifest(list(locus_annotation = ann_path), out_dir)

# --- 1) BSJ filter + replicate selection, then the raw count matrices --------
suppressPackageStartupMessages(library(DESeq2))

keep_label <- function(lbl) as.integer(sub(".*_", "", lbl)) %in% reps

build_matrix <- function(tp) {
  fn <- find_input(paste0("bsj_count_matrix_unfiltered_", sub("h$", "", tp), ".txt"),
                   in_dir = opt$`merge-dir`,
                   what = paste("R04 unfiltered BSJ matrix", tp))
  m <- read_auto(fn)
  ids <- m$circRNA_ID
  cm <- as.matrix(m[, -1, drop = FALSE]); rownames(cm) <- ids

  # replicate selection
  cm <- cm[, vapply(colnames(cm), keep_label, logical(1)), drop = FALSE]

  # BSJ support: a locus is retained if any kept sample reaches the threshold.
  n_before <- nrow(cm)
  cm <- cm[apply(cm, 1, function(r) any(r >= opt$`bsj-min`, na.rm = TRUE)), , drop = FALSE]
  message("[R04b] ", tp, ": BSJ filter (>= ", opt$`bsj-min`, ") kept ",
          nrow(cm), " of ", n_before, " loci across ", ncol(cm), " samples")

  out <- as_tibble(cm, rownames = "circRNA_ID")
  fn_out <- paste0("raw_count_matrix_all_", sub("h$", "", tp), ".txt")
  readr::write_tsv(out, file.path(out_dir, fn_out))
  message("[out] ", fn_out)
  cm
}

cm6  <- build_matrix("6h")
cm24 <- build_matrix("24h")

# --- 2) DE and the complete annotated tables ---------------------------------
run_de <- function(cm, tp) {
  storage.mode(cm) <- "integer"
  cond  <- factor(ifelse(grepl("^n", colnames(cm)), "NE", "Ctrl"), levels = c("Ctrl", "NE"))
  batch <- factor(sub(".*_", "", colnames(cm)))

  if (isTRUE(opt$`use-combat`)) {
    if (!requireNamespace("sva", quietly = TRUE)) {
      stop("--use-combat requires the 'sva' package.", call. = FALSE)
    }
    message("[R04b] ", tp, ": applying ComBat_seq (historical behaviour)")
    cm <- sva::ComBat_seq(counts = cm, batch = batch, group = cond)
  }

  dds <- DESeqDataSetFromMatrix(cm, data.frame(condition = cond), design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", "NE", "Ctrl"),
                               independentFiltering = TRUE))
  res$circRNA_ID <- rownames(res)

  out <- res %>%
    left_join(all_circ, by = "circRNA_ID") %>%
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
de6  <- run_de(cm6,  "6h")
de24 <- run_de(cm24, "24h")

# --- 3) IntaRNA target FASTAs -------------------------------------------------
# CIRI-AS reconstructs the circRNA sequences. If stage 15 output is present we
# use it; otherwise we emit the locus list and stop short of writing FASTA,
# because inventing sequence would be worse than failing.
seq_lookup <- NULL
as_dir <- opt$`ciri-as-dir`
if (dir.exists(as_dir)) {
  as_files <- list.files(as_dir, pattern = "\\.ASout(_jav\\.list)?$",
                         recursive = TRUE, full.names = TRUE)
  if (length(as_files)) {
    message("[R04b] CIRI-AS files found: ", length(as_files))
    seq_lookup <- tryCatch({
      map_dfr(as_files, function(f) {
        x <- try(read.delim(f, header = TRUE, stringsAsFactors = FALSE), silent = TRUE)
        if (inherits(x, "try-error") || !ncol(x)) return(NULL)
        sc <- names(x)[grepl("seq", names(x), ignore.case = TRUE)][1]
        ic <- names(x)[grepl("circ|id", names(x), ignore.case = TRUE)][1]
        if (is.na(sc) || is.na(ic)) return(NULL)
        tibble(circRNA_ID = as.character(x[[ic]]), sequence = as.character(x[[sc]]))
      }) %>% filter(!is.na(sequence), nzchar(sequence)) %>%
        distinct(circRNA_ID, .keep_all = TRUE)
    }, error = function(e) NULL)
  }
}

write_fa <- function(df, path) {
  writeLines(as.vector(rbind(paste0(">", df$circRNA_ID), df$sequence)), path)
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
  message("[R04b] The four circRNA_*_all_records.fa targets were NOT produced,")
  message("[R04b] because no CIRI-AS reconstructed sequences were found under")
  message("[R04b]   ", as_dir)
  message("[R04b] Bash stage 16 (circRNA-miRNA IntaRNA) cannot run without them.")
  message("[R04b] Note that stage 15 only covers SQ13-SQ16 in the original, so")
  message("[R04b] sequences for the other samples may never have existed.")
  message("[R04b] See docs/OPEN_ITEMS.md A5 and A7.")
  message("=====================================================================")
}

writeLines(c(
  "R04b circRNA DE and export - RECONSTRUCTED STAGE",
  "",
  "Implements stage 5 of repository_merge_plan.md.",
  "",
  paste0("BSJ minimum reads     : ", opt$`bsj-min`,
         if (opt$`bsj-min` == 3) "   (matches the notebook code)"
         else if (opt$`bsj-min` == 2) "   (matches the Methods prose)" else ""),
  paste0("Replicates included   : ", paste(reps, collapse = ", ")),
  paste0("ComBat_seq            : ", if (opt$`use-combat`) "ON" else "OFF"),
  paste0("Nominal p threshold   : ", opt$alpha),
  paste0("|log2FC| threshold    : ", opt$`lfc-cut`),
  paste0("Loci after BSJ filter : 6h ", nrow(cm6), ", 24h ", nrow(cm24)),
  "",
  "Every parameter above is an unresolved open item. Confirm the BSJ rule, the",
  "replicate set and the ComBat_seq decision with S. Urquiza and J. P. before",
  "these outputs reach the manuscript. See docs/OPEN_ITEMS.md."),
  file.path(out_dir, "README_RECONSTRUCTION.txt"))

message("[R04b] complete -> ", out_dir)
