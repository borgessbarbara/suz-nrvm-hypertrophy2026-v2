#!/usr/bin/env Rscript
# =============================================================================
# check_io_graph.R — every file a stage reads must have a known source.
#
# The single biggest defect in the material received was orphan inputs: files
# that several scripts read and that nothing anywhere produced. Four families
# of them (raw_count_matrix_all_*, 06/07_circRNA_*_DE_complete_table,
# circRNA_*_all_records.fa, 12_STRICT_DE_coherent_circRNA_triplets.csv) were
# only discovered by reading all 47 scripts side by side.
#
# This check makes that class of defect impossible to reintroduce silently.
# For every filename a stage asks for, it must be one of:
#   * produced by an earlier stage      (metadata/stage_io.tsv)
#   * a frozen input                    (data/frozen/)
#   * declared external and not produced by this pipeline
#                                       (metadata/circRNA_input_manifest.tsv,
#                                        metadata/reference_manifest.tsv)
# Anything else is reported as an orphan.
# =============================================================================

suppressPackageStartupMessages({ library(readr); library(dplyr); library(stringr) })

root <- normalizePath(file.path(dirname(sub("^--file=", "",
          commandArgs(trailingOnly = FALSE)[grepl("^--file=",
          commandArgs(trailingOnly = FALSE))][1])), "..", ".."), winslash = "/")
setwd(root)

cat("== stage I/O graph check ==\n")

stage_files <- sort(list.files("scripts/R", pattern = "^[0-9]+b?_.*\\.R$",
                               full.names = TRUE))
stage_of <- function(p) sub("^([0-9]+b?)_.*$", "\\1", basename(p))

# R stages rank by their number ("18b" sorts between 18 and 19). Bash stages
# are given rank -1: they interleave with the R stages (bash 16 runs after R13
# and R18b), so an ordering comparison against them would be meaningless. Their
# outputs still count as a valid source, which is what matters here.
stage_rank <- function(s) {
  ifelse(grepl("^bash", s), -1, suppressWarnings(as.numeric(sub("b$", ".5", s))))
}

# --- filenames each stage asks for -------------------------------------------
# Only filenames in a READ context count. A naive scan of every quoted string
# would also pick up the stage's own output names and report them as orphans.
#
# A read context starts on a line calling one of the resolver/reader functions
# below and continues until the parentheses balance, so multi-line candidate
# vectors like
#     lnc_file <- find_input(
#       candidates = c("a.csv",
#                      "b.csv"))
# are captured in full.
READ_FNS <- paste0(
  "\\b(find_input|find_file|find_first|find_optional_file|find_exact_file|",
  "find_all_interaction_files|read_auto|read_flex|read\\.csv|read_csv|",
  "read\\.delim|read_tsv|read\\.table|read_frozen|readRDS|",
  "read_excel|import)\\s*\\(")

FILE_RX <- '["\']([A-Za-z0-9_.()\\-]+\\.(?:csv|txt|tsv|gtf|fa|fasta|xlsx))["\']'

collect_inputs <- function(path) {
  txt <- readLines(path, warn = FALSE)
  txt <- sub("#.*$", "", txt)          # strip trailing comments
  found <- character(0)
  depth <- 0L
  for (ln in txt) {
    starts <- grepl(READ_FNS, ln, perl = TRUE)
    if (depth > 0L || starts) {
      m <- str_match_all(ln, FILE_RX)[[1]]
      if (nrow(m)) found <- c(found, m[, 2])
      opens  <- lengths(regmatches(ln, gregexpr("\\(", ln)))
      closes <- lengths(regmatches(ln, gregexpr("\\)", ln)))
      depth <- max(0L, depth + opens - closes)
      if (!starts && depth == 0L) next
      if (starts && depth <= 0L) depth <- 0L
    }
  }
  unique(found)
}

# --- declared outputs ---------------------------------------------------------
io_path <- "metadata/stage_io.tsv"
if (!file.exists(io_path)) {
  cat("FAIL: ", io_path, " not found; cannot resolve stage outputs.\n", sep = "")
  quit(status = 1)
}
io <- suppressWarnings(read_tsv(io_path, show_col_types = FALSE, progress = FALSE))
stopifnot(all(c("stage", "direction", "filename") %in% names(io)))

outputs <- io %>% filter(direction == "out") %>%
  transmute(filename, stage = as.character(stage), rank = stage_rank(stage))

frozen <- basename(list.files("data/frozen", pattern = "\\.(csv|txt|tsv)(\\.gz)?$"))
frozen <- unique(sub("\\.gz$", "", frozen))

external <- character(0)
for (f in c("metadata/circRNA_input_manifest.tsv", "metadata/reference_manifest.tsv",
             "metadata/external_inputs.tsv")) {
  if (!file.exists(f)) next
  x <- suppressWarnings(read_tsv(f, show_col_types = FALSE, progress = FALSE))
  cols <- intersect(c("expected_filename_or_pattern", "final_file_name",
                      "resource_name", "filename"), names(x))
  for (cc in cols) external <- c(external, as.character(x[[cc]]))
}
external <- unique(na.omit(external))
# brace-expand the manifest patterns, e.g. raw_count_matrix_all_{6,24}.txt
expand_braces <- function(s) {
  out <- s
  while (any(grepl("\\{[^}]*\\}", out))) {
    out <- unlist(lapply(out, function(z) {
      m <- regmatches(z, regexpr("\\{[^}]*\\}", z))
      if (!length(m)) return(z)
      alts <- strsplit(gsub("[{}]", "", m), ",")[[1]]
      vapply(alts, function(a) sub("\\{[^}]*\\}", a, z), character(1))
    }), use.names = FALSE)
  }
  out
}
external <- unique(unlist(lapply(external, expand_braces)))

report <- list(); orphans <- 0L; ordering <- 0L

for (p in stage_files) {
  s <- stage_of(p); r <- stage_rank(s)
  ins <- collect_inputs(p)
  own <- outputs$filename[outputs$stage == s]
  ins <- setdiff(ins, own)                     # a stage reading its own output is fine
  if (!length(ins)) next
  for (f in ins) {
    src <- outputs %>% filter(filename == f)
    if (nrow(src)) {
      if (min(src$rank) > r) {
        cat(sprintf("ORDER  R%-4s reads %-58s produced later by R%s\n",
                    s, f, src$stage[which.min(src$rank)]))
        ordering <- ordering + 1L
      }
      next
    }
    if (f %in% frozen)   next
    if (f %in% external) next
    cat(sprintf("ORPHAN R%-4s reads %-58s no declared source\n", s, f))
    orphans <- orphans + 1L
  }
}

cat("\n")
cat("  stages scanned        : ", length(stage_files), "\n", sep = "")
cat("  declared outputs      : ", nrow(outputs), "\n", sep = "")
cat("  frozen inputs         : ", length(frozen), "\n", sep = "")
cat("  declared external     : ", length(external), "\n", sep = "")
cat("  orphan inputs         : ", orphans, "\n", sep = "")
cat("  ordering violations   : ", ordering, "\n", sep = "")

if (orphans > 0 || ordering > 0) {
  cat("\nFAILED - declare the source in metadata/stage_io.tsv, or list the file as\n")
  cat("external in metadata/circRNA_input_manifest.tsv, or fix the stage order.\n")
  quit(status = 1)
}
cat("\nPASS\n")
