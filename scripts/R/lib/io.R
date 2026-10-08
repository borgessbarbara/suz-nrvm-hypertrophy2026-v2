# ---------------------------------------------------------------------------
# lib/io.R — input resolution, delimiter sniffing and column-name matching.
#
# Replaces three bad habits inherited from the original scripts:
#
#  1. Hard-coded absolute paths. Eleven raw scripts and six of the final figure
#     scripts begin with
#         base_dir <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/..."
#     (12_..._ALL_RECORDS_FINAL.R repeats it three times, at :8, :21 and :28),
#     18_curated_Fig4A_Fig6C_from_shortlist.R uses "C:/Users/sebau/Downloads",
#     and 36b_Figure3_panelB...R searches c(".", "/mnt/data").
#
#  2. Browser-download filename variants hard-coded as FIRST preference:
#         "Interact_circUP_miDOWN_6h_all_modeM(1).csv"
#         "DE_miRNA_All_with_status_NE6(1).txt"
#         "NE_vs_Ctrl_6h_lncRNA_DE(2).txt"
#     present in 31_, 32f_, 35_, 36b_, 37c_, 39e_ and 24_FINAL_. These pick a
#     stale copy in preference to the clean one. find_input() rejects any
#     candidate matching the "(N)" pattern outright.
#
#  3. Silent degradation. Several originals carry on with NULL/NA when an input
#     is missing. find_input() stops with the search path printed unless the
#     caller passes optional = TRUE.
# ---------------------------------------------------------------------------

.suz_download_artifact_rx <- "\\(\\d+\\)\\.[A-Za-z0-9]+$"

#' Build the ordered list of directories an input may live in.
suz_search_path <- function(in_dir = NULL, extra = character()) {
  p <- c(in_dir, extra, ".")
  p <- p[!is.na(p) & nzchar(p)]
  unique(normalizePath(p, winslash = "/", mustWork = FALSE))
}

#' Locate an input file by basename beneath a search path.
#'
#' @param candidates One or more basenames, in order of preference.
#' @param in_dir     Primary directory to search (recursively).
#' @param extra_dirs Additional directories to search after in_dir.
#' @param optional   If TRUE return NA_character_ instead of stopping.
#' @param what       Human label used in messages.
find_input <- function(candidates,
                       in_dir     = NULL,
                       extra_dirs = character(),
                       optional   = FALSE,
                       what       = NULL) {
  stopifnot(length(candidates) >= 1)
  what <- what %||% candidates[1]

  bad <- grepl(.suz_download_artifact_rx, candidates)
  if (any(bad)) {
    stop("find_input(): candidate(s) look like browser-download duplicates ",
         "and are not allowed: ", paste(candidates[bad], collapse = ", "),
         call. = FALSE)
  }

  dirs <- suz_search_path(in_dir, extra_dirs)
  for (cand in candidates) {
    # exact path first (lets callers pass a full path through unchanged)
    if (file.exists(cand) && !dir.exists(cand)) return(normalizePath(cand, "/"))
    for (d in dirs) {
      if (!dir.exists(d)) next
      direct <- file.path(d, cand)
      if (file.exists(direct)) return(normalizePath(direct, "/"))
      hits <- list.files(d, pattern = paste0("^", .suz_rx_escape(cand), "$"),
                         recursive = TRUE, full.names = TRUE)
      hits <- hits[!grepl(.suz_download_artifact_rx, basename(hits))]
      if (length(hits)) {
        if (length(hits) > 1) {
          hits <- hits[order(file.mtime(hits), decreasing = TRUE)]
          warning("find_input(): ", length(hits), " copies of '", cand,
                  "' under ", d, "; using the newest: ", hits[1],
                  call. = FALSE)
        }
        return(normalizePath(hits[1], "/"))
      }
    }
  }
  if (optional) {
    message("[io] optional input not found, continuing without it: ", what)
    return(NA_character_)
  }
  stop("Required input not found: ", what,
       "\n  tried basenames : ", paste(candidates, collapse = ", "),
       "\n  searched        : ", paste(dirs, collapse = "\n                    "),
       "\nPass the correct directory with --in-dir, or see ",
       "metadata/ for inputs that must be supplied externally.",
       call. = FALSE)
}

.suz_rx_escape <- function(x) gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", x, perl = TRUE)

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

#' Read a delimited table, sniffing the delimiter from the header line.
#'
#' 12_..._ALL_RECORDS_FINAL.R:40-44 sniffs tab-vs-comma; 24_FINAL_... also
#' handles semicolons. This merges both behaviours in one place.
read_auto <- function(path, ...) {
  if (is.na(path)) return(NULL)
  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  hdr <- readLines(con, n = 1L, warn = FALSE)
  if (!length(hdr)) return(NULL)
  n_tab <- lengths(regmatches(hdr, gregexpr("\t", hdr)))
  n_com <- lengths(regmatches(hdr, gregexpr(",",  hdr)))
  n_sem <- lengths(regmatches(hdr, gregexpr(";",  hdr)))
  delim <- c("\t" = n_tab, "," = n_com, ";" = n_sem)
  sep <- names(delim)[which.max(delim)]
  if (max(delim) == 0) sep <- "\t"
  readr::read_delim(path, delim = sep, show_col_types = FALSE,
                    progress = FALSE, locale = readr::locale(encoding = "UTF-8"),
                    ...)
}

#' Resolve a column by exact name first, then by case-insensitive regex.
#'
#' Generalised from MASTER_results_report_builder_v3_8_1.R:157, which maps
#' gene_id across ^ensembl_gene_id$, ^ensgene$ and friends.
choose_col <- function(df, exact = character(), regex = character(),
                       required = TRUE, what = "column") {
  nms <- names(df)
  for (e in exact) if (e %in% nms) return(e)
  for (r in regex) {
    hit <- nms[grepl(r, nms, ignore.case = TRUE)]
    if (length(hit)) return(hit[1])
  }
  if (required) {
    stop("Could not resolve ", what, " in a table with columns: ",
         paste(nms, collapse = ", "), call. = FALSE)
  }
  NA_character_
}

#' Create an output directory, with the Windows MAX_PATH guard that
#' 32f_export_ceRNA_cytoscape...R:28-31 carried.
ensure_out_dir <- function(path) {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(path)) {
      stop("Could not create the output directory:\n  ", path,
           "\nOn Windows this is usually the 260-character MAX_PATH limit. ",
           "Move the analysis to a shorter writable path, for example ",
           "C:/suz, and rerun.", call. = FALSE)
    }
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

#' Write a CSV and log it, returning the path invisibly.
write_out <- function(x, out_dir, filename, ...) {
  p <- file.path(out_dir, filename)
  readr::write_csv(x, p, na = "", ...)
  message("[out] ", filename, "  (", nrow(x), " rows)")
  invisible(p)
}

#' Write a TSV and log it.
write_out_tsv <- function(x, out_dir, filename, ...) {
  p <- file.path(out_dir, filename)
  readr::write_tsv(x, p, na = "", ...)
  message("[out] ", filename, "  (", nrow(x), " rows)")
  invisible(p)
}

#' Record every resolved input in a manifest, so each stage's output directory
#' carries its own provenance.
write_input_manifest <- function(inputs, out_dir,
                                 filename = "00_INPUT_manifest.csv") {
  inputs <- inputs[!vapply(inputs, function(p) is.null(p) || is.na(p), logical(1))]
  if (!length(inputs)) return(invisible(NULL))
  df <- data.frame(
    role       = names(inputs),
    path       = unlist(inputs, use.names = FALSE),
    exists     = file.exists(unlist(inputs, use.names = FALSE)),
    size_bytes = file.size(unlist(inputs, use.names = FALSE)),
    mtime      = format(file.mtime(unlist(inputs, use.names = FALSE)),
                        "%Y-%m-%dT%H:%M:%S"),
    stringsAsFactors = FALSE
  )
  write_out(df, out_dir, filename)
}

#' Read one of the frozen inputs under data/frozen/ (plain or .gz).
read_frozen <- function(name, repo_root = NULL) {
  repo_root <- repo_root %||% suz_repo_root()
  base <- file.path(repo_root, "data", "frozen", name)
  for (p in c(base, paste0(base, ".gz"))) {
    if (file.exists(p)) return(read_auto(p))
  }
  stop("Frozen input not found: ", name,
       "\n  looked in: ", file.path(repo_root, "data", "frozen"),
       "\nSee data/frozen/README.md.", call. = FALSE)
}

#' Best-effort repository root, from this file's own location.
suz_repo_root <- function() {
  this <- tryCatch(normalizePath(sys.frame(1)$ofile, "/"), error = function(e) NA)
  if (is.na(this)) {
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grepl("^--file=", args)])
    this <- if (length(f)) normalizePath(f[1], "/") else NA_character_
  }
  if (is.na(this)) return(normalizePath(".", "/"))
  # scripts/R/<stage>.R  ->  repo root is two levels up
  normalizePath(file.path(dirname(this), "..", ".."), "/", mustWork = FALSE)
}

#' Load the shared libraries from any stage script.
suz_load_lib <- function(which = c("io", "ensembl", "plotting")) {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  here <- if (length(f)) dirname(normalizePath(f[1], "/")) else "scripts/R"
  for (w in which) {
    p <- file.path(here, "lib", paste0(w, ".R"))
    if (file.exists(p)) source(p, local = FALSE)
  }
}
