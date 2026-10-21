# ---------------------------------------------------------------------------
# lib/ensembl.R - the single biomaRt entry point.
#
# Every Ensembl query goes through connect_ensembl(). There must be no other
# useMart()/useEnsembl() call under scripts/R/;
# scripts/checks/check_biomart_pinning.sh enforces that.
#
# The release is set in metadata/reference_release.tsv and drives both the
# biomaRt archive host and the GTF that bash stage 00 downloads, so they
# cannot drift apart. Override for one run with ENSEMBL_RELEASE=114.
#
#   release 114 -> https://may2025.archive.ensembl.org
#   release 115 -> https://sep2025.archive.ensembl.org
#   release 116 -> https://jun2026.archive.ensembl.org   (active)
#
# Exports: connect_ensembl(), check_ensembl_connection(), suz_getBM(),
#          suz_ensembl_release(), suz_ensembl_host().
# ---------------------------------------------------------------------------

SUZ_ENSEMBL_ARCHIVES <- c("114" = "https://may2025.archive.ensembl.org",
                          "115" = "https://sep2025.archive.ensembl.org",
                          "116" = "https://jun2026.archive.ensembl.org")

#' Read the active Ensembl release from metadata/reference_release.tsv.
suz_ensembl_release <- function(repo_root = NULL) {
  env <- Sys.getenv("ENSEMBL_RELEASE", unset = "")
  if (nzchar(env)) return(trimws(env))
  repo_root <- repo_root %||% tryCatch(suz_repo_root(), error = function(e) ".")
  f <- file.path(repo_root, "metadata", "reference_release.tsv")
  if (file.exists(f)) {
    x <- try(utils::read.delim(f, stringsAsFactors = FALSE), silent = TRUE)
    if (!inherits(x, "try-error") && "active" %in% names(x) && "release" %in% names(x)) {
      hit <- x$release[tolower(as.character(x$active)) %in% c("yes", "true", "1")]
      if (length(hit)) return(trimws(as.character(hit[1])))
    }
  }
  "116"   # fall back to the active lock if the table is missing or unreadable
}

SUZ_ENSEMBL_HOST_DEFAULT <- NULL   # resolved at call time from the release
SUZ_ENSEMBL_DATASET      <- "rnorvegicus_gene_ensembl"

# Each archive was confirmed live against its own registry, so no host is
# mapped to a release on the strength of its name:
#   curl <host>/biomart/martservice?type=registry
#   may2025 -> database="ensembl_mart_114"  displayName="Ensembl Genes 114"
#   sep2025 -> database="ensembl_mart_115"  displayName="Ensembl Genes 115"
#   jun2026 -> database="ensembl_mart_116"  displayName="Ensembl Genes 116"
# All three serve GRCr8, so switching release changes the annotation only.

#' Strip a trailing slash from a host URL.
#'
suz_normalize_host <- function(host) {
  if (is.null(host) || !nzchar(host)) return(host)
  sub("/+$", "", trimws(host))
}

#' Resolve the Ensembl host to use.
#'
#' Precedence: explicit argument > R option `suz.ensembl_host` >
#' environment variable `ENSEMBL_HOST` > the pinned default.
suz_ensembl_host <- function(host = NULL) {
  if (!is.null(host) && nzchar(host)) return(suz_normalize_host(host))
  opt <- getOption("suz.ensembl_host", NULL)
  if (!is.null(opt) && nzchar(opt)) return(suz_normalize_host(opt))
  env <- Sys.getenv("ENSEMBL_HOST", unset = "")
  if (nzchar(env)) return(suz_normalize_host(env))

  rel <- suz_ensembl_release()
  if (!rel %in% names(SUZ_ENSEMBL_ARCHIVES)) {
    stop("Ensembl release '", rel, "' has no archive mapped in ",
         "SUZ_ENSEMBL_ARCHIVES (scripts/R/lib/ensembl.R). Known: ",
         paste(names(SUZ_ENSEMBL_ARCHIVES), collapse = ", "),
         ".\nAdd the archive host there and in metadata/reference_release.tsv.",
         call. = FALSE)
  }
  message("[ensembl] release : ", rel, "  (metadata/reference_release.tsv)")
  suz_normalize_host(unname(SUZ_ENSEMBL_ARCHIVES[[rel]]))
}

#' Connect to the pinned Ensembl archive.
#'
#' @param host    Archive URL. Defaults to the active release's archive.
#' @param dataset biomaRt dataset; defaults to rat.
#' @param retries Number of attempts; archive mirrors occasionally time out.
#' @return A biomaRt Mart object.
connect_ensembl <- function(host    = NULL,
                            dataset = SUZ_ENSEMBL_DATASET,
                            retries = 3L) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required. Install it with:\n",
         "  BiocManager::install('biomaRt')", call. = FALSE)
  }
  host <- suz_ensembl_host(host)
  message("[ensembl] host    : ", host)
  message("[ensembl] dataset : ", dataset)

  last_err <- NULL
  for (attempt in seq_len(retries)) {
    mart <- tryCatch(
      biomaRt::useEnsembl(biomart = "genes", dataset = dataset, host = host),
      error = function(e) { last_err <<- e; NULL }
    )
    if (!is.null(mart)) {
      if (attempt > 1L) message("[ensembl] connected on attempt ", attempt)
      return(mart)
    }
    message("[ensembl] attempt ", attempt, " failed: ",
            conditionMessage(last_err))
    if (attempt < retries) Sys.sleep(5 * attempt)
  }
  stop("Could not reach the pinned Ensembl archive after ", retries,
       " attempts.\n  host: ", host,
       "\n  last error: ", conditionMessage(last_err),
       "\nIf the archive has been retired, record the replacement host in ",
       "metadata/reference_manifest.tsv and docs/CHANGES_FROM_ORIGINALS.md before ",
       "changing it here — the annotation release is part of the methods.",
       call. = FALSE)
}

#' Report which release the pinned host actually serves.
#'
#' Used by the --dry-run-biomart smoke test. Returns invisibly a list with the
#' host, dataset and the dataset version string biomaRt reports.
check_ensembl_connection <- function(host = NULL, dataset = SUZ_ENSEMBL_DATASET) {
  mart <- connect_ensembl(host = host, dataset = dataset)
  ds <- tryCatch(biomaRt::listDatasets(mart), error = function(e) NULL)
  version <- NA_character_
  if (!is.null(ds) && nrow(ds)) {
    hit <- ds[ds$dataset == dataset, , drop = FALSE]
    if (nrow(hit)) version <- as.character(hit$version[1])
  }
  message("[ensembl] dataset version reported by the server: ",
          if (is.na(version)) "<unavailable>" else version)
  invisible(list(host = suz_ensembl_host(host), dataset = dataset,
                 version = version, mart = mart))
}

#' getBM wrapper that refuses the historical `values` without `filters` bug.
#'
#' archive/final_package_v3_original/02_CIRCRNA_UPSTREAM_RAW_PROVENANCE/00_annotation_biomart.Rmd:47 and :129 both call
#' getBM(attributes = ..., values = genes, mart = mart) with no `filters`
#' argument. biomaRt ignores `values` in that case, so those chunks never
#' retrieved what they appear to retrieve, and a rerun would not reproduce
#' their own output. Fail loudly instead of repeating it.
suz_getBM <- function(attributes, filters, values, mart, ...) {
  if (missing(filters) || is.null(filters) || !length(filters) ||
      !any(nzchar(filters))) {
    stop("suz_getBM(): 'filters' is required whenever 'values' is supplied. ",
         "A getBM() call with values= but no filters= silently ignores the ",
         "values (see docs/OPEN_ISSUES.md).", call. = FALSE)
  }
  biomaRt::getBM(attributes = attributes, filters = filters,
                 values = values, mart = mart, ...)
}
