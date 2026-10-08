# ---------------------------------------------------------------------------
# lib/ensembl.R — the single biomaRt entry point for the whole pipeline.
#
# Every Ensembl query in this repository goes through connect_ensembl(). There
# must be no other useMart()/useEnsembl() call anywhere under scripts/R/;
# scripts/checks/check_biomart_pinning.sh enforces that.
#
# WHY THE HOST IS PINNED
# ----------------------
# The historical scripts connected to the live Ensembl service with no version
# pin at all:
#
#   Babi analysis/2_nrvm_transcriptome_metadata_handling.qmd:533
#       useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl",
#                  version = 114)
#   provenance/circRNA_upstream/00_annotation_biomart.Rmd:42,122
#       useDataset("rnorvegicus_gene_ensembl", useMart("ensembl"))
#
# `version=` is resolved by biomaRt at runtime against www.ensembl.org, and that
# resolution has changed between biomaRt releases; `useMart("ensembl")` is not
# pinned at all. Neither is reproducible. An explicit archive host is.
#
#   https://sep2025.archive.ensembl.org/  ==  Ensembl release 115
#
# Verified against the Ensembl FTP: release-115 files are dated 2025-09-03, and
# release 115 ships Rattus_norvegicus.GRCr8.115.gtf — the exact GTF this project
# locks to (metadata/reference_manifest.tsv).
#
# Do NOT also pass version= together with host=. With an archive host the
# release is already determined, and passing both is redundant at best and
# contradictory at worst.
#
# Releases 114 and 115 are both GRCr8; 113 and earlier are mRatBN7.2. See
# docs/REFERENCE_AUDIT.md.
# ---------------------------------------------------------------------------

SUZ_ENSEMBL_HOST_DEFAULT <- "https://sep2025.archive.ensembl.org"
SUZ_ENSEMBL_DATASET      <- "rnorvegicus_gene_ensembl"

# Confirmed live during the rebuild:
#   curl https://sep2025.archive.ensembl.org/biomart/martservice?type=registry
#   -> <MartURLLocation database="ensembl_mart_115" displayName="Ensembl Genes 115" ...>
# so this host is unambiguously release 115.

#' Strip a trailing slash from a host URL.
#'
#' THIS MATTERS. biomaRt pastes "/biomart/martservice" onto the host, so
#' "https://sep2025.archive.ensembl.org/" becomes
#' ".../biomart//martservice" and the call fails with the entirely unhelpful
#'     Error: Can't retrieve empty body.
#' The same URL without the trailing slash works. Both forms are accepted
#' here so nobody loses an afternoon to it.
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
  SUZ_ENSEMBL_HOST_DEFAULT
}

#' Connect to the pinned Ensembl archive.
#'
#' @param host    Archive URL. Defaults to the release-115 archive.
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
       "metadata/reference_manifest.tsv and docs/REFERENCE_AUDIT.md before ",
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
#' provenance/circRNA_upstream/00_annotation_biomart.Rmd:47 and :129 both call
#' getBM(attributes = ..., values = genes, mart = mart) with no `filters`
#' argument. biomaRt ignores `values` in that case, so those chunks never
#' retrieved what they appear to retrieve, and a rerun would not reproduce
#' their own output. Fail loudly instead of repeating it.
suz_getBM <- function(attributes, filters, values, mart, ...) {
  if (missing(filters) || is.null(filters) || !length(filters) ||
      !any(nzchar(filters))) {
    stop("suz_getBM(): 'filters' is required whenever 'values' is supplied. ",
         "A getBM() call with values= but no filters= silently ignores the ",
         "values (see docs/OPEN_ITEMS.md).", call. = FALSE)
  }
  biomaRt::getBM(attributes = attributes, filters = filters,
                 values = values, mart = mart, ...)
}
