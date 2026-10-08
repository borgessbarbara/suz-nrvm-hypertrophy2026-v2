#!/usr/bin/env Rscript
# =============================================================================
# R25 — Validate and repair the Cytoscape import keys produced by R22
#
# PORTED, NOT REWRITTEN.
# Source: archive/.../33_validate_fix_cytoscape_import_keys.R (121 lines)
# Produced by scripts/dev/port_from_archive.sh, which applies exactly five
# classes of change and passes everything else through unmodified:
#   1. this preamble (optparse CLI, lib/ sourcing, base_dir/out_dir/dpi)
#   2. hard-coded C:/Users/sebau/... base_dir and out_dir commented out
#   3. base = "C:/Users/sebau/..." default arguments -> base = base_dir
#   4. browser-download duplicates "name(1).csv" removed from candidate lists
#   5. the /mnt/data container search path removed
#
# Diff this file against the archived original to audit the port:
#   diff archive/.../33_validate_fix_cytoscape_import_keys.R scripts/R/23_cytoscape_key_validation.R
# =============================================================================

suppressPackageStartupMessages(library(optparse))

.here <- dirname(normalizePath(sub("^--file=", "",
          commandArgs(trailingOnly = FALSE)[grepl("^--file=",
          commandArgs(trailingOnly = FALSE))][1]), winslash = "/"))
source(file.path(.here, "lib", "io.R"))
source(file.path(.here, "lib", "plotting.R"))

.opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--in-dir", type = "character", default = "results",
                help = "root directory searched for inputs [default %default]"),
    make_option("--out-dir", type = "character", default = "results/R25_cytoscape_key_validation",
                help = "output directory [default %default]"),
    make_option("--dpi", type = "integer", default = 600,
                help = "raster figure resolution [default %default]")),
  description = "\nR23: Validate and repair the Cytoscape import keys produced by R22\nPorted from 33_validate_fix_cytoscape_import_keys.R."))

base_dir <- normalizePath(.opt$`in-dir`, winslash = "/", mustWork = FALSE)
out_dir  <- ensure_out_dir(.opt$`out-dir`)
.suz_out_root <- out_dir   ## [ported] preserved: the body below rebinds out_dir
png_dpi  <- .opt$dpi
dpi_png  <- .opt$dpi
message("[R25] in-dir : ", base_dir)
message("[R25] out-dir: ", out_dir)

# ----------------------------- ported body ----------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# ============================================================
# Validate and repair Cytoscape node/edge identifiers
#
# Expected inputs (searched recursively):
#   ALL_MFE20_nodes.csv
#   ALL_MFE20_edges.csv
#
# Outputs:
#   Cytoscape_import_safe/ALL_MFE20_edges_CYTOSCAPE_SAFE.tsv
#   Cytoscape_import_safe/ALL_MFE20_nodes_CYTOSCAPE_SAFE.tsv
#   Cytoscape_import_safe/ALL_MFE20_edges_MINIMAL.tsv
#   Cytoscape_import_safe/Cytoscape_ID_match_report.txt
#   Cytoscape_import_safe/missing_edge_endpoints_from_nodes.tsv (only if needed)
#   Cytoscape_import_safe/node_ids_not_used_by_edges.tsv (only if needed)
# ============================================================

find_latest <- function(pattern) {
  hits <- list.files(".", pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) stop("Could not find a file matching: ", pattern)
  hits[order(file.info(hits)$mtime, decreasing = TRUE)][1]
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\u00A0", " ")
  x <- str_replace(x, "^\\ufeff", "")
  x <- str_trim(x)
  x <- str_replace(x, '^"(.*)"$', '\\1')
  x <- str_replace(x, "^'(.*)'$", "\\1")
  x <- str_trim(x)
  x
}

nodes_file <- find_latest("^ALL_MFE20_nodes\\.csv$")
edges_file <- find_latest("^ALL_MFE20_edges\\.csv$")

nodes <- read_csv(nodes_file, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
edges <- read_csv(edges_file, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))

required_nodes <- "node_id"
required_edges <- c("source", "target")
if (!all(required_nodes %in% names(nodes))) stop("Node file does not contain column: node_id")
if (!all(required_edges %in% names(edges))) stop("Edge file does not contain source and target columns")

nodes <- nodes %>%
  mutate(node_id = clean_id(node_id)) %>%
  filter(!is.na(node_id), node_id != "") %>%
  distinct(node_id, .keep_all = TRUE)

edges <- edges %>%
  mutate(
    source = clean_id(source),
    target = clean_id(target),
    interaction = if ("interaction" %in% names(.)) clean_id(interaction) else "interacts_with"
  ) %>%
  filter(!is.na(source), source != "", !is.na(target), target != "")

edge_ids <- unique(c(edges$source, edges$target))
node_ids <- unique(nodes$node_id)
missing_in_nodes <- setdiff(edge_ids, node_ids)
unused_nodes <- setdiff(node_ids, edge_ids)

out_dir <- ensure_out_dir(file.path(.suz_out_root, "Cytoscape_import_safe"))  ## [ported] was a CWD-relative path
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv(nodes, file.path(out_dir, "ALL_MFE20_nodes_CYTOSCAPE_SAFE.tsv"), na = "")
write_tsv(edges, file.path(out_dir, "ALL_MFE20_edges_CYTOSCAPE_SAFE.tsv"), na = "")
write_tsv(edges %>% select(source, target, interaction),
          file.path(out_dir, "ALL_MFE20_edges_MINIMAL.tsv"), na = "")

if (length(missing_in_nodes) > 0) {
  write_tsv(tibble(node_id = missing_in_nodes),
            file.path(out_dir, "missing_edge_endpoints_from_nodes.tsv"))
}
if (length(unused_nodes) > 0) {
  write_tsv(tibble(node_id = unused_nodes),
            file.path(out_dir, "node_ids_not_used_by_edges.tsv"))
}

report <- c(
  "Cytoscape identifier validation report",
  "======================================",
  paste0("Nodes input: ", nodes_file),
  paste0("Edges input: ", edges_file),
  "",
  paste0("Unique node_id values in node file: ", length(node_ids)),
  paste0("Unique endpoint IDs in edge file: ", length(edge_ids)),
  paste0("Edge endpoint IDs absent from node file: ", length(missing_in_nodes)),
  paste0("Node IDs not used by any edge: ", length(unused_nodes)),
  "",
  "First node IDs:",
  paste(head(node_ids, 10), collapse = "\n"),
  "",
  "First edge endpoint IDs:",
  paste(head(edge_ids, 10), collapse = "\n"),
  "",
  if (length(missing_in_nodes) == 0) {
    "RESULT: Every source/target identifier has an exact matching node_id after cleaning."
  } else {
    "RESULT: Some source/target identifiers do not match node_id. Inspect missing_edge_endpoints_from_nodes.tsv."
  },
  "",
  "Recommended Cytoscape import:",
  "1. Import ALL_MFE20_edges_MINIMAL.tsv as a network.",
  "2. Set source = Source Node, target = Target Node, interaction = Interaction Type.",
  "3. Import ALL_MFE20_nodes_CYTOSCAPE_SAFE.tsv as Node Table Columns.",
  "4. Mark node_id as the imported table key and match it to the network column name.",
  "5. If importing into a network collection, try shared name only after confirming name values in the Node Table."
)

writeLines(report, file.path(out_dir, "Cytoscape_ID_match_report.txt"))
message("Done. Clean Cytoscape files written to: ", normalizePath(out_dir))

