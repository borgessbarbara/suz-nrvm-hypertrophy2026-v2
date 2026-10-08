#!/usr/bin/env Rscript

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

out_dir <- "Cytoscape_import_safe"
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

