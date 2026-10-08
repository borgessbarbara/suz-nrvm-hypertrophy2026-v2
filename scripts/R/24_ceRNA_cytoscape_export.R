#!/usr/bin/env Rscript
# =============================================================================
# R24 — Figure 6B-F and S5: Cytoscape node/edge exports, MFE <= -20, curated KEGG modules
#
# PORTED, NOT REWRITTEN.
# Source: archive/.../32f_export_ceRNA_cytoscape_TF11_curatedKEGG_MFE20.R (1163 lines)
# Produced by scripts/dev/port_from_archive.sh, which applies exactly five
# classes of change and passes everything else through unmodified:
#   1. this preamble (optparse CLI, lib/ sourcing, base_dir/out_dir/dpi)
#   2. hard-coded C:/Users/sebau/... base_dir and out_dir commented out
#   3. base = "C:/Users/sebau/..." default arguments -> base = base_dir
#   4. browser-download duplicates "name(1).csv" removed from candidate lists
#   5. the /mnt/data container search path removed
#
# Diff this file against the archived original to audit the port:
#   diff archive/.../32f_export_ceRNA_cytoscape_TF11_curatedKEGG_MFE20.R scripts/R/22_ceRNA_cytoscape_export.R
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
    make_option("--out-dir", type = "character", default = "results/R24_ceRNA_cytoscape_export",
                help = "output directory [default %default]"),
    make_option("--dpi", type = "integer", default = 600,
                help = "raster figure resolution [default %default]")),
  description = "\nR22: Figure 6B-F and S5: Cytoscape node/edge exports, MFE <= -20, curated KEGG modules\nPorted from 32f_export_ceRNA_cytoscape_TF11_curatedKEGG_MFE20.R."))

base_dir <- normalizePath(.opt$`in-dir`, winslash = "/", mustWork = FALSE)
out_dir  <- ensure_out_dir(.opt$`out-dir`)
.suz_out_root <- out_dir   ## [ported] preserved: the body below rebinds out_dir
png_dpi  <- .opt$dpi
dpi_png  <- .opt$dpi
message("[R24] in-dir : ", base_dir)
message("[R24] out-dir: ", out_dir)

# ----------------------------- ported body ----------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================
# Cytoscape exports for high-confidence ceRNA networks
# FIXED: exact KEGG modules are reconstructed from the long-form mapping;
# no parallel separate_rows() expansion of unequal collapsed fields.
# organized by temporal phase, sampling state, ncRNA type,
# expression axis, miRNA hub, and exact curated KEGG pathways
# already used in the manuscript. No biological process categories
# are inferred from gene names or introduced by this script.
#
# Strict filters:
#   - IntaRNA MFE <= -20 kcal/mol
#   - coherent directions only:
#       ncRNA DOWN / miRNA UP / mRNA DOWN
#       ncRNA UP / miRNA DOWN / mRNA UP
#   - miRNA-mRNA support already present in the integrated triplets
#
# Outputs are Cytoscape-ready CSV node/edge tables.
# ============================================================

out_dir <- ensure_out_dir(file.path(.suz_out_root, "CytoKEGG_TF11_MFE20"))  ## [ported] was a CWD-relative path
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) {
  stop("Could not create output directory: ", out_dir,
       ". Move the analysis to a shorter writable path, for example C:/ceRNA, and rerun.")
}

mfe_cutoff <- -20
min_triplets_for_hub_export <- 3
export_curated_kegg_modules <- TRUE
export_miRNA_hubs <- TRUE

# KEGG export policy:
#   "manuscript_selected" = exact 11 KEGG terms displayed in the TF–KEGG figure
#   "all_curated"         = every pathway present in the three curated KEGG files
#   "custom"              = use custom_selected_pathways below
kegg_export_mode <- "manuscript_selected"

# Exact KEGG descriptions displayed in the manuscript TF–KEGG figure.
# This is an allow-list only. A pathway is exported only when a strict ceRNA
# triplet maps to an exact miRNA–mRNA edge and terminal gene associated with
# that curated pathway. Pathways without ceRNA-layer overlap remain in the
# representation audit as mapped_to_ceRNA = FALSE.
manuscript_selected_pathways <- c(
  # Exact curated KEGG pathways displayed in the manuscript TF–KEGG figure
  "Cytokine-cytokine receptor interaction",
  "Focal adhesion",
  "TGF-beta signaling pathway",
  "Hippo signaling pathway",
  "Leukocyte transendothelial migration",
  "Protein processing in endoplasmic reticulum",
  "Tight junction",
  "Cell cycle",
  "Cytoskeleton in muscle cells",
  "Ether lipid metabolism",
  "PPAR signaling pathway"
)

# Used only when kegg_export_mode == "custom".
custom_selected_pathways <- character()

# Suggested visual range for log2FC continuous mapping
log2fc_min <- -2
log2fc_mid <- 0
log2fc_max <- 2

# -----------------------------
# Helpers
# -----------------------------
find_input <- function(exact_names, regex_pattern = NULL, required = TRUE) {
  for (nm in exact_names) {
    if (!is.na(nm) && file.exists(nm)) return(nm)
  }
  if (!is.null(regex_pattern)) {
    hits <- list.files(".", recursive = TRUE, full.names = TRUE,
                       pattern = regex_pattern, ignore.case = TRUE)
    if (length(hits) > 0) {
      info <- file.info(hits)
      return(hits[order(info$mtime, decreasing = TRUE)][1])
    }
  }
  if (required) stop("Required input not found: ", paste(exact_names, collapse = ", "))
  NA_character_
}

as_bool <- function(x) toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
numify <- function(x) suppressWarnings(as.numeric(as.character(x)))
normalize_mir <- function(x) tolower(trimws(as.character(x)))
normalize_gene <- function(x) toupper(trimws(as.character(x)))

col_chr <- function(df, col, default = NA_character_) {
  if (col %in% names(df)) as.character(df[[col]]) else rep(default, nrow(df))
}
col_num <- function(df, col) {
  if (col %in% names(df)) numify(df[[col]]) else rep(NA_real_, nrow(df))
}
collapse_unique <- function(x, sep = ";") {
  y <- unique(trimws(as.character(x)))
  y <- y[!is.na(y) & y != ""]
  if (length(y) == 0) return(NA_character_)
  paste(y, collapse = sep)
}
first_non_na <- function(x) {
  y <- x[!is.na(x) & x != ""]
  if (length(y) == 0) return(NA)
  y[1]
}
most_extreme_by_abs <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) return(NA_real_)
  y[which.max(abs(y))][1]
}
safe_min <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) return(NA_real_)
  min(y)
}

sanitize_name <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9._-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

shorten_stem <- function(x, max_chars = 80) {
  x <- sanitize_name(x)
  ifelse(nchar(x) > max_chars, substr(x, 1, max_chars), x)
}

safe_write_csv <- function(df, path) {
  parent <- dirname(path)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(parent)) {
    stop("Could not create output folder: ", parent)
  }
  # Windows can fail on long paths. Use a short fallback filename when needed.
  if (.Platform$OS.type == "windows" && nchar(normalizePath(parent, winslash = "/", mustWork = TRUE)) + nchar(basename(path)) > 235) {
    ext <- tools::file_ext(path)
    stem <- tools::file_path_sans_ext(basename(path))
    path <- file.path(parent, paste0(substr(stem, 1, 45), "_short.", ext))
  }
  readr::write_csv(df, path, na = "")
  if (!file.exists(path)) stop("The file was not created: ", path)
  path
}

axis_from_directions <- function(nc_dir, mir_dir, mrna_dir) {
  case_when(
    nc_dir == "DOWN" & mir_dir == "UP"   & mrna_dir == "DOWN" ~ "DOWN–UP–DOWN",
    nc_dir == "UP"   & mir_dir == "DOWN" & mrna_dir == "UP"   ~ "UP–DOWN–UP",
    TRUE ~ NA_character_
  )
}

node_id_lnc  <- function(x) paste0("lnc:", x)
node_id_circ <- function(x) paste0("circ:", x)
node_id_mir  <- function(x) paste0("mir:", normalize_mir(x))
node_id_mrna <- function(x) paste0("mrna:", normalize_gene(x))

shape_map <- c(
  "lncRNA" = "ROUND_RECTANGLE",
  "circRNA" = "ELLIPSE",
  "miRNA" = "DIAMOND",
  "mRNA" = "RECTANGLE"
)

# -----------------------------
# Curated KEGG policy
# -----------------------------
# Biological modules are exported only from exact KEGG pathway
# descriptions already present in the manuscript's curated resources.
# No keyword-based or gene-name-based process assignment is performed.

kegg_group_labels <- c(
  "Only6h" = "6 h state-only",
  "Shared" = "Shared",
  "Only24h" = "24 h state-only"
)

# -----------------------------
# Locate inputs
# -----------------------------
lnc_triplet_file <- find_input(
  c("COHERENT_lncRNA_miRNA_mRNA_triplets.csv",
    
    "Fig8_priority_ceRNA_triplets.txt"),
  regex_pattern = "^(COHERENT_lncRNA_miRNA_mRNA_triplets|Fig8_priority_ceRNA_triplets).*(csv|txt)$"
)

circ_triplet_file <- find_input(
  c("12_STRICT_DE_coherent_circRNA_triplets.csv"),
  regex_pattern = "^12_STRICT_DE_coherent_circRNA_triplets.*\\.csv$",
  required = FALSE
)

lnc6_file <- find_input(c("NE_vs_Ctrl_6h_lncRNA_DE.txt"),
                        "^NE_vs_Ctrl_6h_lncRNA_DE.*\\.txt$")
lnc24_file <- find_input(c("NE_vs_Ctrl_24h_lncRNA_DE.txt"),
                         "^NE_vs_Ctrl_24h_lncRNA_DE.*\\.txt$")
mir6_file <- find_input(c("DE_miRNA_All_with_status_NE6.txt"),
                        "^DE_miRNA_All_with_status_NE6.*\\.txt$")
mir24_file <- find_input(c("DE_miRNA_All_with_status_NE24.txt"),
                         "^DE_miRNA_All_with_status_NE24.*\\.txt$")
mir_temporal_file <- find_input(c("miRNA_temporal_summary_table.csv"),
                               "^miRNA_temporal_summary_table.*\\.csv$")

mrna6_file <- find_input(c("NE_vs_Ctrl_6h_mRNA_DE.txt"),
                         "^NE_vs_Ctrl_6h_mRNA_DE.*\\.txt$")
mrna24_file <- find_input(c("NE_vs_Ctrl_24h_mRNA_DE.txt"),
                          "^NE_vs_Ctrl_24h_mRNA_DE.*\\.txt$")

circ6_file <- find_input(c("06_circRNA_sequences_de6_DE_complete_table.txt"),
                         "^06_circRNA_sequences_de6_DE_complete_table.*\\.txt$",
                         required = FALSE)
circ24_file <- find_input(c("07_circRNA_sequences_de24_DE_complete_table.txt"),
                          "^07_circRNA_sequences_de24_DE_complete_table.*\\.txt$",
                          required = FALSE)

kegg_files <- list.files(".", recursive = TRUE, full.names = TRUE,
                         pattern = "^SUPP_edges_.*_KEGG_phase\\.csv$",
                         ignore.case = TRUE)

curated_only6h_file <- find_input(
  c("KEGG_6h_vs_24h_Only6h_long_format_SYMBOLS.txt"),
  "^KEGG_6h_vs_24h_Only6h_long_format_SYMBOLS.*\\.txt$"
)
curated_shared_file <- find_input(
  c("KEGG_6h_vs_24h_Shared_long_format_SYMBOLS.txt"),
  "^KEGG_6h_vs_24h_Shared_long_format_SYMBOLS.*\\.txt$"
)
curated_only24h_file <- find_input(
  c("KEGG_6h_vs_24h_Only24h_long_format_SYMBOLS.txt"),
  "^KEGG_6h_vs_24h_Only24h_long_format_SYMBOLS.*\\.txt$"
)

# -----------------------------
# Expression and phase maps
# -----------------------------
read_lnc_de <- function(path, tp) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  tibble(
    timepoint = tp,
    gene_id = as.character(x$gene_id),
    ncRNA_log2FC = numify(x$log2FoldChange),
    ncRNA_padj = numify(x$padj),
    ncRNA_pvalue = col_num(x, "pvalue"),
    ncRNA_baseMean = col_num(x, "baseMean")
  ) %>% distinct(timepoint, gene_id, .keep_all = TRUE)
}

read_mir_de <- function(path, tp) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  mir_col <- if ("miRNA_id" %in% names(x)) "miRNA_id" else if ("miRNA" %in% names(x)) "miRNA" else names(x)[1]
  tibble(
    timepoint = tp,
    miRNA_key = normalize_mir(x[[mir_col]]),
    miRNA_label_DE = as.character(x[[mir_col]]),
    miRNA_log2FC = numify(x$log2FoldChange),
    miRNA_padj = numify(x$padj),
    miRNA_pvalue = col_num(x, "pvalue"),
    miRNA_baseMean = col_num(x, "baseMean")
  ) %>% distinct(timepoint, miRNA_key, .keep_all = TRUE)
}

read_mrna_de <- function(path, tp) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  symbol_col <- if ("gene_name" %in% names(x)) "gene_name" else if ("gene_symbol" %in% names(x)) "gene_symbol" else NA_character_
  tibble(
    timepoint = tp,
    gene_id = as.character(x$gene_id),
    gene_symbol = if (!is.na(symbol_col)) as.character(x[[symbol_col]]) else as.character(x$gene_id),
    gene_key = normalize_gene(if (!is.na(symbol_col)) x[[symbol_col]] else x$gene_id),
    mRNA_log2FC_map = numify(x$log2FoldChange),
    mRNA_padj_map = numify(x$padj),
    mRNA_direction_map = toupper(as.character(x$direction))
  ) %>% distinct(timepoint, gene_key, .keep_all = TRUE)
}

lnc_de <- bind_rows(read_lnc_de(lnc6_file, "6h"), read_lnc_de(lnc24_file, "24h"))
mir_de <- bind_rows(read_mir_de(mir6_file, "6h"), read_mir_de(mir24_file, "24h"))
mrna_de <- bind_rows(read_mrna_de(mrna6_file, "6h"), read_mrna_de(mrna24_file, "24h"))

# mRNA temporal class map
mrna_phase_map <- full_join(
  mrna_de %>% filter(timepoint == "6h") %>% select(gene_key, dir_6h = mRNA_direction_map),
  mrna_de %>% filter(timepoint == "24h") %>% select(gene_key, dir_24h = mRNA_direction_map),
  by = "gene_key"
) %>%
  mutate(
    mRNA_phase_map = case_when(
      !is.na(dir_6h) & is.na(dir_24h) ~ "Early",
      is.na(dir_6h) & !is.na(dir_24h) ~ "Late",
      !is.na(dir_6h) & !is.na(dir_24h) ~ "Sustained",
      TRUE ~ "Unassigned"
    ),
    mRNA_direction_switch = !is.na(dir_6h) & !is.na(dir_24h) & dir_6h != dir_24h
  )

# miRNA temporal class map
mir_temporal_raw <- read.csv(mir_temporal_file, stringsAsFactors = FALSE, check.names = FALSE)
mir_phase_map <- mir_temporal_raw %>%
  transmute(
    miRNA_key = normalize_mir(miRNA_id),
    miRNA_phase = case_when(
      temporal_class == "Early_6h_only" ~ "Early",
      temporal_class == "Sustained_both" ~ "Sustained",
      temporal_class == "Late_24h_only" ~ "Late",
      TRUE ~ "Non_DE_or_unassigned"
    ),
    miRNA_temporal_pattern = as.character(temporal_pattern)
  ) %>% distinct(miRNA_key, .keep_all = TRUE)

# circRNA temporal class map
read_circ_ids <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- if ("circ_id" %in% names(x)) "circ_id" else if ("circRNA_ID" %in% names(x)) "circRNA_ID" else names(x)[1]
  unique(as.character(x[[id_col]]))
}

circ_phase_map <- tibble()
if (!is.na(circ6_file) && !is.na(circ24_file)) {
  c6 <- read_circ_ids(circ6_file)
  c24 <- read_circ_ids(circ24_file)
  all_c <- union(c6, c24)
  circ_phase_map <- tibble(
    circ_locus = all_c,
    ncRNA_phase = case_when(
      all_c %in% c6 & !(all_c %in% c24) ~ "Early",
      all_c %in% c6 & all_c %in% c24 ~ "Sustained",
      !(all_c %in% c6) & all_c %in% c24 ~ "Late",
      TRUE ~ "Unassigned"
    )
  )
}

# Exact curated KEGG reference and miRNA-mRNA pathway mapping
read_curated_kegg_reference <- function(path, group_raw) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

  if (group_raw == "Shared") {
    tibble(
      pathway_id = as.character(x$ID),
      pathway_description = as.character(x$Description),
      kegg_program_raw = group_raw,
      kegg_program = unname(kegg_group_labels[group_raw]),
      NES_6h = numify(x$NES_6h),
      NES_24h = numify(x$NES_24h),
      pathway_padj_6h_reference = numify(x$p.adjust_6h),
      pathway_padj_24h_reference = numify(x$p.adjust_24h),
      gene_key_reference = normalize_gene(x$gene_symbol),
      curated_reference_file = basename(path)
    )
  } else if (group_raw == "Only6h") {
    tibble(
      pathway_id = as.character(x$ID),
      pathway_description = as.character(x$Description),
      kegg_program_raw = group_raw,
      kegg_program = unname(kegg_group_labels[group_raw]),
      NES_6h = numify(x$NES),
      NES_24h = NA_real_,
      pathway_padj_6h_reference = numify(x$p.adjust),
      pathway_padj_24h_reference = NA_real_,
      gene_key_reference = normalize_gene(x$gene_symbol),
      curated_reference_file = basename(path)
    )
  } else {
    tibble(
      pathway_id = as.character(x$ID),
      pathway_description = as.character(x$Description),
      kegg_program_raw = group_raw,
      kegg_program = unname(kegg_group_labels[group_raw]),
      NES_6h = NA_real_,
      NES_24h = numify(x$NES),
      pathway_padj_6h_reference = NA_real_,
      pathway_padj_24h_reference = numify(x$p.adjust),
      gene_key_reference = normalize_gene(x$gene_symbol),
      curated_reference_file = basename(path)
    )
  }
}

curated_kegg_reference_gene <- bind_rows(
  read_curated_kegg_reference(curated_only6h_file, "Only6h"),
  read_curated_kegg_reference(curated_shared_file, "Shared"),
  read_curated_kegg_reference(curated_only24h_file, "Only24h")
) %>%
  distinct(kegg_program_raw, pathway_id, pathway_description, gene_key_reference, .keep_all = TRUE)


available_curated_pathways <- sort(unique(curated_kegg_reference_gene$pathway_description))
selected_pathway_allowlist <- switch(
  kegg_export_mode,
  "all_curated" = available_curated_pathways,
  "manuscript_selected" = manuscript_selected_pathways,
  "custom" = custom_selected_pathways,
  stop("Unknown kegg_export_mode: ", kegg_export_mode)
)
selected_pathway_allowlist <- unique(trimws(selected_pathway_allowlist))
missing_selected_terms <- setdiff(selected_pathway_allowlist, available_curated_pathways)
if (length(missing_selected_terms) > 0) {
  warning("The following selected pathway descriptions were not found exactly in the curated KEGG files: ",
          paste(missing_selected_terms, collapse = "; "))
}
selected_pathway_allowlist <- intersect(selected_pathway_allowlist, available_curated_pathways)

curated_kegg_reference_gene_selected <- curated_kegg_reference_gene %>%
  filter(pathway_description %in% selected_pathway_allowlist)

curated_kegg_reference_pathway <- curated_kegg_reference_gene %>%
  group_by(kegg_program_raw, kegg_program, pathway_id, pathway_description) %>%
  summarise(
    NES_6h = first_non_na(NES_6h),
    NES_24h = first_non_na(NES_24h),
    pathway_padj_6h_reference = safe_min(pathway_padj_6h_reference),
    pathway_padj_24h_reference = safe_min(pathway_padj_24h_reference),
    curated_reference_files = collapse_unique(curated_reference_file),
    n_reference_genes = n_distinct(gene_key_reference),
    .groups = "drop"
  )

kegg_map <- tibble()
if (length(kegg_files) > 0) {
  kegg_map_long <- map_dfr(kegg_files, function(f) {
    x <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    tibble(
      timepoint = as.character(x$time),
      miRNA_key = normalize_mir(x$miRNA),
      gene_key = normalize_gene(x$gene_key),
      kegg_program_raw = as.character(x$KEGG_group),
      pathway_id = as.character(x$ID),
      pathway_description = as.character(x$Description),
      edge_pathway_padj = numify(x$padj),
      kegg_edge_source_file = basename(f)
    )
  }) %>%
    inner_join(
      curated_kegg_reference_gene_selected,
      by = c(
        "kegg_program_raw",
        "pathway_id",
        "pathway_description",
        "gene_key" = "gene_key_reference"
      )
    )

  kegg_map <- kegg_map_long %>%
    group_by(timepoint, miRNA_key, gene_key) %>%
    summarise(
      curated_pathway_ids = collapse_unique(pathway_id),
      curated_pathway_descriptions = collapse_unique(pathway_description),
      curated_kegg_programs_raw = collapse_unique(kegg_program_raw),
      curated_kegg_programs = collapse_unique(kegg_program),
      curated_pathway_min_padj = safe_min(edge_pathway_padj),
      curated_NES_6h = collapse_unique(NES_6h),
      curated_NES_24h = collapse_unique(NES_24h),
      curated_kegg_edge_source_files = collapse_unique(kegg_edge_source_file),
      curated_kegg_reference_files = collapse_unique(curated_reference_file),
      .groups = "drop"
    )
}

# -----------------------------
# Standardize lncRNA triplets
# -----------------------------
read_lnc_triplets <- function(path) {
  if (grepl("Fig8_priority_ceRNA_triplets", basename(path), ignore.case = TRUE)) {
    x <- read.delim(path, skip = 2, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  tibble(
    source_scope = ifelse(grepl("priority", basename(path), ignore.case = TRUE), "priority", "complete"),
    ncRNA_type = "lncRNA",
    timepoint = as.character(x$timepoint),
    ncRNA_id = as.character(x$gene_id),
    ncRNA_label = if ("ncRNA_label" %in% names(x)) as.character(x$ncRNA_label) else as.character(x$gene_id),
    transcript_id = col_chr(x, "transcript_id"),
    ncRNA_phase = stringr::str_to_title(trimws(col_chr(x, "phase"))),
    ncRNA_direction = as.character(x$ncRNA_direction),
    miRNA = as.character(x$miRNA),
    miRNA_key = if ("miRNA_key" %in% names(x)) normalize_mir(x$miRNA_key) else normalize_mir(x$miRNA),
    miRNA_direction = as.character(x$miRNA_direction),
    MFE = if ("E" %in% names(x)) numify(x$E) else col_num(x, "best_MFE"),
    mRNA = if ("gene_symbol" %in% names(x)) as.character(x$gene_symbol) else col_chr(x, "mRNA_gene_final"),
    gene_key = normalize_gene(if ("gene_symbol" %in% names(x)) x$gene_symbol else col_chr(x, "mRNA_gene_final")),
    mRNA_phase_input = stringr::str_to_title(trimws(col_chr(x, "phase_mRNA"))),
    mRNA_direction = if ("mRNA_direction" %in% names(x)) as.character(x$mRNA_direction) else col_chr(x, "mRNA_DE_direction"),
    mRNA_log2FC = col_num(x, "mRNA_log2FC"),
    mRNA_padj = col_num(x, "mRNA_padj"),
    coherent_triplet = if ("coherent_triplet" %in% names(x)) as_bool(x$coherent_triplet) else rep(TRUE, nrow(x)),
    target_region = col_chr(x, "region"),
    target_score = col_chr(x, "target_score"),
    support_files = col_chr(x, "support_file")
  ) %>%
    mutate(axis = axis_from_directions(ncRNA_direction, miRNA_direction, mRNA_direction)) %>%
    filter(
      coherent_triplet,
      !is.na(MFE), MFE <= mfe_cutoff,
      ncRNA_phase %in% c("Early", "Sustained", "Late"),
      !is.na(axis)
    ) %>%
    distinct(timepoint, ncRNA_id, transcript_id, miRNA_key, gene_key, .keep_all = TRUE) %>%
    left_join(lnc_de, by = c("timepoint", "ncRNA_id" = "gene_id")) %>%
    left_join(mir_de, by = c("timepoint", "miRNA_key")) %>%
    left_join(mir_phase_map, by = "miRNA_key") %>%
    left_join(mrna_phase_map, by = "gene_key") %>%
    mutate(
      mRNA_phase = case_when(
        mRNA_phase_input %in% c("Early", "Sustained", "Late") ~ mRNA_phase_input,
        TRUE ~ mRNA_phase_map
      )
    )
}

lnc_trip <- read_lnc_triplets(lnc_triplet_file)

# -----------------------------
# Standardize circRNA triplets
# -----------------------------
read_circ_triplets <- function(path) {
  if (is.na(path) || !file.exists(path)) return(tibble())
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  tibble(
    source_scope = "strict_circRNA",
    ncRNA_type = "circRNA",
    timepoint = as.character(x$timepoint),
    ncRNA_id = as.character(x$circ_locus),
    ncRNA_label = as.character(x$circ_locus),
    transcript_id = col_chr(x, "candidate_variant_ids"),
    ncRNA_direction = as.character(x$circ_direction),
    ncRNA_log2FC = numify(x$circ_log2FC),
    ncRNA_padj = numify(x$circ_padj),
    ncRNA_pvalue = numify(x$circ_pvalue),
    ncRNA_baseMean = NA_real_,
    miRNA = as.character(x$miRNA),
    miRNA_key = normalize_mir(x$miRNA_key),
    miRNA_direction = as.character(x$miRNA_DE_direction),
    miRNA_log2FC = numify(x$miRNA_log2FC),
    miRNA_padj = numify(x$miRNA_padj),
    miRNA_pvalue = NA_real_,
    miRNA_baseMean = NA_real_,
    MFE = numify(x$best_MFE),
    mRNA = as.character(x$mRNA_gene_final),
    gene_key = normalize_gene(x$mRNA_gene_final),
    mRNA_direction = as.character(x$mRNA_DE_direction),
    mRNA_log2FC = numify(x$mRNA_log2FC),
    mRNA_padj = numify(x$mRNA_padj),
    coherent_triplet = as_bool(x$strict_DE_coherent_circRNA_triplet),
    target_region = as.character(x$target_regions),
    target_score = as.character(x$target_scores),
    support_files = as.character(x$target_source_files)
  ) %>%
    left_join(circ_phase_map, by = c("ncRNA_id" = "circ_locus")) %>%
    mutate(
      axis = axis_from_directions(ncRNA_direction, miRNA_direction, mRNA_direction)
    ) %>%
    filter(
      coherent_triplet,
      !is.na(MFE), MFE <= mfe_cutoff,
      ncRNA_phase %in% c("Early", "Sustained", "Late"),
      !is.na(axis)
    ) %>%
    distinct(timepoint, ncRNA_id, miRNA_key, gene_key, .keep_all = TRUE) %>%
    left_join(mir_phase_map, by = "miRNA_key") %>%
    left_join(mrna_phase_map, by = "gene_key") %>%
    mutate(mRNA_phase = mRNA_phase_map)
}

circ_trip <- read_circ_triplets(circ_triplet_file)

# -----------------------------
# Combine triplets and attach exact curated KEGG annotations
# -----------------------------
all_trip <- bind_rows(lnc_trip, circ_trip) %>%
  mutate(
    triplet_id = paste(timepoint, ncRNA_type, ncRNA_id, miRNA_key, gene_key, sep = "__"),
    phase_relation = paste0("ncRNA_", ncRNA_phase, "__miRNA_", miRNA_phase, "__mRNA_", mRNA_phase),
    node_ncrna_id = ifelse(ncRNA_type == "lncRNA", node_id_lnc(ncRNA_id), node_id_circ(ncRNA_id)),
    node_mir_id = node_id_mir(miRNA),
    node_mrna_id = node_id_mrna(mRNA)
  )

if (nrow(kegg_map) > 0) {
  all_trip <- all_trip %>%
    left_join(kegg_map, by = c("timepoint", "miRNA_key", "gene_key"))
} else {
  all_trip <- all_trip %>%
    mutate(
      curated_pathway_ids = NA_character_,
      curated_pathway_descriptions = NA_character_,
      curated_kegg_programs_raw = NA_character_,
      curated_kegg_programs = NA_character_,
      curated_pathway_min_padj = NA_real_,
      curated_NES_6h = NA_character_,
      curated_NES_24h = NA_character_,
      curated_kegg_edge_source_files = NA_character_,
      curated_kegg_reference_files = NA_character_
    )
}

all_trip <- all_trip %>%
  mutate(
    curated_KEGG_status = ifelse(
      !is.na(curated_pathway_descriptions) & curated_pathway_descriptions != "",
      "Mapped_to_curated_KEGG",
      "Not_mapped_to_curated_KEGG"
    ),
    miRNA_module = paste0("module_", normalize_mir(miRNA))
  )

safe_write_csv(all_trip,
               file.path(out_dir, "01_ALL_filtered_triplets_with_phase_curatedKEGG_annotations.csv"))
safe_write_csv(filter(all_trip, curated_KEGG_status == "Not_mapped_to_curated_KEGG"),
               file.path(out_dir, "01B_triplets_not_mapped_to_selected_curated_KEGG.csv"))
safe_write_csv(curated_kegg_reference_pathway,
               file.path(out_dir, "01C_all_curated_KEGG_pathway_reference.csv"))
safe_write_csv(curated_kegg_reference_pathway %>% filter(pathway_description %in% selected_pathway_allowlist),
               file.path(out_dir, "01C2_selected_curated_KEGG_pathway_reference.csv"))
safe_write_csv(tibble(selected_pathway_description = selected_pathway_allowlist),
               file.path(out_dir, "01C3_selected_KEGG_allowlist_used.csv"))

# -----------------------------
# Node/edge builders
# -----------------------------
make_nodes <- function(trip) {
  if (nrow(trip) == 0) return(tibble())

  nc_nodes <- trip %>%
    transmute(
      timepoint,
      node_id = node_ncrna_id,
      label = ncRNA_label,
      node_type = ncRNA_type,
      node_shape = ifelse(ncRNA_type == "lncRNA", shape_map[["lncRNA"]], shape_map[["circRNA"]]),
      direction = ncRNA_direction,
      log2FC = ncRNA_log2FC,
      padj = ncRNA_padj,
      pvalue = ncRNA_pvalue,
      temporal_phase = ncRNA_phase,
      axis,
      phase_relation,
      curated_KEGG_status,
      curated_pathway_ids,
      curated_pathways = curated_pathway_descriptions,
      curated_KEGG_programs = curated_kegg_programs,
      curated_NES_6h,
      curated_NES_24h,
      miRNA_modules = miRNA_module,
      triplet_id
    )

  mir_nodes <- trip %>%
    transmute(
      timepoint,
      node_id = node_mir_id,
      label = miRNA,
      node_type = "miRNA",
      node_shape = shape_map[["miRNA"]],
      direction = miRNA_direction,
      log2FC = miRNA_log2FC,
      padj = miRNA_padj,
      pvalue = miRNA_pvalue,
      temporal_phase = miRNA_phase,
      axis,
      phase_relation,
      curated_KEGG_status,
      curated_pathway_ids,
      curated_pathways = curated_pathway_descriptions,
      curated_KEGG_programs = curated_kegg_programs,
      curated_NES_6h,
      curated_NES_24h,
      miRNA_modules = miRNA_module,
      triplet_id
    )

  mrna_nodes <- trip %>%
    transmute(
      timepoint,
      node_id = node_mrna_id,
      label = mRNA,
      node_type = "mRNA",
      node_shape = shape_map[["mRNA"]],
      direction = mRNA_direction,
      log2FC = mRNA_log2FC,
      padj = mRNA_padj,
      pvalue = NA_real_,
      temporal_phase = mRNA_phase,
      axis,
      phase_relation,
      curated_KEGG_status,
      curated_pathway_ids,
      curated_pathways = curated_pathway_descriptions,
      curated_KEGG_programs = curated_kegg_programs,
      curated_NES_6h,
      curated_NES_24h,
      miRNA_modules = miRNA_module,
      triplet_id
    )

  bind_rows(nc_nodes, mir_nodes, mrna_nodes) %>%
    group_by(node_id, label, node_type, node_shape) %>%
    summarise(
      timepoints = collapse_unique(timepoint),
      direction = collapse_unique(direction),
      log2FC = most_extreme_by_abs(log2FC),
      padj = safe_min(padj),
      pvalue = safe_min(pvalue),
      temporal_phases = collapse_unique(temporal_phase),
      axes = collapse_unique(axis),
      phase_relations = collapse_unique(phase_relation),
      curated_KEGG_status = collapse_unique(curated_KEGG_status),
      curated_pathway_ids = collapse_unique(curated_pathway_ids),
      curated_pathways = collapse_unique(curated_pathways),
      curated_KEGG_programs = collapse_unique(curated_KEGG_programs),
      curated_NES_6h = collapse_unique(curated_NES_6h),
      curated_NES_24h = collapse_unique(curated_NES_24h),
      miRNA_modules = collapse_unique(miRNA_modules),
      n_supported_triplets = n_distinct(triplet_id),
      .groups = "drop"
    ) %>%
    mutate(
      preferred_log2FC_range = paste(log2fc_min, log2fc_mid, log2fc_max, sep = ";")
    )
}

make_edges <- function(trip) {
  if (nrow(trip) == 0) return(tibble())

  nc_mir <- trip %>%
    transmute(
      source = node_ncrna_id,
      target = node_mir_id,
      interaction = paste0(ncRNA_type, "-miRNA"),
      edge_type = "ncRNA-miRNA",
      timepoint,
      ncRNA_type,
      ncRNA_phase,
      axis,
      phase_relation,
      MFE,
      target_region = NA_character_,
      target_score = NA_character_,
      curated_KEGG_status,
      curated_pathway_ids,
      curated_pathways = curated_pathway_descriptions,
      curated_KEGG_programs = curated_kegg_programs,
      curated_NES_6h,
      curated_NES_24h,
      support_files,
      triplet_id
    )

  mir_mrna <- trip %>%
    transmute(
      source = node_mir_id,
      target = node_mrna_id,
      interaction = "miRNA-mRNA",
      edge_type = "miRNA-mRNA",
      timepoint,
      ncRNA_type,
      ncRNA_phase,
      axis,
      phase_relation,
      MFE = NA_real_,
      target_region,
      target_score,
      curated_KEGG_status,
      curated_pathway_ids,
      curated_pathways = curated_pathway_descriptions,
      curated_KEGG_programs = curated_kegg_programs,
      curated_NES_6h,
      curated_NES_24h,
      support_files,
      triplet_id
    )

  bind_rows(nc_mir, mir_mrna) %>%
    group_by(source, target, interaction, edge_type) %>%
    summarise(
      timepoints = collapse_unique(timepoint),
      ncRNA_types = collapse_unique(ncRNA_type),
      ncRNA_phases = collapse_unique(ncRNA_phase),
      axes = collapse_unique(axis),
      phase_relations = collapse_unique(phase_relation),
      MFE = safe_min(MFE),
      target_regions = collapse_unique(target_region),
      target_scores = collapse_unique(target_score),
      curated_KEGG_status = collapse_unique(curated_KEGG_status),
      curated_pathway_ids = collapse_unique(curated_pathway_ids),
      curated_pathways = collapse_unique(curated_pathways),
      curated_KEGG_programs = collapse_unique(curated_KEGG_programs),
      curated_NES_6h = collapse_unique(curated_NES_6h),
      curated_NES_24h = collapse_unique(curated_NES_24h),
      support_files = collapse_unique(support_files),
      n_supported_triplets = n_distinct(triplet_id),
      supported_triplets = collapse_unique(triplet_id),
      .groups = "drop"
    )
}

# -----------------------------
# Export function and manifest
# -----------------------------
manifest <- list()

export_scope <- function(trip, scope_name, folder_name, file_stem = NULL) {
  if (nrow(trip) == 0) return(invisible(NULL))
  folder <- file.path(out_dir, folder_name)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(folder)) stop("Could not create module folder: ", folder)

  nodes <- make_nodes(trip)
  edges <- make_edges(trip)

  safe_scope <- if (is.null(file_stem)) shorten_stem(scope_name) else shorten_stem(file_stem)
  node_file <- file.path(folder, paste0(safe_scope, "_nodes.csv"))
  edge_file <- file.path(folder, paste0(safe_scope, "_edges.csv"))
  trip_file <- file.path(folder, paste0(safe_scope, "_triplets.csv"))

  node_file <- safe_write_csv(nodes, node_file)
  edge_file <- safe_write_csv(edges, edge_file)
  trip_file <- safe_write_csv(trip, trip_file)

  manifest[[length(manifest) + 1]] <<- tibble(
    scope = scope_name,
    folder = folder_name,
    n_triplet_rows = nrow(trip),
    n_unique_triplets = n_distinct(trip$triplet_id),
    n_nodes = nrow(nodes),
    n_edges = nrow(edges),
    node_file = node_file,
    edge_file = edge_file,
    triplet_file = trip_file
  )
}

# All network
export_scope(all_trip, "ALL_MFE20", "00_all")

# By timepoint
for (tp in c("6h", "24h")) {
  export_scope(filter(all_trip, timepoint == tp), paste0("timepoint_", tp), "01_by_timepoint")
}

# By ncRNA phase
for (ph in c("Early", "Sustained", "Late")) {
  export_scope(filter(all_trip, ncRNA_phase == ph), paste0("phase_", ph), "02_by_ncRNA_phase")
}

# By timepoint and phase
phase_tp <- distinct(all_trip, timepoint, ncRNA_phase)
for (i in seq_len(nrow(phase_tp))) {
  tp <- phase_tp$timepoint[i]
  ph <- phase_tp$ncRNA_phase[i]
  export_scope(
    filter(all_trip, timepoint == tp, ncRNA_phase == ph),
    paste0("timepoint_", tp, "__phase_", ph),
    "03_by_timepoint_and_phase"
  )
}

# By ncRNA type
for (typ in c("lncRNA", "circRNA")) {
  export_scope(filter(all_trip, ncRNA_type == typ), paste0("ncRNA_type_", typ), "04_by_ncRNA_type")
}

# By expression axis
for (ax in c("DOWN–UP–DOWN", "UP–DOWN–UP")) {
  export_scope(filter(all_trip, axis == ax), paste0("axis_", ax), "05_by_expression_axis")
}

# By phase and expression axis
phase_axis <- distinct(all_trip, ncRNA_phase, axis)
for (i in seq_len(nrow(phase_axis))) {
  ph <- phase_axis$ncRNA_phase[i]
  ax <- phase_axis$axis[i]
  export_scope(
    filter(all_trip, ncRNA_phase == ph, axis == ax),
    paste0("phase_", ph, "__axis_", ax),
    "06_by_phase_and_axis"
  )
}

# By exact selected curated KEGG pathway
if (export_curated_kegg_modules && exists("kegg_map_long") && nrow(kegg_map_long) > 0) {

  # Collapse the exact KEGG map to one row per timepoint + miRNA + mRNA + pathway.
  # Duplicate target regions/support files are retained as collapsed evidence,
  # rather than creating repeated many-to-many rows.
  kegg_exact_unique <- kegg_map_long %>%
    transmute(
      timepoint = as.character(timepoint),
      miRNA_key = normalize_mir(miRNA_key),
      gene_key = normalize_gene(gene_key),
      pathway_id = as.character(pathway_id),
      pathway_description = as.character(pathway_description),
      kegg_program_raw = as.character(kegg_program_raw),
      kegg_program = as.character(kegg_program),
      edge_pathway_padj = numify(edge_pathway_padj),
      NES_6h = numify(NES_6h),
      NES_24h = numify(NES_24h),
      kegg_edge_source_file = as.character(kegg_edge_source_file),
      curated_reference_file = as.character(curated_reference_file)
    ) %>%
    filter(pathway_description %in% selected_pathway_allowlist) %>%
    group_by(timepoint, miRNA_key, gene_key, pathway_id, pathway_description,
             kegg_program_raw, kegg_program) %>%
    summarise(
      edge_pathway_padj = safe_min(edge_pathway_padj),
      NES_6h = first_non_na(NES_6h),
      NES_24h = first_non_na(NES_24h),
      kegg_edge_source_file = collapse_unique(kegg_edge_source_file),
      curated_reference_file = collapse_unique(curated_reference_file),
      .groups = "drop"
    )

  # One row per unique ceRNA triplet before pathway expansion.
  triplet_core <- all_trip %>%
    distinct(triplet_id, .keep_all = TRUE)

  # A triplet may legitimately map to several exact curated pathways, and a
  # miRNA-mRNA pair may occur in triplets with several ncRNAs. The relationship
  # is therefore intentionally many-to-many at the key level, followed by an
  # explicit one-row-per-triplet-per-pathway deduplication.
  triplet_to_pathway <- triplet_core %>%
    select(triplet_id, timepoint, miRNA_key, gene_key) %>%
    distinct() %>%
    inner_join(
      kegg_exact_unique,
      by = c("timepoint", "miRNA_key", "gene_key"),
      relationship = "many-to-many"
    ) %>%
    distinct(triplet_id, pathway_id, kegg_program, .keep_all = TRUE) %>%
    # Keep only the KEGG-side attributes plus triplet_id.  The biological
    # triplet identifiers (timepoint, miRNA_key, gene_key, mRNA, etc.) are
    # already present in triplet_core.  Removing them here prevents dplyr from
    # creating gene_key.x/gene_key.y during the next join.
    select(
      triplet_id,
      pathway_id,
      pathway_description,
      kegg_program_raw,
      kegg_program,
      edge_pathway_padj,
      NES_6h,
      NES_24h,
      kegg_edge_source_file,
      curated_reference_file
    )

  pathway_expanded <- triplet_core %>%
    select(
      -curated_pathway_ids,
      -curated_pathway_descriptions,
      -curated_kegg_programs_raw,
      -curated_kegg_programs,
      -curated_pathway_min_padj,
      -curated_NES_6h,
      -curated_NES_24h,
      -curated_kegg_edge_source_files,
      -curated_kegg_reference_files,
      -curated_KEGG_status
    ) %>%
    inner_join(triplet_to_pathway, by = "triplet_id", relationship = "one-to-many") %>%
    mutate(
      curated_KEGG_status = "Mapped_to_selected_curated_KEGG",
      curated_pathway_ids = pathway_id,
      curated_pathway_descriptions = pathway_description,
      curated_kegg_programs_raw = kegg_program_raw,
      curated_kegg_programs = kegg_program,
      curated_pathway_min_padj = edge_pathway_padj,
      curated_NES_6h = NES_6h,
      curated_NES_24h = NES_24h,
      curated_kegg_edge_source_files = kegg_edge_source_file,
      curated_kegg_reference_files = curated_reference_file
    ) %>%
    distinct(triplet_id, curated_pathway_ids, curated_kegg_programs, .keep_all = TRUE)

  required_pathway_cols <- c("triplet_id", "gene_key", "mRNA", "curated_pathway_ids",
                             "curated_pathway_descriptions", "curated_kegg_programs")
  missing_pathway_cols <- setdiff(required_pathway_cols, names(pathway_expanded))
  if (length(missing_pathway_cols) > 0) {
    stop(
      "Internal pathway-table error: missing columns after KEGG join: ",
      paste(missing_pathway_cols, collapse = ", "),
      ". Check for duplicated .x/.y column names in the join."
    )
  }

  safe_write_csv(
    pathway_expanded,
    file.path(out_dir, "01D_triplet_to_selected_curated_KEGG.csv")
  )

  # Audit which selected pathways are represented by terminal mRNAs.
  pathway_representation <- tibble(
    pathway_description = selected_pathway_allowlist
  ) %>%
    left_join(
      pathway_expanded %>%
        group_by(curated_pathway_descriptions, curated_pathway_ids, curated_kegg_programs) %>%
        summarise(
          n_triplets = n_distinct(triplet_id),
          n_terminal_mRNAs = n_distinct(gene_key),
          terminal_mRNAs = collapse_unique(mRNA),
          .groups = "drop"
        ),
      by = c("pathway_description" = "curated_pathway_descriptions")
    ) %>%
    mutate(mapped_to_ceRNA = !is.na(n_triplets))

  safe_write_csv(
    pathway_representation,
    file.path(out_dir, "01E_selected_KEGG_pathway_representation_audit.csv")
  )

  # Export one network per exact pathway. Use the KEGG ID as the short file stem
  # to avoid Windows path-length failures; the exact description remains in the
  # manifest and network tables.
  pathway_index <- pathway_expanded %>%
    distinct(curated_pathway_ids, curated_pathway_descriptions) %>%
    arrange(curated_pathway_ids)

  for (i in seq_len(nrow(pathway_index))) {
    pid <- pathway_index$curated_pathway_ids[i]
    pdesc <- pathway_index$curated_pathway_descriptions[i]
    export_scope(
      filter(pathway_expanded, curated_pathway_ids == pid,
             curated_pathway_descriptions == pdesc),
      paste0("curated_KEGG_", pid, "__", pdesc),
      "07_KEGG_pathways",
      file_stem = paste0("KEGG_", pid)
    )
  }

  program_file_stems <- c(
    "6 h state-only" = "KEGG_program_6h",
    "Shared" = "KEGG_program_shared",
    "24 h state-only" = "KEGG_program_24h"
  )
  for (program_name in names(program_file_stems)) {
    export_scope(
      filter(pathway_expanded, curated_kegg_programs == program_name),
      paste0("curated_KEGG_program_", program_name),
      "08_KEGG_programs",
      file_stem = unname(program_file_stems[program_name])
    )
  }
} else if (export_curated_kegg_modules) {
  warning("Selected curated KEGG modules were requested, but no exact long-form KEGG mapping was available.")
}

# By miRNA hub
if (export_miRNA_hubs) {
  hub_counts <- all_trip %>%
    count(miRNA_key, miRNA, name = "n_triplets") %>%
    filter(n_triplets >= min_triplets_for_hub_export)

  for (i in seq_len(nrow(hub_counts))) {
    key <- hub_counts$miRNA_key[i]
    label <- hub_counts$miRNA[i]
    export_scope(
      filter(all_trip, miRNA_key == key),
      paste0("miRNA_hub_", label),
      "09_by_miRNA_hub"
    )
  }
}

manifest_tbl <- bind_rows(manifest)
safe_write_csv(manifest_tbl, file.path(out_dir, "02_EXPORT_manifest.csv"))

# -----------------------------
# Cytoscape visual-style template
# -----------------------------
style_template <- tribble(
  ~visual_property, ~mapping_column, ~mapping_type, ~key_or_range, ~visual_value,
  "NODE_SHAPE", "node_type", "discrete", "lncRNA", shape_map[["lncRNA"]],
  "NODE_SHAPE", "node_type", "discrete", "circRNA", shape_map[["circRNA"]],
  "NODE_SHAPE", "node_type", "discrete", "miRNA", shape_map[["miRNA"]],
  "NODE_SHAPE", "node_type", "discrete", "mRNA", shape_map[["mRNA"]],
  "NODE_FILL_COLOR", "log2FC", "continuous", as.character(log2fc_min), "#2166AC",
  "NODE_FILL_COLOR", "log2FC", "continuous", as.character(log2fc_mid), "#F7F7F7",
  "NODE_FILL_COLOR", "log2FC", "continuous", as.character(log2fc_max), "#B2182B",
  "EDGE_STROKE_UNSELECTED_PAINT", "edge_type", "discrete", "ncRNA-miRNA", "#555555",
  "EDGE_STROKE_UNSELECTED_PAINT", "edge_type", "discrete", "miRNA-mRNA", "#2C7FB8",
  "EDGE_WIDTH", "n_supported_triplets", "continuous", "minimum", "2",
  "EDGE_WIDTH", "n_supported_triplets", "continuous", "maximum", "8",
  "NODE_BORDER_WIDTH", "n_supported_triplets", "continuous", "minimum", "1.5",
  "NODE_BORDER_WIDTH", "n_supported_triplets", "continuous", "maximum", "6"
)
safe_write_csv(style_template, file.path(out_dir, "03_Cytoscape_visual_style_template.csv"))

# Exact curated KEGG pathway reference used for all pathway modules
safe_write_csv(
  curated_kegg_reference_pathway %>% filter(pathway_description %in% selected_pathway_allowlist),
  file.path(out_dir, "04_selected_curated_KEGG_pathways_used_for_modules.csv")
)

# -----------------------------
# README
# -----------------------------
readme <- c(
  "High-confidence ceRNA Cytoscape exports organized by phase and exact curated KEGG pathways",
  "",
  paste0("MFE cutoff: <= ", mfe_cutoff),
  paste0("lncRNA triplet input: ", lnc_triplet_file),
  paste0("circRNA triplet input: ", ifelse(is.na(circ_triplet_file), "not used", circ_triplet_file)),
  "",
  "IMPORTANT CORRECTION",
  "No biological process modules are inferred from gene names or keyword rules.",
  "Every pathway-specific network is named with the exact KEGG Description already present in the curated manuscript resources.",
  "Triplets without an exact curated KEGG mapping remain in temporal and ncRNA-type exports but are not assigned to a pathway module.",
  "",
  "NETWORK ORGANIZATION",
  "00_all: all retained triplets",
  "01_by_timepoint: 6 h and 24 h",
  "02_by_ncRNA_phase: Early, Sustained, Late",
  "03_by_timepoint_and_phase: state-specific phase networks",
  "04_by_ncRNA_type: lncRNA and circRNA",
  "05_by_expression_axis: DOWN-UP-DOWN and UP-DOWN-UP",
  "06_by_phase_and_axis: temporal phase combined with coherence axis",
  "07_KEGG_pathways: one network per exact selected curated KEGG pathway",
  "08_by_curated_KEGG_program: 6 h state-only, Shared, and 24 h state-only pathway programs",
  "09_by_miRNA_hub: miRNA-centered subnetworks",
  "",
  "KEGG SOURCES",
  basename(curated_only6h_file),
  basename(curated_shared_file),
  basename(curated_only24h_file),
  "SUPP_edges_*_KEGG_phase.csv files connect supported miRNA-mRNA edges to those exact pathways.",
  "",
  "NODE VISUAL MAPPINGS",
  "lncRNA = round rectangle",
  "circRNA = ellipse",
  "miRNA = diamond",
  "mRNA = rectangle",
  "continuous log2FC fill: blue -> white -> red",
  "",
  "Use 02_EXPORT_manifest.csv to inspect the size and location of every exported network."
)
readme <- c(readme,
  "",
  "TF–KEGG FIGURE CONSISTENCY",
  "Only the 11 exact KEGG pathways shown in the manuscript TF–KEGG figure are eligible for pathway-module export.",
  "A pathway absent from 07_KEGG_pathways did not overlap the strict ceRNA layer after all filters; it was not deleted from the upstream KEGG analysis.",
  "Salivary secretion, Neuroactive ligand signaling, and other pathways outside this 11-term list are intentionally excluded from this run."
)

writeLines(readme, file.path(out_dir, "README_Cytoscape_TF11_curatedKEGG_exports.txt"))

message("Done. Phase- and TF11-curated-KEGG Cytoscape exports written to: ", out_dir)
