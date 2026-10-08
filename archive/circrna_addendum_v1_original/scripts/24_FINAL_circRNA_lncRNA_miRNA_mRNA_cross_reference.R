suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(scales)
})

base_dir <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis"
out_dir <- file.path(base_dir, "CircRNA_ceRNA_cross_reference_FINAL")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

energy_cutoff <- -20
alpha <- 0.05
mirna_lfc_cutoff <- 0.5
mrna_lfc_cutoff <- 1
png_dpi <- 600

all_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE)

find_first <- function(candidates) {
  for (nm in candidates) {
    hit <- all_files[basename(all_files) == nm]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

find_pattern <- function(pattern) {
  all_files[str_detect(basename(all_files), regex(pattern, ignore_case = TRUE))]
}

read_auto <- function(path) {
  if (is.na(path) || !file.exists(path)) return(tibble())
  first_line <- readLines(path, n = 1, warn = FALSE)
  if (length(first_line) == 0) return(tibble())
  if (str_count(first_line, fixed(";")) > str_count(first_line, fixed(","))) {
    return(read_delim(path, delim = ";", show_col_types = FALSE, trim_ws = TRUE, escape_double = FALSE))
  }
  if (str_detect(first_line, "\\t")) return(read_tsv(path, show_col_types = FALSE))
  read_csv(path, show_col_types = FALSE)
}

read_table_with_header_search <- function(path) {
  if (is.na(path) || !file.exists(path)) return(tibble())
  lines <- readLines(path, n = 80, warn = FALSE)
  header_idx <- which(str_detect(lines, regex("(^|\\t)(source_file|ncrna_type|gene_id|miRNA)(\\t|$)", ignore_case = TRUE)))[1]
  if (is.na(header_idx)) return(read_auto(path))
  read_tsv(path, skip = header_idx - 1, show_col_types = FALSE)
}

choose_col <- function(df, candidates, fallback = NULL) {
  if (ncol(df) == 0) return(fallback)
  low <- str_to_lower(names(df))
  for (cand in candidates) {
    idx <- which(low == str_to_lower(cand))[1]
    if (!is.na(idx)) return(names(df)[idx])
  }
  for (cand in candidates) {
    idx <- which(str_detect(low, regex(cand, ignore_case = TRUE)))[1]
    if (!is.na(idx)) return(names(df)[idx])
  }
  fallback
}

pull_chr <- function(df, candidates, default = NA_character_) {
  col <- choose_col(df, candidates)
  if (is.null(col)) return(rep(default, nrow(df)))
  as.character(df[[col]])
}

pull_num <- function(df, candidates, default = NA_real_) {
  col <- choose_col(df, candidates)
  if (is.null(col)) return(rep(default, nrow(df)))
  suppressWarnings(as.numeric(df[[col]]))
}

norm_mir <- function(x) {
  x %>% as.character() %>% str_trim() %>% str_to_lower() %>% str_replace_all("\\s+", "")
}

norm_gene <- function(x) {
  x %>% as.character() %>% str_trim() %>% str_to_upper() %>% str_replace_all("\\s+", "")
}

canonical_direction <- function(x, logfc = NULL) {
  y <- str_to_upper(str_trim(as.character(x)))
  if (is.null(logfc)) logfc <- rep(NA_real_, length(y))
  logfc <- suppressWarnings(as.numeric(logfc))
  case_when(
    y %in% c("UP", "UPREGULATED", "POSITIVE") ~ "UP",
    y %in% c("DOWN", "DOWNREGULATED", "NEGATIVE") ~ "DOWN",
    !is.na(logfc) & logfc > 0 ~ "UP",
    !is.na(logfc) & logfc < 0 ~ "DOWN",
    TRUE ~ NA_character_
  )
}

canonical_locus <- function(x) {
  x <- as.character(x)
  m <- str_match(x, "([A-Za-z0-9_.-]+)[:_]([0-9]+)[|_]([0-9]+)")
  if_else(is.na(m[, 1]), NA_character_, paste0(m[, 2], ":", m[, 3], "|", m[, 4]))
}

collapse_chr <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & as.character(x) != ""])))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

collapse_num <- function(x, digits = 3) {
  x <- sort(unique(round(suppressWarnings(as.numeric(x[!is.na(x)])), digits)))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

first_non_missing <- function(x, default = NA) {
  z <- x[!is.na(x)]
  if (length(z) == 0) default else z[1]
}

collapse_direction <- function(x) {
  z <- unique(canonical_direction(x))
  z <- z[!is.na(z)]
  if (length(z) == 0) NA_character_ else if (length(z) == 1) z else "MIXED_VARIANTS"
}

save_plot <- function(p, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = height, bg = "white", useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = width, height = height, dpi = png_dpi, bg = "white")
}

empty_lnc_triplets <- function() tibble(
  lnc_source_file = character(), lncRNA_gene_id = character(), lncRNA_transcript_id = character(),
  lncRNA_label = character(), lncRNA_phase = character(), timepoint = character(),
  lncRNA_direction = character(), lncRNA_MFE = numeric(), lncRNA_interaction_tier = character(),
  miRNA = character(), miRNA_key = character(), miRNA_direction = character(),
  mRNA_gene = character(), gene_key = character(), mRNA_direction = character(),
  mRNA_log2FC = numeric(), mRNA_padj = numeric(), mRNA_phase = character(),
  target_region = character(), target_support_file = character(), target_score = numeric(),
  coherent_triplet_original = logical(), lnc_triplet_input_scope = character()
)

empty_targets <- function() tibble(
  timepoint = character(), miRNA = character(), miRNA_key = character(),
  mRNA_gene = character(), gene_key = character(), target_region = character(),
  target_score = numeric(), target_source_file = character(), target_source_class = character()
)

# 1. IntaRNA all-record circRNA-miRNA duplexes
interaction_specs <- tribble(
  ~timepoint, ~circ_direction_file, ~miRNA_direction_file, ~candidate_names,
  "6h",  "UP",   "DOWN", list(c("Interact_circUP_miDOWN_6h_all_modeM(1).csv", "Interact_circUP_miDOWN_6h_all_modeM.csv")),
  "6h",  "DOWN", "UP",   list(c("Interact_circDOWN_miUP_6h_all_modeM(1).csv", "Interact_circDOWN_miUP_6h_all_modeM.csv")),
  "24h", "UP",   "DOWN", list(c("Interact_circUP_miDOWN_24h_all_modeM(1).csv", "Interact_circUP_miDOWN_24h_all_modeM.csv")),
  "24h", "DOWN", "UP",   list(c("Interact_circDOWN_miUP_24h_all_modeM(1).csv", "Interact_circDOWN_miUP_24h_all_modeM.csv"))
) %>%
  mutate(path = map_chr(candidate_names, find_first))

if (all(is.na(interaction_specs$path))) stop("No circRNA all-record IntaRNA files were found.")

read_intarna <- function(path, timepoint, circ_direction, miRNA_direction) {
  if (is.na(path)) return(tibble())
  x <- read_auto(path)
  required <- c("id1", "start1", "end1", "id2", "start2", "end2", "E")
  if (!all(required %in% names(x))) stop("Unrecognized IntaRNA columns in ", basename(path))
  tibble(
    source_file = basename(path),
    timepoint = timepoint,
    circ_direction_file = circ_direction,
    miRNA_direction_file = miRNA_direction,
    circ_id1_raw = as.character(x$id1),
    circ_locus = canonical_locus(x$id1),
    circ_start = suppressWarnings(as.integer(x$start1)),
    circ_end = suppressWarnings(as.integer(x$end1)),
    miRNA = as.character(x$id2),
    miRNA_key = norm_mir(x$id2),
    miRNA_start = suppressWarnings(as.integer(x$start2)),
    miRNA_end = suppressWarnings(as.integer(x$end2)),
    subseqDP = if ("subseqDP" %in% names(x)) as.character(x$subseqDP) else NA_character_,
    hybridDP = if ("hybridDP" %in% names(x)) as.character(x$hybridDP) else NA_character_,
    E = suppressWarnings(as.numeric(x$E))
  ) %>%
    filter(!is.na(circ_locus), !is.na(miRNA_key), miRNA_key != "", !is.na(E))
}

intarna_records <- pmap_dfr(
  interaction_specs,
  function(timepoint, circ_direction_file, miRNA_direction_file, candidate_names, path) {
    read_intarna(path, timepoint, circ_direction_file, miRNA_direction_file)
  }
)

hc_records <- intarna_records %>% filter(E <= energy_cutoff)

hc_duplexes <- hc_records %>%
  group_by(timepoint, circ_locus, miRNA_key, circ_direction_file, miRNA_direction_file) %>%
  summarise(
    miRNA = first(miRNA),
    best_MFE = min(E, na.rm = TRUE),
    n_all_records = n(),
    n_unique_binding_sites = n_distinct(circ_start, circ_end, miRNA_start, miRNA_end),
    circ_binding_coordinates = collapse_chr(paste0(circ_start, "-", circ_end)),
    miRNA_binding_coordinates = collapse_chr(paste0(miRNA_start, "-", miRNA_end)),
    all_MFE_values = collapse_num(E),
    source_files = collapse_chr(source_file),
    .groups = "drop"
  ) %>%
  mutate(
    duplex_id = paste(timepoint, circ_locus, miRNA_key, sep = "__"),
    expected_mRNA_direction = if_else(miRNA_direction_file == "UP", "DOWN", "UP")
  )

# 2. circRNA locus, DE direction, host annotation and candidate variants
circ6_path <- find_first(c("06_circRNA_sequences_de6_DE_complete_table.txt"))
circ24_path <- find_first(c("07_circRNA_sequences_de24_DE_complete_table.txt"))
if (is.na(circ6_path) || is.na(circ24_path)) stop("Both circRNA DE complete tables are required.")

read_circ_de <- function(path, timepoint) {
  x <- read_auto(path)
  locus_col <- choose_col(x, c("circ_id", "circRNA_ID", "circRNA_ID2"), names(x)[1])
  id2_col <- choose_col(x, c("circRNA_ID2"))
  variant_col <- choose_col(x, c("circRNA_id_AS"))
  deid_col <- choose_col(x, c("circRNA_DE_id"))
  length_col <- choose_col(x, c("length"))
  direction_col <- choose_col(x, c("regulations", "direction"))
  lfc_col <- choose_col(x, c("log2FC", "log2FoldChange"))
  padj_col <- choose_col(x, c("padj", "FDR", "adj.P.Val"))
  pvalue_col <- choose_col(x, c("pvalue", "PValue"))
  tibble(
    timepoint = timepoint,
    circ_locus = canonical_locus(x[[locus_col]]),
    circ_annotation = if (!is.null(id2_col)) as.character(x[[id2_col]]) else NA_character_,
    variant_id = if (!is.null(variant_col)) as.character(x[[variant_col]]) else NA_character_,
    circRNA_DE_id = if (!is.null(deid_col)) as.character(x[[deid_col]]) else NA_character_,
    variant_length = if (!is.null(length_col)) suppressWarnings(as.numeric(x[[length_col]])) else NA_real_,
    circ_direction_variant = canonical_direction(
      if (!is.null(direction_col)) x[[direction_col]] else NA_character_,
      if (!is.null(lfc_col)) x[[lfc_col]] else NA_real_
    ),
    circ_log2FC_variant = if (!is.null(lfc_col)) suppressWarnings(as.numeric(x[[lfc_col]])) else NA_real_,
    circ_padj_variant = if (!is.null(padj_col)) suppressWarnings(as.numeric(x[[padj_col]])) else NA_real_,
    circ_pvalue_variant = if (!is.null(pvalue_col)) suppressWarnings(as.numeric(x[[pvalue_col]])) else NA_real_
  ) %>%
    filter(!is.na(circ_locus))
}

circ_variant_records <- bind_rows(
  read_circ_de(circ6_path, "6h"),
  read_circ_de(circ24_path, "24h")
)

circ_locus_map <- circ_variant_records %>%
  group_by(timepoint, circ_locus) %>%
  summarise(
    circ_direction = collapse_direction(circ_direction_variant),
    circ_log2FC = first_non_missing(circ_log2FC_variant, NA_real_),
    circ_padj = {
      z <- circ_padj_variant[!is.na(circ_padj_variant)]
      if (length(z) == 0) NA_real_ else min(z)
    },
    circ_pvalue = {
      z <- circ_pvalue_variant[!is.na(circ_pvalue_variant)]
      if (length(z) == 0) NA_real_ else min(z)
    },
    n_candidate_variants = n_distinct(variant_id[!is.na(variant_id) & variant_id != ""]),
    candidate_variant_ids = collapse_chr(variant_id),
    candidate_variant_lengths = collapse_num(variant_length, digits = 0),
    candidate_circRNA_DE_ids = collapse_chr(circRNA_DE_id),
    circ_annotations = collapse_chr(circ_annotation),
    .groups = "drop"
  )

hc_duplexes <- hc_duplexes %>%
  left_join(circ_locus_map, by = c("timepoint", "circ_locus")) %>%
  mutate(
    circ_mapping_status = case_when(
      is.na(circ_direction) ~ "unmapped_circRNA_locus",
      circ_direction == circ_direction_file ~ "mapped_direction_confirmed",
      circ_direction == "MIXED_VARIANTS" ~ "mapped_mixed_variant_directions",
      TRUE ~ "mapped_direction_discordant"
    ),
    variant_resolution_status = case_when(
      is.na(n_candidate_variants) ~ "unmapped",
      n_candidate_variants == 0 ~ "locus_mapped_variant_id_missing",
      n_candidate_variants == 1 ~ "single_candidate_variant",
      TRUE ~ "multiple_candidate_variants_unresolved"
    )
  )

# 3. Verify the selected miRNAs with the full miRNA DE files
mir6_path <- find_first(c("DE_miRNA_All_with_status_NE6(1).txt", "DE_miRNA_All_with_status_NE6.txt"))
mir24_path <- find_first(c("DE_miRNA_All_with_status_NE24(1).txt", "DE_miRNA_All_with_status_NE24.txt"))

read_mirna_de <- function(path, timepoint) {
  if (is.na(path)) return(tibble())
  x <- read_auto(path)
  mir_col <- choose_col(x, c("miRNA_id", "miRNA", "mirna"), names(x)[1])
  lfc <- pull_num(x, c("log2FoldChange", "log2FC"))
  padj <- pull_num(x, c("padj", "FDR", "adj.P.Val"))
  tibble(
    timepoint = timepoint,
    miRNA = as.character(x[[mir_col]]),
    miRNA_key = norm_mir(x[[mir_col]]),
    miRNA_log2FC = lfc,
    miRNA_padj = padj,
    miRNA_DE_direction = case_when(
      !is.na(padj) & padj < alpha & lfc >= mirna_lfc_cutoff ~ "UP",
      !is.na(padj) & padj < alpha & lfc <= -mirna_lfc_cutoff ~ "DOWN",
      TRUE ~ "NS"
    )
  ) %>%
    distinct(timepoint, miRNA_key, .keep_all = TRUE)
}

mirna_de <- bind_rows(
  read_mirna_de(mir6_path, "6h"),
  read_mirna_de(mir24_path, "24h")
)

hc_duplexes <- hc_duplexes %>%
  left_join(
    mirna_de %>% select(timepoint, miRNA_key, miRNA_log2FC, miRNA_padj, miRNA_DE_direction),
    by = c("timepoint", "miRNA_key")
  ) %>%
  mutate(
    miRNA_verification = case_when(
      is.na(miRNA_DE_direction) ~ "miRNA_not_mapped",
      miRNA_DE_direction == miRNA_direction_file ~ "DE_direction_confirmed",
      miRNA_DE_direction == "NS" ~ "not_DE_by_current_threshold",
      TRUE ~ "DE_direction_discordant"
    ),
    opposite_circ_miRNA_directions = circ_direction_file != miRNA_direction_file,
    duplex_passes_validation = best_MFE <= energy_cutoff &
      circ_mapping_status %in% c("mapped_direction_confirmed", "mapped_mixed_variant_directions") &
      miRNA_verification == "DE_direction_confirmed" &
      opposite_circ_miRNA_directions
  )

selected_miRNAs <- hc_duplexes %>%
  distinct(timepoint, miRNA, miRNA_key, miRNA_direction_file, miRNA_log2FC, miRNA_padj, miRNA_verification)

# 4. Read the complete lncRNA-miRNA-mRNA triplets and the previous priority/figure table
lnc_complete_path <- find_first(c(
  "COHERENT_lncRNA_miRNA_mRNA_triplets.csv",
  "COHERENT_lncRNA_miRNA_mRNA_triplets(1).csv"
))
lnc_priority_path <- find_first(c(
  "PRIORITY_ceRNA_triplets_best_per_ncRNA_miRNA_mRNA.csv",
  "PRIORITY_ceRNA_triplets_best_per_ncRNA_miRNA_mRNA(1).csv"
))
figure_priority_path <- find_first(c(
  "Fig8_priority_ceRNA_triplets(1).txt",
  "Fig8_priority_ceRNA_triplets.txt",
  "Figure7_priority_ceRNA_triplets.csv",
  "MAIN_FIGURE_ceRNA_triplets.csv"
))

standardize_lnc_triplets <- function(path, input_scope) {
  if (is.na(path)) return(empty_lnc_triplets())
  x <- if (str_detect(path, regex("\\.txt$", ignore_case = TRUE))) read_table_with_header_search(path) else read_auto(path)
  if (nrow(x) == 0) return(empty_lnc_triplets())

  timepoint <- pull_chr(x, c("timepoint", "time", "state"))
  miRNA <- pull_chr(x, c("miRNA", "mirna"))
  gene_symbol <- pull_chr(x, c("gene_symbol", "gene_symbol_miRWalk", "mRNA_label", "target_gene", "gene"))
  mrna_lfc <- pull_num(x, c("mRNA_log2FC", "mRNA_log2FoldChange", "target_log2FC"))
  mrna_padj <- pull_num(x, c("mRNA_padj", "target_padj", "padj_mRNA"))
  mrna_dir_raw <- pull_chr(x, c("mRNA_direction", "mRNA_dir", "target_direction"))
  coherent_raw <- str_to_upper(pull_chr(x, c("coherent_triplet", "coherent"), "TRUE"))

  tibble(
    lnc_source_file = pull_chr(x, c("source_file"), basename(path)),
    lncRNA_gene_id = pull_chr(x, c("gene_id", "lncRNA_gene_id")),
    lncRNA_transcript_id = pull_chr(x, c("transcript_id", "lncRNA_transcript_id")),
    lncRNA_label = pull_chr(x, c("ncRNA_label", "lncRNA_label", "gene_id", "transcript_id")),
    lncRNA_phase = pull_chr(x, c("phase", "lncRNA_phase")),
    timepoint = timepoint,
    lncRNA_direction = canonical_direction(pull_chr(x, c("ncRNA_direction", "lncRNA_direction"))),
    lncRNA_MFE = pull_num(x, c("E", "MFE", "lncRNA_MFE")),
    lncRNA_interaction_tier = pull_chr(x, c("interaction_tier", "tier")),
    miRNA = miRNA,
    miRNA_key = norm_mir(miRNA),
    miRNA_direction = canonical_direction(pull_chr(x, c("miRNA_direction", "miRNA_dir"))),
    mRNA_gene = gene_symbol,
    gene_key = norm_gene(gene_symbol),
    mRNA_direction = canonical_direction(mrna_dir_raw, mrna_lfc),
    mRNA_log2FC = mrna_lfc,
    mRNA_padj = mrna_padj,
    mRNA_phase = pull_chr(x, c("phase_mRNA", "mRNA_phase", "phase_target")),
    target_region = pull_chr(x, c("region", "target_region", "region_support")),
    target_support_file = pull_chr(x, c("support_file", "target_source", "target_support_file")),
    target_score = pull_num(x, c("target_score", "score")),
    coherent_triplet_original = coherent_raw %in% c("TRUE", "T", "1", "YES"),
    lnc_triplet_input_scope = input_scope
  ) %>%
    filter(!is.na(timepoint), timepoint != "", !is.na(miRNA_key), miRNA_key != "", !is.na(gene_key), gene_key != "") %>%
    distinct()
}

lnc_complete <- standardize_lnc_triplets(lnc_complete_path, "complete_coherent_lncRNA_triplets")
lnc_priority <- standardize_lnc_triplets(lnc_priority_path, "previous_priority_table")
figure_priority <- standardize_lnc_triplets(figure_priority_path, "previous_figure_priority_table")

if (nrow(lnc_complete) == 0) {
  if (nrow(lnc_priority) > 0) {
    warning("Complete lncRNA triplet file was not found; using the previous priority table as the lncRNA comparison universe.")
    lnc_complete <- lnc_priority
  } else if (nrow(figure_priority) > 0) {
    warning("Complete lncRNA triplet file was not found; using the previous figure-priority table as the lncRNA comparison universe.")
    lnc_complete <- figure_priority
  } else {
    warning("No lncRNA triplet table was found. The script will continue with miRWalk-only reconstruction.")
  }
}

lnc_selected <- bind_rows(empty_lnc_triplets(), lnc_complete) %>%
  semi_join(selected_miRNAs, by = c("timepoint", "miRNA_key")) %>%
  filter(coherent_triplet_original | is.na(coherent_triplet_original))

priority_pair_keys <- bind_rows(empty_lnc_triplets(), lnc_priority, figure_priority) %>%
  transmute(timepoint, miRNA_key, gene_key, in_previous_priority = TRUE) %>%
  distinct()

figure_pair_keys <- bind_rows(empty_lnc_triplets(), figure_priority) %>%
  transmute(timepoint, miRNA_key, gene_key, in_previous_figure_priority = TRUE) %>%
  distinct()

# 5. Reconstruct circRNA triplets from existing lncRNA triplets using only miRNA + timepoint
cross_ncRNA_modules <- hc_duplexes %>%
  inner_join(lnc_selected, by = c("timepoint", "miRNA_key"), suffix = c("_circ", "_lnc")) %>%
  mutate(
    miRNA = coalesce(miRNA_circ, miRNA_lnc),
    circ_expected_mRNA_direction = expected_mRNA_direction,
    circ_mRNA_direction_coherent = mRNA_direction == circ_expected_mRNA_direction,
    strict_circ_triplet_from_lnc = duplex_passes_validation & circ_mRNA_direction_coherent,
    cross_module_id = paste(timepoint, circ_locus, lncRNA_label, miRNA, mRNA_gene, sep = "__")
  ) %>%
  left_join(priority_pair_keys, by = c("timepoint", "miRNA_key", "gene_key")) %>%
  left_join(figure_pair_keys, by = c("timepoint", "miRNA_key", "gene_key")) %>%
  mutate(
    in_previous_priority = replace_na(in_previous_priority, FALSE),
    in_previous_figure_priority = replace_na(in_previous_figure_priority, FALSE),
    cross_module_class = case_when(
      strict_circ_triplet_from_lnc & in_previous_figure_priority ~ "Strict_cross_ncRNA_axis_already_in_figure_priority",
      strict_circ_triplet_from_lnc & in_previous_priority ~ "Strict_cross_ncRNA_axis_in_previous_priority",
      strict_circ_triplet_from_lnc ~ "Strict_cross_ncRNA_axis_from_complete_lncRNA_triplets",
      TRUE ~ "Shared_miRNA_but_not_directionally_coherent_for_circRNA"
    )
  )

# 6. Read raw miRWalk targets for the selected circRNA-associated miRNAs
raw_mirwalk_files <- find_pattern("^miRWalk_miRNA_Targets_.*\\.csv$")

infer_timepoint_from_filename <- function(filename) {
  case_when(
    str_detect(filename, regex("24h|UP24|DOWN24|_24", ignore_case = TRUE)) ~ "24h",
    str_detect(filename, regex("6h|UP6|DOWN6|_6", ignore_case = TRUE)) ~ "6h",
    TRUE ~ NA_character_
  )
}

read_mirwalk <- function(path) {
  x <- read_auto(path)
  if (nrow(x) == 0) return(empty_targets())
  mir_col <- choose_col(x, c("miRNA", "mirna", "miRNA_name", "mature_miRNA"))
  gene_col <- choose_col(x, c("gene_symbol", "genesymbol", "target_gene", "gene_name", "Symbol", "gene"))
  if (is.null(mir_col) || is.null(gene_col)) return(empty_targets())
  fn <- basename(path)
  tibble(
    timepoint = infer_timepoint_from_filename(fn),
    miRNA = as.character(x[[mir_col]]),
    miRNA_key = norm_mir(x[[mir_col]]),
    mRNA_gene = as.character(x[[gene_col]]),
    gene_key = norm_gene(x[[gene_col]]),
    target_region = case_when(
      str_detect(fn, regex("3UTR", ignore_case = TRUE)) ~ "3UTR",
      str_detect(fn, regex("5UTR", ignore_case = TRUE)) ~ "5UTR",
      str_detect(fn, regex("CDS", ignore_case = TRUE)) ~ "CDS",
      TRUE ~ pull_chr(x, c("region", "binding_region"), "Unknown")
    ),
    target_score = pull_num(x, c("score", "binding_score", "probability", "energy")),
    target_source_file = fn,
    target_source_class = "raw_miRWalk"
  ) %>%
    filter(!is.na(timepoint), !is.na(miRNA_key), miRNA_key != "", !is.na(gene_key), gene_key != "")
}

raw_mirwalk_targets <- if (length(raw_mirwalk_files) > 0) {
  map_dfr(raw_mirwalk_files, read_mirwalk)
} else {
  empty_targets()
}

# Fallback target evidence from previous KEGG edge files and lncRNA triplets
supp_edge_files <- find_pattern("^SUPP_edges_(6h|24h)_miRNA(Up|Down)_mRNA(Up|Down)_KEGG_phase\\.csv$")

read_supp_edges <- function(path) {
  x <- read_auto(path)
  if (nrow(x) == 0) return(empty_targets())
  tibble(
    timepoint = pull_chr(x, c("time", "timepoint")),
    miRNA = pull_chr(x, c("miRNA", "mirna")),
    miRNA_key = norm_mir(pull_chr(x, c("miRNA_key", "miRNA", "mirna"))),
    mRNA_gene = pull_chr(x, c("gene_key", "gene_symbol", "gene")),
    gene_key = norm_gene(pull_chr(x, c("gene_key", "gene_symbol", "gene"))),
    target_region = pull_chr(x, c("region_support", "region")),
    target_score = NA_real_,
    target_source_file = basename(path),
    target_source_class = "previous_DE_coherent_KEGG_edge"
  ) %>%
    filter(!is.na(timepoint), !is.na(miRNA_key), miRNA_key != "", !is.na(gene_key), gene_key != "")
}

supp_targets <- if (length(supp_edge_files) > 0) map_dfr(supp_edge_files, read_supp_edges) else empty_targets()

lnc_target_fallback <- bind_rows(empty_lnc_triplets(), lnc_complete) %>%
  transmute(
    timepoint, miRNA, miRNA_key, mRNA_gene, gene_key,
    target_region, target_score,
    target_source_file = coalesce(target_support_file, lnc_source_file),
    target_source_class = "existing_lncRNA_triplet"
  )

all_target_evidence <- bind_rows(empty_targets(), raw_mirwalk_targets, supp_targets, lnc_target_fallback) %>%
  semi_join(selected_miRNAs, by = c("timepoint", "miRNA_key")) %>%
  distinct()

candidate_miRNA_mRNA <- all_target_evidence %>%
  group_by(timepoint, miRNA_key, gene_key) %>%
  summarise(
    miRNA = first_non_missing(miRNA, NA_character_),
    mRNA_gene = first_non_missing(mRNA_gene, NA_character_),
    target_regions = collapse_chr(target_region),
    target_scores = collapse_num(target_score),
    target_source_files = collapse_chr(target_source_file),
    target_source_classes = collapse_chr(target_source_class),
    supported_by_raw_miRWalk = any(target_source_class == "raw_miRWalk"),
    supported_by_existing_lncRNA_triplet = any(target_source_class == "existing_lncRNA_triplet"),
    supported_by_previous_KEGG_edge = any(target_source_class == "previous_DE_coherent_KEGG_edge"),
    .groups = "drop"
  ) %>%
  left_join(priority_pair_keys, by = c("timepoint", "miRNA_key", "gene_key")) %>%
  left_join(figure_pair_keys, by = c("timepoint", "miRNA_key", "gene_key")) %>%
  mutate(
    in_previous_priority = replace_na(in_previous_priority, FALSE),
    in_previous_figure_priority = replace_na(in_previous_figure_priority, FALSE)
  )

# 7. Full mRNA DESeq2 results for DE-coherent, non-DE, discordant and unmapped classification
annotation_path <- find_first(c("hisat2_6_24_meta_ann.csv"))

read_gene_annotation <- function(path) {
  if (is.na(path)) return(tibble(gene_id = character(), gene_symbol_annotation = character()))
  x <- read_auto(path)
  id_col <- choose_col(x, c("ensembl_gene_id", "gene_id"))
  symbol_col <- choose_col(x, c("external_gene_name", "gene_name", "gene_symbol", "symbol"))
  if (is.null(id_col) || is.null(symbol_col)) return(tibble(gene_id = character(), gene_symbol_annotation = character()))
  tibble(
    gene_id = as.character(x[[id_col]]),
    gene_symbol_annotation = as.character(x[[symbol_col]])
  ) %>%
    filter(!is.na(gene_id), gene_id != "") %>%
    distinct(gene_id, .keep_all = TRUE)
}

gene_annotation <- read_gene_annotation(annotation_path)

mRNA6_full_path <- find_first(c("res_NE_vs_Ctrl_6h.csv", "res_NE_vs_Ctrl_6h.txt", "NE_vs_Ctrl_6h_mRNA_full_with_NS.txt"))
mRNA24_full_path <- find_first(c("res_NE_vs_Ctrl_24h.csv", "res_NE_vs_Ctrl_24h.txt", "NE_vs_Ctrl_24h_mRNA_full_with_NS.txt"))
mRNA6_de_path <- find_first(c("NE_vs_Ctrl_6h_mRNA_DE(2).txt", "NE_vs_Ctrl_6h_mRNA_DE.txt"))
mRNA24_de_path <- find_first(c("NE_vs_Ctrl_24h_mRNA_DE(2).txt", "NE_vs_Ctrl_24h_mRNA_DE.txt"))

read_mrna_results <- function(path, timepoint, is_full_table) {
  if (is.na(path)) return(tibble())
  x <- read_auto(path)
  if (nrow(x) == 0) return(tibble())
  id_col <- choose_col(x, c("gene_id", "ensembl_gene_id", "...1", "X"), names(x)[1])
  symbol_col <- choose_col(x, c("gene_name", "gene_symbol", "external_gene_name", "symbol"))
  gene_id <- as.character(x[[id_col]])
  gene_symbol_input <- if (!is.null(symbol_col)) as.character(x[[symbol_col]]) else rep(NA_character_, nrow(x))
  lfc <- pull_num(x, c("log2FoldChange", "log2FC"))
  padj <- pull_num(x, c("padj", "FDR", "adj.P.Val"))
  base_mean <- pull_num(x, c("baseMean", "mean_normalized_count", "mean_count"))
  direction_input <- canonical_direction(pull_chr(x, c("direction", "regulations")), lfc)

  tibble(
    gene_id = gene_id,
    gene_symbol_input = gene_symbol_input,
    mRNA_log2FC = lfc,
    mRNA_padj = padj,
    mRNA_baseMean = base_mean,
    direction_input = direction_input
  ) %>%
    left_join(gene_annotation, by = "gene_id") %>%
    mutate(
      timepoint = timepoint,
      mRNA_gene = coalesce(gene_symbol_input, gene_symbol_annotation, gene_id),
      gene_key = norm_gene(mRNA_gene),
      mRNA_DE_direction = case_when(
        !is.na(mRNA_padj) & mRNA_padj < alpha & mRNA_log2FC >= mrna_lfc_cutoff ~ "UP",
        !is.na(mRNA_padj) & mRNA_padj < alpha & mRNA_log2FC <= -mrna_lfc_cutoff ~ "DOWN",
        is_full_table ~ "NS",
        !is.na(direction_input) ~ direction_input,
        TRUE ~ NA_character_
      ),
      mRNA_input_scope = if_else(is_full_table, "full_DESeq2_results", "DE_only_results")
    ) %>%
    filter(!is.na(gene_key), gene_key != "") %>%
    distinct(timepoint, gene_key, .keep_all = TRUE) %>%
    select(timepoint, gene_id, mRNA_gene, gene_key, mRNA_log2FC, mRNA_padj, mRNA_baseMean, mRNA_DE_direction, mRNA_input_scope)
}

mRNA_results <- bind_rows(
  read_mrna_results(if_else(!is.na(mRNA6_full_path), mRNA6_full_path, mRNA6_de_path), "6h", !is.na(mRNA6_full_path)),
  read_mrna_results(if_else(!is.na(mRNA24_full_path), mRNA24_full_path, mRNA24_de_path), "24h", !is.na(mRNA24_full_path))
)

# 8. Reconstruct all putative circRNA-miRNA-mRNA relationships
all_putative <- hc_duplexes %>%
  left_join(candidate_miRNA_mRNA, by = c("timepoint", "miRNA_key"), suffix = c("_circ", "_target")) %>%
  left_join(mRNA_results, by = c("timepoint", "gene_key"), suffix = c("_target", "_expression")) %>%
  mutate(
    miRNA = coalesce(miRNA_circ, miRNA_target),
    supported_by_raw_miRWalk = replace_na(supported_by_raw_miRWalk, FALSE),
    supported_by_existing_lncRNA_triplet = replace_na(supported_by_existing_lncRNA_triplet, FALSE),
    supported_by_previous_KEGG_edge = replace_na(supported_by_previous_KEGG_edge, FALSE),
    in_previous_priority = replace_na(in_previous_priority, FALSE),
    in_previous_figure_priority = replace_na(in_previous_figure_priority, FALSE),
    mRNA_gene_final = coalesce(mRNA_gene_expression, mRNA_gene_target),
    strict_direction_coherence = mRNA_DE_direction == expected_mRNA_direction,
    target_evidence_class = case_when(
      is.na(gene_key) ~ "No_candidate_miRNA_mRNA_target_found",
      !is.na(mRNA_DE_direction) & mRNA_DE_direction == expected_mRNA_direction ~ "DE_coherent",
      mRNA_DE_direction %in% c("UP", "DOWN") & mRNA_DE_direction != expected_mRNA_direction ~ "DE_discordant",
      mRNA_DE_direction == "NS" ~ "Expressed_nonDE",
      is.na(mRNA_DE_direction) & supported_by_existing_lncRNA_triplet ~ "Existing_lncRNA_triplet_target_not_remapped",
      is.na(mRNA_DE_direction) & supported_by_previous_KEGG_edge ~ "Previous_DE_coherent_target_not_remapped",
      TRUE ~ "Not_mapped_to_expression"
    ),
    strict_DE_coherent_circRNA_triplet = duplex_passes_validation & target_evidence_class == "DE_coherent",
    relationship_origin = case_when(
      strict_DE_coherent_circRNA_triplet & supported_by_existing_lncRNA_triplet ~ "Shared_with_existing_lncRNA_triplet",
      strict_DE_coherent_circRNA_triplet & supported_by_raw_miRWalk ~ "New_strict_circRNA_triplet_from_miRWalk",
      supported_by_existing_lncRNA_triplet ~ "Existing_lncRNA_target_without_strict_circRNA_support",
      supported_by_raw_miRWalk ~ "Additional_miRWalk_target_context",
      TRUE ~ "Fallback_target_context"
    ),
    recommended_placement = case_when(
      strict_DE_coherent_circRNA_triplet & in_previous_figure_priority ~ "Potential_main_figure_extension_or_cross_ncRNA_inset",
      strict_DE_coherent_circRNA_triplet & supported_by_existing_lncRNA_triplet ~ "Supplementary_cross_ncRNA_convergence_module",
      strict_DE_coherent_circRNA_triplet ~ "Supplementary_new_circRNA_triplet",
      target_evidence_class == "Expressed_nonDE" ~ "Supplementary_predicted_target_context_only",
      TRUE ~ "Audit_or_supplementary_table_only"
    ),
    putative_triplet_id = paste(timepoint, circ_locus, miRNA_key, gene_key, sep = "__")
  )

strict_triplets <- all_putative %>% filter(strict_DE_coherent_circRNA_triplet)
shared_strict_triplets <- strict_triplets %>% filter(supported_by_existing_lncRNA_triplet)
new_strict_triplets <- strict_triplets %>% filter(!supported_by_existing_lncRNA_triplet)
expressed_nonDE_targets <- all_putative %>% filter(target_evidence_class == "Expressed_nonDE")

# 9. Audit each filtering stage
audit <- all_putative %>%
  mutate(
    failure_reason = case_when(
      best_MFE > energy_cutoff ~ "MFE_above_cutoff",
      circ_mapping_status == "unmapped_circRNA_locus" ~ "circRNA_locus_not_mapped",
      circ_mapping_status == "mapped_direction_discordant" ~ "circRNA_direction_discordant",
      miRNA_verification == "miRNA_not_mapped" ~ "miRNA_not_found_in_DE_table",
      miRNA_verification == "not_DE_by_current_threshold" ~ "miRNA_not_DE_by_current_threshold",
      miRNA_verification == "DE_direction_discordant" ~ "miRNA_direction_discordant",
      is.na(gene_key) ~ "no_candidate_miRNA_mRNA_target",
      target_evidence_class == "Not_mapped_to_expression" ~ "target_not_mapped_to_mRNA_results",
      target_evidence_class == "Expressed_nonDE" ~ "target_expressed_but_not_DE",
      target_evidence_class == "DE_discordant" ~ "target_DE_but_direction_discordant",
      strict_DE_coherent_circRNA_triplet ~ "PASS_strict_DE_coherent_triplet",
      TRUE ~ target_evidence_class
    )
  )

# 10. circRNA temporal overlap at genomic-locus level
circ_temporal <- circ_locus_map %>%
  select(timepoint, circ_locus, circ_direction) %>%
  distinct() %>%
  pivot_wider(names_from = timepoint, values_from = circ_direction) %>%
  mutate(
    `6h` = replace_na(`6h`, "Not_in_DE_table"),
    `24h` = replace_na(`24h`, "Not_in_DE_table"),
    temporal_class = case_when(
      `6h` %in% c("UP", "DOWN") & `24h` == "Not_in_DE_table" ~ paste0("6h_only_", `6h`),
      `6h` == "Not_in_DE_table" & `24h` %in% c("UP", "DOWN") ~ paste0("24h_only_", `24h`),
      `6h` == `24h` & `6h` %in% c("UP", "DOWN") ~ paste0("Sustained_", `6h`),
      `6h` %in% c("UP", "DOWN") & `24h` %in% c("UP", "DOWN") & `6h` != `24h` ~ "Direction_reversal",
      str_detect(`6h`, "MIXED") | str_detect(`24h`, "MIXED") ~ "Mixed_variant_directions",
      TRUE ~ "Other"
    )
  )

overlap_counts <- circ_temporal %>% count(direction_6h = `6h`, direction_24h = `24h`, name = "n_circRNA_loci")
transition_counts <- circ_temporal %>% count(temporal_class, name = "n_circRNA_loci") %>% arrange(desc(n_circRNA_loci))

# 11. Cytoscape-ready cross-ncRNA nodes and edges
cyto_cross <- cross_ncRNA_modules %>% filter(strict_circ_triplet_from_lnc)

cyto_nodes <- bind_rows(
  cyto_cross %>% transmute(node_id = paste0("circ:", circ_locus), label = circ_locus, node_type = "circRNA", timepoint, direction = circ_direction_file, log2FC = circ_log2FC),
  cyto_cross %>% transmute(node_id = paste0("lnc:", lncRNA_label), label = lncRNA_label, node_type = "lncRNA", timepoint, direction = lncRNA_direction, log2FC = NA_real_),
  cyto_cross %>% transmute(node_id = paste0("mir:", miRNA_key), label = miRNA, node_type = "miRNA", timepoint, direction = miRNA_direction_file, log2FC = miRNA_log2FC),
  cyto_cross %>% transmute(node_id = paste0("mrna:", gene_key), label = mRNA_gene, node_type = "mRNA", timepoint, direction = mRNA_direction, log2FC = mRNA_log2FC)
) %>% distinct(node_id, timepoint, .keep_all = TRUE)

cyto_edges <- bind_rows(
  cyto_cross %>% transmute(source = paste0("circ:", circ_locus), target = paste0("mir:", miRNA_key), edge_type = "predicted_circRNA_miRNA", timepoint, MFE = best_MFE, support = paste0(n_all_records, " all-record hit(s)")),
  cyto_cross %>% transmute(source = paste0("lnc:", lncRNA_label), target = paste0("mir:", miRNA_key), edge_type = "predicted_lncRNA_miRNA", timepoint, MFE = lncRNA_MFE, support = lnc_source_file),
  cyto_cross %>% transmute(source = paste0("mir:", miRNA_key), target = paste0("mrna:", gene_key), edge_type = "predicted_miRNA_mRNA", timepoint, MFE = target_score, support = target_support_file)
) %>% distinct()

# 12. Figures
p_overlap <- ggplot(overlap_counts, aes(direction_24h, direction_6h, fill = n_circRNA_loci)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = comma(n_circRNA_loci)), fontface = "bold", size = 4.5) +
  scale_fill_gradient(low = "grey95", high = "grey25", labels = comma) +
  labs(title = "Directional overlap of circRNA loci", subtitle = "Locus-level classification from the 6 h and 24 h circRNA input tables", x = "24 h state", y = "6 h state", fill = "Loci") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.grid = element_blank())
save_plot(p_overlap, "FIG1_circRNA_directional_overlap_heatmap", 7.5, 6.2)

p_transition <- ggplot(transition_counts, aes(fct_reorder(temporal_class, n_circRNA_loci), n_circRNA_loci)) +
  geom_col(fill = "grey35") +
  geom_text(aes(label = comma(n_circRNA_loci)), hjust = -0.15, size = 3.8) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "circRNA temporal classes", x = NULL, y = "Number of circRNA loci") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank(), plot.margin = margin(5.5, 35, 5.5, 5.5))
save_plot(p_transition, "FIG2_circRNA_temporal_classes", 9, 6)

if (nrow(hc_duplexes) > 0) {
  p_duplex <- hc_duplexes %>%
    mutate(module_label = paste0(circ_locus, " | ", miRNA), timepoint = factor(timepoint, c("6h", "24h"))) %>%
    ggplot(aes(timepoint, fct_rev(factor(module_label)), size = n_all_records, color = best_MFE)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(5, 11)) +
    scale_color_gradient(low = "grey10", high = "grey70") +
    labs(title = "High-confidence circRNA-miRNA duplexes", subtitle = paste0("IntaRNA MFE <= ", energy_cutoff, " kcal/mol"), x = "State", y = NULL, size = "All-record hits", color = "Best MFE") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.grid.minor = element_blank())
  save_plot(p_duplex, "FIG3_high_confidence_circRNA_miRNA_duplexes", 9.5, 5.5)
}

if (nrow(all_putative) > 0) {
  evidence_counts <- all_putative %>%
    filter(!is.na(gene_key)) %>%
    count(timepoint, circ_locus, miRNA, target_evidence_class, name = "n_targets") %>%
    mutate(module = paste(timepoint, circ_locus, miRNA, sep = " | "))

  p_evidence <- ggplot(evidence_counts, aes(fct_reorder(module, n_targets, .fun = sum), n_targets, fill = target_evidence_class)) +
    geom_col() +
    coord_flip() +
    labs(title = "Evidence classes for circRNA-associated miRNA targets", x = NULL, y = "Candidate mRNA targets", fill = "Evidence class") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank())
  save_plot(p_evidence, "FIG4_circRNA_miRNA_target_evidence", 11, 6.5)
}

if (nrow(cross_ncRNA_modules) > 0) {
  cross_counts <- cross_ncRNA_modules %>%
    filter(strict_circ_triplet_from_lnc) %>%
    distinct(timepoint, circ_locus, lncRNA_label, miRNA_key, gene_key) %>%
    count(timepoint, circ_locus, miRNA, cross_module_class, name = "n_lncRNA_mRNA_axes") %>%
    mutate(module = paste(timepoint, circ_locus, miRNA, sep = " | "))

  if (nrow(cross_counts) > 0) {
    p_cross <- ggplot(cross_counts, aes(fct_reorder(module, n_lncRNA_mRNA_axes, .fun = sum), n_lncRNA_mRNA_axes, fill = cross_module_class)) +
      geom_col() +
      coord_flip() +
      labs(title = "Cross-ncRNA convergence through shared miRNA-mRNA axes", x = NULL, y = "Distinct lncRNA-mRNA axes", fill = "Classification") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank())
    save_plot(p_cross, "FIG5_cross_ncRNA_convergence_summary", 11, 6.5)
  }
}

# 13. Export tables
input_manifest <- tribble(
  ~input_role, ~path, ~required,
  "IntaRNA circUP-miDOWN 6h", interaction_specs$path[1], TRUE,
  "IntaRNA circDOWN-miUP 6h", interaction_specs$path[2], TRUE,
  "IntaRNA circUP-miDOWN 24h", interaction_specs$path[3], TRUE,
  "IntaRNA circDOWN-miUP 24h", interaction_specs$path[4], TRUE,
  "circRNA DE table 6h", circ6_path, TRUE,
  "circRNA DE table 24h", circ24_path, TRUE,
  "miRNA full DE 6h", mir6_path, TRUE,
  "miRNA full DE 24h", mir24_path, TRUE,
  "complete lncRNA triplets", lnc_complete_path, FALSE,
  "previous lncRNA priority table", lnc_priority_path, FALSE,
  "previous figure-priority table", figure_priority_path, FALSE,
  "mRNA full results 6h", mRNA6_full_path, FALSE,
  "mRNA full results 24h", mRNA24_full_path, FALSE,
  "mRNA DE-only fallback 6h", mRNA6_de_path, FALSE,
  "mRNA DE-only fallback 24h", mRNA24_de_path, FALSE,
  "gene annotation", annotation_path, FALSE
) %>%
  mutate(found = !is.na(path) & file.exists(path), filename = if_else(found, basename(path), NA_character_))

write_csv(input_manifest, file.path(out_dir, "00_INPUT_manifest.csv"))
write_csv(intarna_records, file.path(out_dir, "01_ALL_IntaRNA_circRNA_miRNA_records.csv"))
write_csv(hc_records, file.path(out_dir, "02_HC_IntaRNA_records_MFE_le_minus20.csv"))
write_csv(hc_duplexes, file.path(out_dir, "03_HC_circRNA_miRNA_duplexes_locus_collapsed.csv"))
write_csv(circ_variant_records, file.path(out_dir, "04_circRNA_variant_level_annotation.csv"))
write_csv(circ_locus_map, file.path(out_dir, "05_circRNA_locus_annotation_and_candidate_variants.csv"))
write_csv(selected_miRNAs, file.path(out_dir, "06_selected_circRNA_associated_miRNAs.csv"))
write_csv(lnc_selected, file.path(out_dir, "07_existing_lncRNA_triplets_for_selected_miRNAs.csv"))
write_csv(cross_ncRNA_modules, file.path(out_dir, "08_CROSS_ncRNA_convergent_modules.csv"))
write_csv(raw_mirwalk_targets %>% semi_join(selected_miRNAs, by = c("timepoint", "miRNA_key")), file.path(out_dir, "09_raw_miRWalk_targets_for_selected_miRNAs.csv"))
write_csv(candidate_miRNA_mRNA, file.path(out_dir, "10_candidate_miRNA_mRNA_target_universe.csv"))
write_csv(all_putative, file.path(out_dir, "11_ALL_putative_circRNA_miRNA_mRNA_relationships.csv"))
write_csv(strict_triplets, file.path(out_dir, "12_STRICT_DE_coherent_circRNA_triplets.csv"))
write_csv(shared_strict_triplets, file.path(out_dir, "13_SHARED_with_lncRNA_triplet_axes.csv"))
write_csv(new_strict_triplets, file.path(out_dir, "14_NEW_strict_circRNA_triplets_not_in_lncRNA_analysis.csv"))
write_csv(expressed_nonDE_targets, file.path(out_dir, "15_EXPRESSED_nonDE_circRNA_miRNA_targets.csv"))
write_csv(audit, file.path(out_dir, "16_AUDIT_circRNA_target_filtering.csv"))
write_csv(cyto_nodes, file.path(out_dir, "17_CYTOSCAPE_cross_ncRNA_nodes.csv"))
write_csv(cyto_edges, file.path(out_dir, "18_CYTOSCAPE_cross_ncRNA_edges.csv"))
write_csv(circ_temporal, file.path(out_dir, "19_circRNA_temporal_direction_by_locus.csv"))
write_csv(overlap_counts, file.path(out_dir, "20_circRNA_directional_overlap_counts.csv"))
write_csv(transition_counts, file.path(out_dir, "21_circRNA_temporal_class_counts.csv"))

summary_table <- tibble(
  metric = c(
    "IntaRNA records read",
    "IntaRNA records with MFE <= -20",
    "Collapsed high-confidence circRNA-miRNA duplexes",
    "Validated high-confidence duplexes",
    "Selected circRNA-associated miRNAs",
    "Existing lncRNA triplet rows sharing selected miRNAs",
    "Strict cross-ncRNA convergent rows",
    "Raw miRWalk target rows for selected miRNAs",
    "Candidate miRNA-mRNA pairs",
    "All putative circRNA-miRNA-mRNA rows",
    "Strict DE-coherent circRNA triplets",
    "Strict triplets shared with lncRNA analysis",
    "New strict circRNA triplets",
    "Expressed non-DE predicted targets"
  ),
  n = c(
    nrow(intarna_records),
    nrow(hc_records),
    nrow(hc_duplexes),
    sum(hc_duplexes$duplex_passes_validation, na.rm = TRUE),
    nrow(selected_miRNAs),
    nrow(lnc_selected),
    sum(cross_ncRNA_modules$strict_circ_triplet_from_lnc, na.rm = TRUE),
    nrow(raw_mirwalk_targets %>% semi_join(selected_miRNAs, by = c("timepoint", "miRNA_key"))),
    nrow(candidate_miRNA_mRNA),
    nrow(all_putative),
    nrow(strict_triplets),
    nrow(shared_strict_triplets),
    nrow(new_strict_triplets),
    nrow(expressed_nonDE_targets)
  )
)
write_csv(summary_table, file.path(out_dir, "22_RUN_summary_counts.csv"))

sink(file.path(out_dir, "README_FINAL_analysis_summary.txt"))
cat("Final circRNA-lncRNA-miRNA-mRNA cross-reference\n")
cat("================================================\n\n")
cat("IntaRNA threshold: MFE <= ", energy_cutoff, " kcal/mol\n", sep = "")
cat("miRNA DE threshold: padj < ", alpha, " and |log2FC| >= ", mirna_lfc_cutoff, "\n", sep = "")
cat("mRNA DE threshold: padj < ", alpha, " and |log2FC| >= ", mrna_lfc_cutoff, "\n\n", sep = "")
print(summary_table, n = Inf)
cat("\nInterpretation classes\n")
cat("- Shared_with_existing_lncRNA_triplet: the same miRNA-mRNA axis was already present in a coherent lncRNA triplet.\n")
cat("- New_strict_circRNA_triplet_from_miRWalk: raw miRWalk target support plus a directionally coherent DE mRNA, absent from the existing lncRNA triplet universe.\n")
cat("- Expressed_nonDE: predicted miRNA target detected in the full mRNA results but not DE by the current cutoff; supplementary context only.\n")
cat("- Not_mapped_to_expression: target prediction could not be matched to the expression table.\n\n")
cat("Important limitations\n")
cat("- Repeated all-record IntaRNA hits are collapsed at circRNA genomic-locus and miRNA level.\n")
cat("- When several splice variants share the same locus, the exact variant carrying the predicted binding site remains unresolved unless variant-specific FASTA headers were used.\n")
cat("- These are computationally reconstructed candidate modules and do not demonstrate direct miRNA sponging.\n")
cat("- A candidate should be considered for the main figure only after biological prioritization; the script flags prior figure/priority overlap but does not make that editorial decision automatically.\n\n")
cat("Primary files to inspect\n")
cat("1. 03_HC_circRNA_miRNA_duplexes_locus_collapsed.csv\n")
cat("2. 08_CROSS_ncRNA_convergent_modules.csv\n")
cat("3. 12_STRICT_DE_coherent_circRNA_triplets.csv\n")
cat("4. 14_NEW_strict_circRNA_triplets_not_in_lncRNA_analysis.csv\n")
cat("5. 15_EXPRESSED_nonDE_circRNA_miRNA_targets.csv\n")
cat("6. 16_AUDIT_circRNA_target_filtering.csv\n")
cat("7. 17_CYTOSCAPE_cross_ncRNA_nodes.csv and 18_CYTOSCAPE_cross_ncRNA_edges.csv\n")
sink()

message("Done. Results written to: ", normalizePath(out_dir))
message("Review first: 03_HC_circRNA_miRNA_duplexes_locus_collapsed.csv")
message("Then: 08_CROSS_ncRNA_convergent_modules.csv and 12_STRICT_DE_coherent_circRNA_triplets.csv")

