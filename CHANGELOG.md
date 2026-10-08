# Changelog

## v2 — rebuilt from the raw original scripts (7 Oct 2026)

This is a rebuild, not a patch. The previous tree (`main` @ `8d50e82`, README
dated 14 Jul 2026) was derived from the 7 Jul starter package alone. Three
things arrived after it and were never integrated, and two of them changed the
shape of the repository.

### Why rebuild

| Arrived | Package | Consequence |
|---|---|---|
| 31 Jul 2026 | `NE_NRVM_circRNA_repo_addendum_v1` | circRNA reintroduced as an independent module. The old README stated circRNA was archived and "not included in the final 7-figure manuscript". |
| 18 Aug 2026 | `NE_NRVM_FINAL_repository_update_v3` | 16 *final* figure scripts. The old `scripts/R/03-14` were ported from the scripts those superseded. |
| — | `Babi analysis` | The chain that generates `hisat2_6_24_meta_*`, listed as a pending item in the old README and absent from the pipeline. |

The manuscript is also now 6 main + 5 supplementary figures, not the 7 the old
figure map described.

### What changed

**Structure**

- R stages renumbered into true dependency order, 01–26 (was 01–14).
  `scripts/checks/check_io_graph.R` now enforces it.
- Four new bash stages, 09b–09e, carrying the transcript-metadata chain.
- Four new bash stages, 13–16, for the circRNA module, plus
  `envs/environment_circrna.yml`.
- `data/frozen/` for inputs that cannot be regenerated.
- `provenance/` for upstream scripts that are evidence, not code to run.
- `archive/superseded_duplicates/` dropped: the superseded scripts were already
  in `archive/raw_scripts_original/`, and duplicating them doubled the bytes for
  nothing. `archive/ARCHIVE_NOTES.md` says where each one went instead.

**Correctness**

- biomaRt pinned to `https://sep2025.archive.ensembl.org` (release 115) through
  a single `connect_ensembl()`. Previously four unpinned call sites:
  `useEnsembl(version = 114)`, `useEnsembl(mirror = "uswest")` and two bare
  `useMart("ensembl")`.
- `Rattus_norvegicus.mRatBN7.2.115.gtf` corrected to `GRCr8.115.gtf`. The former
  does not exist — release 115 has no mRatBN7.2 build.
- `12_STRICT_DE_coherent_circRNA_triplets.csv` is no longer an orphan input: it
  is output `12_` of the cross-reference stage, which now runs before the two
  stages that read it.
- Browser-download duplicates (`name(1).csv`) removed from seven scripts, where
  they were listed as *first-preference* candidates and would silently select a
  stale copy.
- `getBM(values = ...)` without `filters = ...` now raises instead of silently
  returning the wrong table.
- Stage 09e checks for GNU awk. Under mawk the three-argument `match()` never
  populates its capture group, every line passes the filter, and the output is
  an unfiltered GTF that looks correct.
- Every stage writes to its own `results/R<NN>_*/`. Previously the raw `12_*`
  script and all three `Figure8_*` variants wrote identical filenames into one
  directory.
- Hard-coded paths removed from 17 scripts — five different machine roots
  between them. `check_hardcoded_paths.sh` enforces it.
- The curated KEGG allowlist moved out of an inline `tribble` into
  `metadata/curated_kegg_allowlist.tsv`, where both stages that need it read
  the same copy.

**Two reconstructions**

- `09b_feelnc_split_exon_gtfs.sh` — the step producing the `exons/` layout was
  in neither the Babi folder nor the bash notebook. Identified from the data as
  the FEELnc_codpot split and verified against the historical transcript counts
  (46,030 / 42,127 / 3,822, with 87 classifier orphans and 81 noORF).
- `04_circRNA_merge_export.R` — the step producing
  `raw_count_matrix_all_*`, the `0{6,7}_*_DE_complete_table.txt` pair and the
  four IntaRNA target FASTAs. Rebuilt from `03_tool_comparation.Rmd`, including
  its asymmetric `+1` start correction for CIRCexplorer only.

Both are labelled RECONSTRUCTED in their headers, in their output, and in
`docs/OPEN_ITEMS.md`.

**Two new panels**

Neither had a standalone script in any package; both were flagged `OPEN` in the
original audit.

- `13_figure3C_circRNA_heatmap.R` (FIG3C_CIRC)
- `17_volcanoes_mRNA_miRNA.R` (S3_MRNA_MIRNA)

**Findings recorded**

- Sample size is **n = 3**, not 4. Replicate 3 is excluded from the
  mRNA/lncRNA DESeq2, the miRNA matrices and the circRNA volcanoes alike. The
  original audit framed the circRNA 3v3 as an inconsistency; it is the 4v4
  historical notebooks that are the outlier.
- `sep2025.archive.ensembl.org` **with a trailing slash** breaks biomaRt with
  `Error: Can't retrieve empty body`. Normalised in `lib/ensembl.R`.
- Ensembl 114 and 115 are both GRCr8; 113 and earlier are mRatBN7.2. Most of
  the apparent reference conflict in the historical material is therefore not a
  conflict at all.
- Every row of the shipped FEELnc classifier output has `isBest == 1`, so that
  filter is a no-op on this data.

**Carried forward unchanged**

- `metadata/fastq_manifest.tsv` with the real server paths, sizes and MD5s for
  all 32 long-RNA files, plus `fastq_md5_checksums.txt` and
  `fastq_file_sizes.tsv`. Those were real work and are kept as-is; the 16
  small-RNA rows were merged in and remain `TO_FILL`.
- `metadata/reference_manifest.tsv` and `scripts/bash/00_download_references.sh`,
  which already had the reference lock right.
