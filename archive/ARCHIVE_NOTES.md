# Archive notes

`archive/` holds the material as received. Nothing here is executed, nothing
here has been path-corrected, and nothing here should be edited. It exists so
that every ported stage can be diffed against its source.

**One deliberate change: line endings.** Every file as received had CRLF
terminators — they came off a Windows machine. `.gitattributes` normalises the
whole repository to LF. The content is otherwise untouched, and the
normalisation is what makes the audit diffs readable: against a CRLF original,
an LF port shows *every* line as changed and the five real differences are
invisible. The original CRLF copies remain in
`SUZ2026/NE_NRVM_all_scripts_raw/` and
`SUZ2026/NE_NRVM_FINAL_repository_update_v3/` if a byte-exact comparison is
ever needed; `archive/final_package_v3_original/04_METADATA/package_md5.tsv`
holds that package's own MD5 manifest.

## What is here

| Directory | Contents | Received |
|---|---|---|
| `raw_scripts_original/` | `NE_NRVM_all_scripts_raw` — 23 R scripts, the bash command notebook, the IntaRNA history, the miRge commands | 2 Jul 2026 |
| `final_package_v3_original/` | `NE_NRVM_FINAL_repository_update_v3` — 16 final figure scripts, the circRNA upstream provenance, the eight IntaRNA commands, the audit metadata | 18 Aug 2026 |

Two further packages were received and are represented elsewhere rather than
duplicated here:

- `NE_NRVM_coauthor_reproducibility_starter_v1` (7 Jul 2026) — its sample
  sheets and manifests are the basis of `metadata/`.
- `NE_NRVM_circRNA_repo_addendum_v1` (31 Jul 2026) — its cross-reference script
  is ported as `scripts/R/22_circRNA_cross_reference.R`; its manifests fed
  `metadata/circRNA_input_manifest.tsv`.

## Where each archived script went

### Ported into the pipeline

| Archived original | Stage |
|---|---|
| `Integrated_Bulk_v3_6_24h.R` | `scripts/R/02_deseq2_mRNA_lncRNA.R` |
| `miRNAs_6_24.R` | `scripts/R/03_mirna_deseq2.R` |
| `script_final_depurado_TF_KEGG.R` | `scripts/R/05_tf_kegg_integration.R` |
| `26_transcriptome_overview_with_circRNA.R` | `scripts/R/06_transcriptome_overview.R` |
| `15_lncRNA_characterization_from_GTF.R` | `scripts/R/07_lncRNA_characterization_gtf.R` |
| `14_FEELnc_supplement_panels.R` | `scripts/R/08_feelnc_panels.R` |
| `35_FINAL_missing_figures_miRNA_heatmap_DE_GBA.R` | `scripts/R/09_de_counts_heatmaps_gba.R` |
| `36b_Figure3_panelB_DOWN_bottom_UP_top.R` | `scripts/R/10_figure3B_temporal_programs.R` |
| `05_Figure3_mRNA_programs.R` | `scripts/R/11_figure3C_mRNA_heatmap.R` |
| `37c_lncRNA_temporal_heatmaps_annotated_novel_MSTRG_FIXED.R` | `scripts/R/12_figure3C_lncRNA_heatmap.R` |
| `script_miRNA_mRNA_6_24.R` | `scripts/R/14_mirna_mRNA_integration.R` |
| `06_Figure4_refine_heatmap_hubs.R` | `scripts/R/15_lncRNA_hubs.R` |
| `07_Figure4_clean_volcanoes.R` | `scripts/R/16_volcanoes_lncRNA.R` |
| `39e_lncRNA_directional_overlap_barplot_UP_top_DOWN_bottom_FINAL.R` | `scripts/R/18_directional_overlap.R` |
| `09_extract_lncRNA_candidates_for_IntaRNA.R` | `scripts/R/19_extract_lncRNA_for_intarna.R` |
| `25_recreate_circRNA_volcanoes_with_all_NS.R` | `scripts/R/20_circRNA_volcanoes.R` |
| `12_ceRNA_coherent_pairs_and_triplets_ALL_RECORDS_FINAL.R` | `scripts/R/21_ceRNA_coherent_pairs_triplets.R` |
| `24_FINAL_circRNA_lncRNA_miRNA_mRNA_cross_reference.R` (addendum) | `scripts/R/22_circRNA_cross_reference.R` |
| `31_phase_aware_ceRNA_coherence_plots.R` | `scripts/R/23_ceRNA_phase_aware_plots.R` |
| `32f_export_ceRNA_cytoscape_TF11_curatedKEGG_MFE20.R` | `scripts/R/24_ceRNA_cytoscape_export.R` |
| `33_validate_fix_cytoscape_import_keys.R` | `scripts/R/25_cytoscape_key_validation.R` |
| `MASTER_results_report_builder_v3_8_1.R` | `scripts/R/26_results_report.R` |

Those ports were produced mechanically by `scripts/dev/port_from_archive.sh`,
which applies five classes of change and nothing else. To audit any of them:

```bash
diff archive/raw_scripts_original/R_scripts/miRNAs_6_24.R scripts/R/03_mirna_deseq2.R
```

Expect to see only: the new preamble, commented-out `base_dir`/`out_dir` lines,
`base = base_dir` substitutions, removed `name(1).csv` duplicates, and the
removed `/mnt/data` search path. Three files needed small hand edits afterwards
because a duplicate filename was the last element of a multi-line vector, which
a line-based transform cannot see; those edits are marked `## [ported]` in place.

### Superseded, not ported

These stay here and only here. They are not duplicated into a
`superseded/` folder — that would just double the bytes for no gain.

| Script | Why |
|---|---|
| `Figure8_ceRNA_rebuild_with_QC.R`, `_v2.R`, `_v3.R` | Three generations of the same analysis, all superseded by `12_..._ALL_RECORDS_FINAL.R`. v1→v2 is 88 changed lines; v2→v3 is effectively a rewrite. All three write into the same output directory as `12_*`, so running them together would overwrite each other. |
| `12_ceRNA_coherent_pairs_and_triplets_cleanFinal.R` | Superseded by the ALL_RECORDS version (601 changed lines, adds locus collapsing and the AUDIT_* outputs). |
| `08_Figure7_ceRNA_IntaRNA_scaffold.R` | Prototype scaffold; superseded by the same. |
| `13_generate_paper_results_report_v2.R` | Superseded by `MASTER_results_report_builder_v3_8_1.R`. |
| `16_lncRNA_characterization_curatedCSV_only.R` | The CSV-driven alternative to the GTF-driven characterization. Dropped from the final handoff. Produces **different** inventory numbers from a different source, so it must not be mixed with R07. |
| `17_redraw_remaining_panels_Fig4_Fig5_Fig6.R` | Re-plotting pass whose figure numbering is historical. Its panels are covered by R10/R12/R15/R18. |
| `18_curated_Fig4A_Fig6C_from_shortlist.R` | Curated shortlist plots driven from `C:/Users/sebau/Downloads`. The curation now lives in `metadata/curated_kegg_allowlist.tsv`. |
| `Comparación_Cinética_6_24.R` | Earlier multifactorial DESeq2 pass. Contains an unpinned biomaRt call. Superseded by `Integrated_Bulk_v3_6_24h.R`. |
| `Versión_Nueva_Nov.R` | Another earlier pass. Unpinned biomaRt with `mirror = "uswest"` — a mirror is not a version pin. |
| `miRNA_DEG_ceRNAs.R` | Historical miRNA/Cytoscape script. Imports `biomaRt` but never calls it. Reads an Excel file with a doubled `.xls.xlsx` extension. Superseded by R03 and R24. |

### Not carried over at all

The `Babi analysis` Quarto notebooks are **not** in this repository, by
decision taken during the rebuild: the metadata-generation logic was ported to
a clean Rscript (`scripts/R/01_transcriptome_metadata_annotation.R`) instead.

The originals remain at:

```
SUZ2026/Babi analysis/2_nrvm_transcriptome_metadata_handling.qmd     (v2, 671 lines - the one that ran)
SUZ2026/Babi analysis/2_nrvm_transcriptome_metadata_handling.html    (rendered)
SUZ2026/Babi analysis/nrvm_transcriptome_metadata_handling.qmd       (v1, 999 lines - superseded, outputs absent)
SUZ2026/Babi analysis/nrvm_transcriptome_visualization.qmd           (40+ panels, no ggsave calls)
SUZ2026/Babi analysis/nrvm_transcriptome_visualization.html          (rendered)
```

The three shell scripts from that folder **are** in the repository, as bash
stages 09c, 09d and 09e, because R01 is unreproducible without them. Their
verbatim originals are in `provenance/babi_preprocessing/`.

The STAR half of that analysis is not carried over either: S. Urquiza confirmed
only HISAT2 was used. That removes eight of the twelve `*_meta_*` /
`*_feelnc_*` outputs and roughly 426 MB, and moots the unexplained STAR novel-lncRNA
deficit (48 and 63 novel transcripts for STAR, against 2,597 and 3,786 for HISAT2).

## Relationship to `provenance/`

`archive/` is "everything as received". `provenance/` is the subset that
documents how upstream data was produced and that is referred to by the active
stages — the circRNA shell scripts and notebooks, the IntaRNA command history,
the bash command notebook, the miRge commands, and the three Babi shell
scripts. Some files appear in both; that is intentional and costs a few hundred
kilobytes.
