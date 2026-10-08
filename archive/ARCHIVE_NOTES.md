# Archive notes

`archive/` holds the material as received. Nothing here is executed, nothing
here has been path-corrected, and nothing here should be edited. It exists so
that every ported stage can be diffed against its source, and so that every
claim in the documentation can be checked against the original.

Everything received lives here **once**. Nothing else in the repository
duplicates it, and the runnable pipeline under `scripts/` never reads from it —
script headers only cite it.

**One deliberate change: line endings.** Every file as received had CRLF
terminators. `.gitattributes` normalises the whole repository to LF. The
content is otherwise untouched, and the normalisation is what makes the audit
diffs readable: against a CRLF original, an LF port shows *every* line as
changed and the real differences are invisible. The original CRLF copies remain
outside the repository if a byte-exact comparison is ever needed, and
`final_package_v3_original/04_METADATA/package_md5.tsv` holds that package's
own MD5 manifest. `metadata/frozen_inputs_checksums.tsv` records the current
size and MD5 of every file here.

## What is here

| Directory | Contents | Received |
|---|---|---|
| `raw_scripts_original/` | `NE_NRVM_all_scripts_raw` — 23 R scripts, the bash command notebook, the IntaRNA history, the miRge commands | 2 Jul 2026 |
| `starter_v1_original/` | `NE_NRVM_coauthor_reproducibility_starter_v1` — sample sheets, reference and software manifests, the script inventory, the draft message | 7 Jul 2026 |
| `circrna_addendum_v1_original/` | `NE_NRVM_circRNA_repo_addendum_v1` — the cross-reference script and the circRNA manifests | 31 Jul 2026 |
| `final_package_v3_original/` | `NE_NRVM_FINAL_repository_update_v3` — 16 final figure scripts, the circRNA upstream scripts and notebooks, the eight IntaRNA commands, the audit metadata | 18 Aug 2026 |
| `transcriptome_metadata_original/` | The three shell scripts of the transcriptome-metadata chain that R01 depends on | B. Borges |

The v3 package's own folder names are preserved verbatim, including
`02_CIRCRNA_UPSTREAM_RAW_PROVENANCE`. When a script header or a document cites

```
archive/final_package_v3_original/02_CIRCRNA_UPSTREAM_RAW_PROVENANCE/03_tool_comparation.Rmd:16
```

that is the only copy of that file in the repository.

## Where each archived script went

### Ported into the pipeline

| Archived original | Stage |
|---|---|
| `extract_transcript_ids.sh` | `scripts/bash/09c_extract_transcript_ids.sh` |
| `subset_merged_gtfs.sh` | `scripts/bash/09d_subset_merged_gtfs.sh` |
| `filter_out_lncRNAs.sh` | `scripts/bash/09e_filter_out_lncRNAs.sh` |
| `CIRCexp3_run_all_2023_SebaNECardio{,_2,_3}.sh` | `scripts/bash/13_circrna_circexplorer.sh` |
| `CIRIquant_all_run_2023_Seba_NE_Cardio_{1,2,3}.sh` | `scripts/bash/14_circrna_ciriquant.sh` |
| `CIRI_2_run_all_seba.sh` | `scripts/bash/15_circrna_ciri_as.sh` |
| `IntaRNA_commands_as_used.sh` | `scripts/bash/12_intarna_run.sh` (4 commands) and `16_circrna_intarna.sh` (4 commands) |
| `Integrated_Bulk_v3_6_24h.R` | `scripts/R/02_deseq2_mRNA_lncRNA.R` |
| `miRNAs_6_24.R` | `scripts/R/03_mirna_deseq2.R` |
| `03_tool_comparation.Rmd` | `scripts/R/04_circRNA_annotation_merge.R` + `04b_circRNA_de_temporal.R` (reconstruction) |
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
| `24_FINAL_circRNA_lncRNA_miRNA_mRNA_cross_reference.R` | `scripts/R/22_circRNA_cross_reference.R` |
| `31_phase_aware_ceRNA_coherence_plots.R` | `scripts/R/23_ceRNA_phase_aware_plots.R` |
| `32f_export_ceRNA_cytoscape_TF11_curatedKEGG_MFE20.R` | `scripts/R/24_ceRNA_cytoscape_export.R` |
| `33_validate_fix_cytoscape_import_keys.R` | `scripts/R/25_cytoscape_key_validation.R` |
| `MASTER_results_report_builder_v3_8_1.R` | `scripts/R/26_results_report.R` |

To audit a port, diff it against its original:

```bash
diff archive/raw_scripts_original/R_scripts/miRNAs_6_24.R \
     scripts/R/03_mirna_deseq2.R
```

Every class of difference you will see is listed in
[`../docs/CHANGES_FROM_ORIGINALS.md`](../docs/CHANGES_FROM_ORIGINALS.md), which
is ordered the way the pipeline runs.

### Written new, with no original

| Stage | Why |
|---|---|
| `scripts/bash/09b_feelnc_split_exon_gtfs.sh` | RECONSTRUCTED — the command producing the `exons/` layout was never recorded |
| `scripts/R/04_circRNA_annotation_merge.R`, `04b_circRNA_de_temporal.R` | RECONSTRUCTED from `03_tool_comparation.Rmd` |
| `scripts/bash/16_circrna_intarna.sh` | RECONSTRUCTED — the four commands exist, the runner did not |
| `scripts/R/13_figure3C_circRNA_heatmap.R` | NEW — no standalone script for this panel in any package |
| `scripts/R/17_volcanoes_mRNA_miRNA.R` | NEW — the plotting code was embedded in two large scripts and never isolated |
| `scripts/R/01_transcriptome_metadata_annotation.R` | ported from a Quarto notebook that stays outside the repository (see below) |

### Superseded, not ported

These stay here and only here. They are not duplicated into a `superseded/`
folder, which would double the bytes for no gain.

| Script | Why |
|---|---|
| `Figure8_ceRNA_rebuild_with_QC.R`, `_v2.R`, `_v3.R` | Three generations of the same analysis, all superseded by `12_..._ALL_RECORDS_FINAL.R`. v1→v2 is 88 changed lines; v2→v3 is effectively a rewrite. All three write into the same output directory as `12_*`, so running them together would overwrite each other. |
| `12_ceRNA_coherent_pairs_and_triplets_cleanFinal.R` | Superseded by the ALL_RECORDS version: 601 changed lines, adding locus collapsing and the `AUDIT_*` outputs. |
| `08_Figure7_ceRNA_IntaRNA_scaffold.R` | Prototype scaffold; superseded by the same. |
| `13_generate_paper_results_report_v2.R` | Superseded by `MASTER_results_report_builder_v3_8_1.R`. |
| `16_lncRNA_characterization_curatedCSV_only.R` | The CSV-driven alternative to the GTF-driven characterization, dropped from the final handoff. Produces **different** inventory numbers from a different source, so it must not be mixed with R07. |
| `17_redraw_remaining_panels_Fig4_Fig5_Fig6.R` | Re-plotting pass whose figure numbering is historical. Its panels are covered by R10/R12/R15/R18. |
| `18_curated_Fig4A_Fig6C_from_shortlist.R` | Curated shortlist plots driven from a downloads folder. The curation now lives in `metadata/curated_kegg_allowlist.tsv`. |
| `Comparación_Cinética_6_24.R` | Earlier multifactorial DESeq2 pass, with an unpinned biomaRt call. Superseded by `Integrated_Bulk_v3_6_24h.R`. |
| `Versión_Nueva_Nov.R` | Another earlier pass. Unpinned biomaRt with `mirror = "uswest"` — a mirror is not a version pin. |
| `miRNA_DEG_ceRNAs.R` | Historical miRNA/Cytoscape script. Imports `biomaRt` but never calls it, and reads an Excel file with a doubled `.xls.xlsx` extension. Superseded by R03 and R24. |
| `00_annotation_biomart.Rmd` | Carries `#NO VOLVER A EJECUTAR` at line 12, and its `getBM` calls pass `values` with no `filters`. Not in the runnable path. |
| `vis_run_ciri2_seba.sh` | Visualization pass for CIRI2. Its sample list is misspelled and it deletes the PDFs it generates. Not ported. |

## Not carried over at all

**The Quarto notebooks of the transcriptome-metadata analysis.** By decision
taken during the rebuild, no `.qmd` enters this repository; the
metadata-generation logic was ported to a clean Rscript,
`scripts/R/01_transcriptome_metadata_annotation.R`. The originals remain
outside the repository, in the author's analysis folder:

```
2_nrvm_transcriptome_metadata_handling.qmd     v2, 671 lines — the one that ran
2_nrvm_transcriptome_metadata_handling.html    rendered
nrvm_transcriptome_metadata_handling.qmd       v1, 999 lines — superseded, outputs absent
nrvm_transcriptome_visualization.qmd           40+ panels, no ggsave calls
nrvm_transcriptome_visualization.html          rendered
```

The three **shell** scripts from that folder are in the repository, as bash
stages 09c-09e, because R01 is unreproducible without them. Their verbatim
originals are in `transcriptome_metadata_original/`.

**The STAR branch of that analysis.** S. Urquiza confirmed only HISAT2 was
used. That removes eight of the twelve `*_meta_*` / `*_feelnc_*` outputs and
roughly 426 MB, and moots the unexplained STAR novel-lncRNA deficit — 48 and 63
novel transcripts for STAR, against 2,597 and 3,786 for HISAT2.
