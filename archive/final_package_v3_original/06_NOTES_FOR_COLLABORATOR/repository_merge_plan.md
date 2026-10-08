# Repository merge plan

## A. Compare/update final manuscript scripts
Use `01_FINAL_FIGURE_AND_ANALYSIS_SCRIPTS/` as the manuscript-aligned source set. If an equivalent parameterized script already exists in the private repository, compare analytical logic and update the existing stage rather than adding a duplicate.

## B. Build the circRNA upstream module from raw provenance
Use `02_CIRCRNA_UPSTREAM_RAW_PROVENANCE/` as historical evidence, not as production-ready code.

Recommended consolidated stages:
1. `circRNA_01_circexplorer.sh` — one sample-sheet-driven CIRCexplorer workflow replacing the three batch-specific shell scripts.
2. `circRNA_02_ciriquant.sh` — one sample-sheet-driven CIRIquant workflow replacing the three batch-specific shell scripts.
3. `circRNA_03_ciri_as.sh` — parameterized BWA-MEM/CIRI2/CIRI-AS reconstruction.
4. `circRNA_04_annotation_merge.R` — fixed-resource annotation and CIRCexplorer/CIRIquant locus merge based on the Rmd logic.
5. `circRNA_05_final_DE_temporal.R` — final manuscript circRNA DE/temporal classification using the verified sample set and verified BSJ rule.
6. `circRNA_06_intarna.sh` — parameterized version of the exact eight IntaRNA commands.
7. Existing repo stage `scripts/R/15_circRNA_cross_reference.R`.

## C. Do not require pixel-identical figure reproduction
Affinity was used for final figure assembly and Cytoscape layouts were manually adjusted. The repo should reproduce:
- numerical data,
- selected features,
- plot source panels,
- network nodes/edges/attributes.

For exact layout provenance, add Affinity source files and the final Cytoscape `.cys`/node-coordinate files if available.

## D. Final validation
Before freezing the release, validate:
- reference version/checksum,
- BSJ threshold,
- circRNA sample inclusion,
- whether ComBat_seq belongs in the final circRNA DE branch,
- final counts in Figures 2-3 and Table 4,
- IntaRNA output filenames and MFE <= -20 filtering.
