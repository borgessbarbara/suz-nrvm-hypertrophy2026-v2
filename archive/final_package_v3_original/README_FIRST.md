# NE-NRVM FINAL repository update v3

This is the final handoff package prepared after recovery of the missing circRNA scripts and exact IntaRNA commands.

## Important
This package is intentionally split between:
- **final manuscript-aligned scripts** to compare/update in the repository, and
- **raw circRNA provenance scripts** that should be consolidated and parameterized rather than run unchanged.

## Folder guide

### `01_FINAL_FIGURE_AND_ANALYSIS_SCRIPTS/`
Final manuscript-aligned R scripts previously reviewed against the 6 main figures and 5 supplementary figures.

### `02_CIRCRNA_UPSTREAM_RAW_PROVENANCE/`
All circRNA scripts received from the circRNA co-author, preserved unchanged.

### `03_INTARNA_COMMANDS/`
Exact eight IntaRNA commands reported as used, plus a state/direction manifest.

### `04_METADATA/`
Panel-to-script map, circRNA script audit, software/reference audit, manual-editing provenance, open issues, and status of the already-parameterized R15 stage.

### `05_METHODS_AND_AUDIT/`
CircRNA methods summary and a critical provenance audit. Read this before changing references, BSJ filters, or DE logic.

### `06_NOTES_FOR_COLLABORATOR/`
Recommended merge/consolidation plan.

## Already in repository
`24_FINAL_circRNA_lncRNA_miRNA_mRNA_cross_reference.R` was already parameterized as `scripts/R/15_circRNA_cross_reference.R`; it is intentionally not duplicated here.

## Manual figure editing
Some final figures were assembled/edited in Affinity and Cytoscape layouts were manually adjusted. Script outputs therefore need to reproduce the analytical content, but not necessarily the pixel-identical final composition.
