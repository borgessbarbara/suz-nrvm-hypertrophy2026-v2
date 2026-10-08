# Conformance with the v3 handoff instructions

The `NE_NRVM_FINAL_repository_update_v3` package (18 Aug 2026) came with five
instruction documents. This file records, item by item, what each asked for and
where it landed — so the handoff can be checked without reading the code.

Source documents, preserved verbatim under
`archive/final_package_v3_original/`:

| Document | Role |
|---|---|
| `README_FIRST.md` | Folder guide and the final/provenance split |
| `05_METHODS_AND_AUDIT/circRNA_provenance_audit_READ_ME.md` | 7 discrepancies that "must remain visible" |
| `05_METHODS_AND_AUDIT/circRNA_methods_repository_summary.md` | circRNA methodology as supplied |
| `06_NOTES_FOR_COLLABORATOR/handoff_summary.txt` | Top 5 validation priorities |
| `06_NOTES_FOR_COLLABORATOR/repository_merge_plan.md` | Recommended consolidation plan |

---

## README_FIRST.md

| Instruction | Status |
|---|---|
| Split final manuscript-aligned scripts from raw circRNA provenance | **Done.** `scripts/` holds the runnable stages; `archive/` holds the circRNA shell scripts and notebooks verbatim, never executed. |
| "Raw circRNA provenance scripts … should be consolidated and parameterized rather than run unchanged" | **Done.** Bash 13-16 and R04/R04b consolidate them; `archive/` is excluded from every runner. |
| Read the provenance audit before changing references, BSJ filters or DE logic | **Done.** All seven items are in `OPEN_ITEMS.md` §A, and each parameter is a CLI flag rather than an edit. |
| `24_FINAL_..._cross_reference.R` already exists as `scripts/R/15_circRNA_cross_reference.R`; do not duplicate | **Superseded.** The repository was rebuilt from zero at the user's instruction, so the old R15 no longer exists. That script is now `scripts/R/22_circRNA_cross_reference.R`, ported once. No duplicate. |
| Figures need not be pixel-identical (Affinity, Cytoscape) | **Done.** `MANUAL_ARTIFACTS.md`. |

## circRNA_provenance_audit_READ_ME.md — the seven discrepancies

All seven were independently re-verified against the code during the rebuild.
All seven hold.

| # | Discrepancy | Where it is now |
|---|---|---|
| 1 | Reference files: GRCr8.114 / rn8.refflat / rn8.yml vs the stated lock. "Do not silently replace the historical paths in the provenance files." | `REFERENCE_AUDIT.md`; `OPEN_ITEMS.md` C1–C3. Provenance paths untouched; the runnable stages take the reference as an argument. Also escalated: the lock as written in the Methods names a GTF that does not exist. |
| 2 | BSJ filter: prose ≥2 vs code `> 2` | `OPEN_ITEMS.md` A1, now **BLOCKING** — manuscript v14 commits to ≥2. `R04b --bsj-min`, logged on every run. |
| 3 | Historical DE: ComBat_seq + \|LFC\|>1 vs final nominal p + \|LFC\|≥0.5 | `OPEN_ITEMS.md` A2. `R04b --use-combat`, off by default. |
| 4 | Replicates: notebooks 4+4 vs final 3+3 | `OPEN_ITEMS.md` A3. **Largely resolved**: every transcriptomic figure legend in v14 states N = 3, and the mRNA/lncRNA and miRNA analyses are 3v3 too — the notebooks are the outlier. `R04b --replicates`, default `1,2,4`. A4 (which physical library "rep3" was) is still open. |
| 5 | CIRI-AS covers SQ13–SQ16 only | `OPEN_ITEMS.md` A5. Bash 15 defaults to those four and refuses `--all` without `--yes-widen`. |
| 6 | CIRCexplorer3 in prose vs CIRCexplorer2 + CLEAR in code | `OPEN_ITEMS.md` A6, `METHODS_NOTES.md` §2. Actual commands preserved in bash 13; both pinned in `environment_circrna.yml`. |
| 7 | Final table lineage not explicit | `OPEN_ITEMS.md` A7. R04 + R04b reconstruct it, labelled RECONSTRUCTED, asserting the recorded 46,302 / 3,419 / 23,739 counts. |
| — | Manual figure provenance (Affinity, Cytoscape) | `MANUAL_ARTIFACTS.md`. |

## circRNA_methods_repository_summary.md

Every tool and version it records is in `metadata/software_tool_inventory.tsv`
with its evidence: HISAT2 2.2.1, STAR 2.7.11b, CIRIquant 1.1, BWA-MEM `-T 19`,
CIRI2 2.0.6, CIRI-AS 1.2, CIRCexplorer2 + CLEAR. Its own caveat — that the BSJ
threshold "requires final audit" — is open item A1. Its statement that circRNA
DE is exploratory at nominal p < 0.05 and \|log2FC\| ≥ 0.5 matches R20 and
manuscript v14.

## handoff_summary.txt — the five priorities

| # | Priority | Status |
|---|---|---|
| 1 | GRCr8/release-115 vs historical rn8/GRCr8.114 | **Resolved as far as the evidence allows.** 114 and 115 are both GRCr8, so this is an annotation-release gap, not an assembly conflict. `REFERENCE_AUDIT.md`. |
| 2 | BSJ threshold | Open, BLOCKING. A1. |
| 3 | circRNA DE provenance | Open. A2, `--use-combat`. |
| 4 | 3 vs 4 replicates | Largely resolved; the manuscript says N = 3. A3. |
| 5 | Final circRNA count/sequence table lineage | Reconstructed in R04 + R04b. A7. |

## repository_merge_plan.md

### A. Compare/update rather than duplicate

Followed in spirit. The rebuild started from zero, so there was nothing to
compare against — but the sound bash stages 00–12 were carried over from the
previous repository rather than rewritten, and each ported R stage exists
exactly once, diffable against its archived original.

### B. The seven recommended consolidated stages

| Seba's recommendation | This repository |
|---|---|
| 1. `circRNA_01_circexplorer.sh` | `scripts/bash/13_circrna_circexplorer.sh` |
| 2. `circRNA_02_ciriquant.sh` | `scripts/bash/14_circrna_ciriquant.sh` |
| 3. `circRNA_03_ciri_as.sh` | `scripts/bash/15_circrna_ciri_as.sh` |
| 4. `circRNA_04_annotation_merge.R` | `scripts/R/04_circRNA_annotation_merge.R` |
| 5. `circRNA_05_final_DE_temporal.R` | `scripts/R/04b_circRNA_de_temporal.R` (+ temporal classification in `13_figure3C_circRNA_heatmap.R`) |
| 6. `circRNA_06_intarna.sh` | `scripts/bash/16_circrna_intarna.sh` |
| 7. Existing `scripts/R/15_circRNA_cross_reference.R` | `scripts/R/22_circRNA_cross_reference.R` |

**One deliberate departure: numbering.** The stages are numbered into the single
global pipeline sequence rather than a separate `circRNA_0N_` series, so that
`check_io_graph.R` can enforce one dependency order across all 26 R stages and
21 bash stages. The content of each stage follows the plan; only the filename
prefix differs. That is what this table is for.

**One correction made in response to this check.** Stages 4 and 5 were
initially collapsed into a single `04_circRNA_merge_export.R`. They are now
split as the plan specifies, and the split is load-bearing rather than
cosmetic: R04 is deterministic (parse, reconcile coordinates, merge loci, emit
an *unfiltered* BSJ matrix), while every contested parameter — the BSJ rule,
the replicate set, ComBat_seq — lives in R04b. Because A1, A2 and A3 are all
still open, DE has to be re-runnable under either BSJ rule without re-parsing
every caller output. Seba's split gives exactly that.

### C. Do not require pixel-identical figures

Done. `MANUAL_ARTIFACTS.md` states what the pipeline guarantees (numerical
data, selected features, source panels, network nodes/edges/attributes) and
what it does not (page composition, Cytoscape coordinates), and asks for the
`.afdesign` and `.cys` files.

### D. Final validation before freezing

| Item | Status |
|---|---|
| Reference version/checksum | Version resolved (`REFERENCE_AUDIT.md`). **Checksums still `TO_FILL`** in `reference_manifest.tsv`. |
| BSJ threshold | Open, BLOCKING. A1. |
| circRNA sample inclusion | Largely resolved (N = 3). A3; A4 open. |
| Whether ComBat_seq belongs in the final branch | Open. A2. |
| Final counts in Figures 2–3 **and Table 4** | Figure 2B **verified exactly**. Figure 3A/3C and S3C reconcile; **Figure 3B's lncRNA panel does not** — see `MANUSCRIPT_ALIGNMENT.md` §4.1. "Table 4" does not exist in v14; the workbook is S1–S8. Assumed to mean Supplementary Table S4 (lncRNA summary), which matches Figure 2C. **Confirm.** |
| IntaRNA output filenames and MFE ≤ −20 filtering | Verified. The eight output names match the `Interact_*.csv` glob in R21; −20 kcal/mol is applied in R21/R24/R25. One caveat: the Methods say "no additional interaction tiers were used" while the ceRNA script defines three. `MANUSCRIPT_ALIGNMENT.md` §3.6. |

---

## Not followed, and why

| Instruction | Reason |
|---|---|
| Keep the `circRNA_0N_` stage naming | Numbered into the global sequence instead, so the I/O dependency check covers every stage uniformly. Crosswalk above. |
| "Do not duplicate `scripts/R/15_circRNA_cross_reference.R`" | The repository was rebuilt from zero at the user's explicit instruction, so there was no existing stage to preserve. It exists once, as R22. |
