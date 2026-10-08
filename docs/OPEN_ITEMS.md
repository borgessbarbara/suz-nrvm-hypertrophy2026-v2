# Open items

Everything that must be confirmed before this repository is frozen, with an
owner and the evidence behind it. Nothing here is a guess: each item cites the
file and line it comes from.

Items are grouped by who can resolve them.

---

## A. For Sebastián and Juan Pablo — circRNA branch

These are the seven points from the original provenance audit
(`archive/final_package_v3_original/05_METHODS_AND_AUDIT/circRNA_provenance_audit_READ_ME.md`).
All seven were re-verified against the code during this rebuild and all seven
hold. Two of them have changed character; see the notes.

### A1. BSJ read-support threshold — prose says ≥2, code says ≥3

- Methods prose: "two or more BSJ reads".
- `provenance/circRNA_upstream/01_circexp3.Rmd:91` — `count_matrix[count_matrix$readNumber > 2,]`
- `provenance/circRNA_upstream/02_ciriquant.Rmd:85` — `df[df$readNumber > 2,]`

`> 2` means **at least 3**. These are not the same rule and they give different
inventories. `scripts/R/04_circRNA_merge_export.R --bsj-min` makes the choice
explicit and logs it; the default of 3 reproduces the code, not the prose.

**Needed:** which rule produced the deposited circRNA inventory, and a
correction to whichever of the two is wrong.

### A2. circRNA DE derivation — ComBat_seq and the fold-change threshold

- Historical notebooks: `ComBat_seq` then `pvalue < 0.05 & |log2FoldChange| > 1`
  (`01_circexp3.Rmd:169-170, 253-254`; `02_ciriquant.Rmd:158-159, 273-274`).
- Final figure script: nominal `p < 0.05` and `|log2FC| >= 0.5`, **no ComBat_seq**
  (`25_recreate_circRNA_volcanoes_with_all_NS.R:26-27`, ported as R20).

**Needed:** confirmation that the final tables came from the second route.
`--use-combat` in R04 exists so the historical route can be reproduced for
comparison; it is off by default.

### A3. Replicate inclusion — **this item has changed**

The original audit framed this as "the notebooks use 4+4, the final circRNA
script uses 3+3; reconcile". Re-checking the whole corpus during the rebuild
shows the final circRNA script is **consistent with everything else**, and the
circRNA notebooks are the outlier:

| Analysis | Scheme | Evidence |
|---|---|---|
| mRNA / lncRNA DESeq2 | **3 v 3**, rep 1/2/4 | `Integrated_Bulk_v3_6_24h.R:30-33` — `Ctrl6_rep1, Ctrl6_rep2, Ctrl6_rep4, ...` |
| miRNA DESeq2 | **3 v 3** | `miRNAs_6_24.R:21-22` — `miR.Counts_C6_NE6_3v3_new.csv` |
| circRNA, final volcano | **3 v 3**, `_3` dropped | `25_...R:284-291` |
| circRNA, historical notebooks | 4 v 4 + ComBat_seq over Batch1-4 | `01_circexp3.Rmd:122-124, 210-212` |

So **replicate 3 is excluded from every final analysis**, and the manuscript is
internally consistent at 3 v 3. Two consequences:

1. The sample sheets describe 16 long-RNA samples, n=4 per group. The analyses
   used 12. `metadata/longRNA_samples.tsv` now carries
   `include_final_analysis = no_excluded_from_final_DE` for the four replicate-3
   rows. **Confirm this is right**, and make sure the Methods text says n=3.
2. In the circRNA notebooks `batch` maps 1:1 to replicate index, so dropping
   replicate 3 removes **Batch3** entirely. If ComBat_seq is retained anywhere,
   that matters.

**Needed:** confirmation that n=3 is intended and correctly described, and the
mapping below.

### A4. The `Ctrl6_rep1` aliases are not mapped to sample IDs

`Integrated_Bulk_v3_6_24h.R` reads Salmon output from directories named
`Ctrl6_rep1_quant` … `NE24_rep4_quant`. **No file in any package maps those
aliases to `H1_ln`, `SQ13_ln`, …** So although we know replicate 3 was dropped,
we cannot say with certainty *which physical library* that was.
`metadata/longRNA_samples.tsv` has a `salmon_quant_alias` column, currently
`TO_CONFIRM` for every row.

**Needed:** the alias → sample_id mapping.

### A5. CIRI-AS covers SQ13–SQ16 only

`provenance/circRNA_upstream/CIRI_2_run_all_seba.sh:2` — `for m in {13..16}`.
There is no CIRI2/CIRI-AS loop for H1–H4 or SQ1–SQ8 anywhere.
`scripts/bash/15_circrna_ciri_as.sh` therefore defaults to those four samples
and refuses `--all` without `--yes-widen`.

**Needed:** either the missing batch commands, or confirmation that
reconstruction genuinely only ever covered SQ13–SQ16. This gates the IntaRNA
circRNA targets (see A7).

### A6. CIRCexplorer2 vs CIRCexplorer3

The methods say "CIRCexplorer3 v1.0.0". The scripts invoke `CIRCexplorer2 parse`
and `CIRCexplorer2 annotate` plus CLEAR's `run3.py` and `circ_quant.py`; the
string `circexplorer3` appears only as a directory name. The commands are
preserved as-is in `scripts/bash/13_circrna_circexplorer.sh`.

**Needed:** a decision on what the Methods section should say. See
`docs/METHODS_NOTES.md`.

### A7. Final circRNA table lineage — four orphan file families

Nothing in any package produces:

| File | Read by |
|---|---|
| `raw_count_matrix_all_{6,24}.txt` | R06, R20 |
| `0{6,7}_circRNA_sequences_de{6,24}_DE_complete_table.txt` | R13, R20, R21, R22 |
| `circRNA_{UP,DOWN}_{6,24}h_all_records.fa` | bash 16 |
| `12_STRICT_DE_coherent_circRNA_triplets.csv` | R23, R24 |

The last one turned out **not** to be an orphan: it is output `12_` of the
cross-reference stage. Reordering so R22 runs before R23/R24 fixed it.

The other three are rebuilt by `scripts/R/04_circRNA_merge_export.R`, marked
RECONSTRUCTED. The FASTA export additionally depends on CIRI-AS sequences,
which A5 puts in doubt.

**Needed:** confirmation that the reconstruction matches what was actually run.
It asserts the counts recorded in `03_tool_comparation.Rmd` (intersection
46,302; CIRCexplorer-only 3,419; CIRIquant-only 23,739) and warns on mismatch.

### A8. Missing external resources

None of these is in any package received:

| Resource | Needed by |
|---|---|
| `rn8.yml` (CIRIquant config) | bash 14 — the stage refuses to guess |
| `rn8.refflat.nochr.txt` | bash 13 |
| `rat_primary_assembly_orderr8.fa` | bash 13, 15 — verify it is the GRCr8 assembly by checksum |
| STAR index `STAR_rn8/` | bash 13 |
| `CIRI2.pl` v2.0.6, `CIRI_AS_v1.2.pl` | bash 15 — not conda packages |

---

## B. For Bárbara and Sebastián — the preprocessing chain

### B1. Step [B] of the metadata chain is unrecovered

The `Babi analysis` chain runs
`extract_transcript_ids.sh` → `subset_merged_gtfs.sh` → `filter_out_lncRNAs.sh`
→ the metadata notebook. But those start by reading an `exons/` directory whose
*generating command is in neither the folder nor `bash_scripts.txt`*.

It was identified from the data, and the evidence is strong:

- FEELnc's 42,040 classified lncRNA transcripts are **100 % contained** in
  `exons/*.lncRNA.gtf` (set difference in that direction is exactly 0);
- `.lncRNA.gtf` and `.otherRNAs.gtf` are **disjoint**;
- all 46,030 transcript IDs are gffcompare `.tmap` `qry_id`s;
- the subset spans all 13 class codes, so it is not class-code driven;
- the files are exon lines only with the minimal attribute set.

That is the FEELnc_codpot signature: `.annotated.gtf` is the candidate set,
`.lncRNA.gtf` the lncRNA calls, `.otherRNAs.gtf` the complement. Two residuals
confirm it: 42,127 − 42,040 = **87** lncRNA transcripts the classifier drops as
orphans (no partner in the 1 kb/100 kb window), and
46,030 − 42,127 − 3,822 = **81** transcripts in the `noORF` bucket.

`scripts/bash/09b_feelnc_split_exon_gtfs.sh` implements this, and
`scripts/checks/check_feelnc_split_counts.sh` verifies it against the recorded
counts. **The check currently passes against the historical `exons/` directory.**

**Needed:** confirmation of the actual command.

### B2. Which FEELnc_codpot run was final

`provenance/longRNA_command_notebook/bash_scripts.txt` records four successive
runs on four different candidate GTFs:
`candidate_lncRNA_NE6_NE24.gtf` → `candidate_lncRNAs_NE6_NE24_oct.gtf` →
`cand_lncRNAs_NE6_NE24_annotated.gtf` → `cand_lncRNAs_ALLPE_annotated.gtf`.

**Needed:** which one produced the deposited results.

### B3. The metadata notebook on disk is not what produced the shipped tables

`Babi analysis/.Rhistory` shows `run_full_analysis()` being called with four
arguments (`config, ensembl, tx_attributes, exon_attributes`); the committed
`.qmd` defines it with three. Exon coordinates came from biomaRt in the earlier
version and from the GTF in the committed one.

This is one of the reasons `data/frozen/` exists. See `data/frozen/README.md`.

### B4. Two different server roots in one chain

Steps C and D ran under `/home/vmaracaja/shared-projects/SUZ_sr_2025`;
step E under `/home/vmaracaja/shared-projects/sebauz_sr`. The chain as received
never ran end to end in one place.

**Needed:** confirmation both pointed at the same data generation.

### B5. `isBest == 1` is a no-op on the shipped FEELnc file

Every row of `lncRNA_classes_CNE_6_24_HISAT2_1kb_100kb.txt` has `isBest == 1`
(42,040 of 42,040). The filter in R08/R14 is harmless but is **not** evidence
that filtering occurred. Do not describe it as a filtering step in the Methods.

---

## C. Annotation and references

### C1. The frozen metadata tables are Ensembl 114; the pipeline is pinned to 115

Resolved by decision, recorded here for the record. See
`docs/REFERENCE_AUDIT.md` and `data/frozen/README.md`. Both releases are GRCr8,
so coordinates are unaffected; gene names and biotypes are not.

**Needed before submission:** decide whether to regenerate at 115 and re-check
Figure 2B's counts, or keep 114 and say so in the Methods.

### C2. `Rattus_norvegicus.mRatBN7.2.115.gtf` does not exist

That filename appears in the starter package's `reference_manifest.tsv` and in
the circRNA addendum. Release 115 has no mRatBN7.2 build — the correct file is
`Rattus_norvegicus.GRCr8.115.gtf`. Already corrected throughout.

### C3. The circRNA branch used GRCr8.114

`CIRI_2_run_all_seba.sh:6` and `vis_run_ciri2_seba.sh:2` use
`Rattus_norvegicus.GRCr8.114.gtf`. Same assembly as the project lock, one
release earlier. An annotation-release gap, not an assembly conflict.

**Needed:** decide whether to rerun the circRNA branch at 115 or document the gap.

---

## D. Small RNA

### D1. Only the 24 h miRge command survives

`provenance/mirge/miRge_3.txt` has three variants of one 24 h invocation and no
6 h command, although `miR.Counts_C6_NE6_3v3_new.csv` is a required input.

### D2. `DESmetadata_24h.csv` is missing

Referenced as `-mdt`. Not in any package.

### D3. miRBase, not miRGeneDB

The run used `-db mirbase`. Confirmed by S. Urquiza on 2 Jul 2026. **The
manuscript Methods text must be corrected.** The exact miRBase release is still
unknown.

---

## E. Figure provenance

### E1. Manual assembly

Multi-panel figures were composed in Affinity and the Cytoscape v3.10.3 layouts
were adjusted by hand. The pipeline reproduces the analytical content, not the
page composition. See `docs/MANUAL_ARTIFACTS.md`.

**Needed:** the Affinity sources and the Cytoscape `.cys` or node coordinates,
if exact layout reproduction is wanted.

### E2. Two panels were written new

Neither had a standalone script in any package:

- Figure 3C circRNA → `scripts/R/13_figure3C_circRNA_heatmap.R` (audit item FIG3C_CIRC)
- Figure S3A-B mRNA/miRNA → `scripts/R/17_volcanoes_mRNA_miRNA.R` (audit item S3_MRNA_MIRNA)

They follow the conventions of their sibling panels, but they are new code and
their output has never been compared against the manuscript figures.

**Needed:** visual confirmation against the current drafts.

### E3. The lncRNA volcanoes have no fold-change cutoff

`07_Figure4_clean_volcanoes.R` sets `absfc_cutoff <- 0`: significance is
`padj` alone. R17 defaults to the same rule so the three S3 panels match.

**Needed:** confirmation that this is intended, and that the figure legend says so.

---

## F. Known bugs in the provenance — preserved, not reproduced

These are real defects in files under `provenance/`. They are kept verbatim
because that directory is the historical record; the active stages do not
repeat them.

| Where | Defect |
|---|---|
| `03_tool_comparation.Rmd:28-30` | `c3_cq <- length(c3_cq$chr)` overwrites the data frame with a scalar, then `c3_cq$circRNA_ID <- ...` is attempted on an atomic vector. **The chunk cannot run past line 30 as committed.** |
| `03_tool_comparation.Rmd:169,173` | `circRNA_sequences_de{6,24}_2.fa` are written with `write.table(sep="\t", col.names=TRUE)` — they are TSVs with a `.fa` extension, and are read back with `read.table(header=TRUE)`. |
| `01_circexp3.Rmd:159` | Writes the CIRCexplorer 6 h normalised matrix into the **CIRIquant** results folder. |
| `01_circexp3.Rmd:244` | Bypasses `paste0(path, ...)` with a hard-coded `D:` path. |
| `00_annotation_biomart.Rmd:47,129` | `getBM(values = ...)` with **no `filters =`**. biomaRt ignores `values`, so these chunks never retrieved what they appear to. `suz_getBM()` refuses this pattern. |
| `vis_run_ciri2_seba.sh:3` | Sample list reads `Q13 Q14 Q15 Q16` instead of `SQ13…SQ16`. |
| `vis_run_ciri2_seba.sh` | `rm ${sample}_vis/*.pdf` — deletes the PDFs it just generated. |
| `subset_merged_gtfs.sh` | No `case` arm for `CNE_24_STAR.annotated`; a correctly-named file would fall through to `*)` and be skipped with a warning. |
| `filter_out_lncRNAs.sh:26` | Requires GNU awk's 3-argument `match()`. Under mawk the filter silently passes every line through. Stage 09e now checks for gawk and aborts. |

---

## G. Resolved during this rebuild

| Was | Now |
|---|---|
| `12_STRICT_DE_coherent_circRNA_triplets.csv` listed as an orphan input | It is output `12_` of the cross-reference stage; R22 now runs before R23/R24 |
| Browser-download duplicates `name(1).csv` as first-preference inputs in 7 scripts | Removed; `find_input()` rejects the pattern outright |
| `base_dir` hard-coded in 17 scripts, 3× in one of them | optparse `--in-dir` / `--out-dir` everywhere; `check_hardcoded_paths.sh` enforces it |
| Unpinned biomaRt in 4 call sites | One `connect_ensembl()` pinned to the release-115 archive; `check_biomart_pinning.sh` enforces it |
| Raw `12_*` and three `Figure8_*` writing identical filenames into one directory | Every stage writes to its own `results/R<NN>_*/` |
| Curated KEGG allowlist duplicated between an inline tribble and a CSV | `metadata/curated_kegg_allowlist.tsv` |
| 4 orphan input families, discovered only by reading all 47 scripts side by side | `check_io_graph.R` — currently 0 orphans, 0 ordering violations |
