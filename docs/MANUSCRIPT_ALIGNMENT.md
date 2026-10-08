# Manuscript alignment check

Verification of this repository against the manuscript as of **v14**:

- `NE_NRVM_CDDIS_adapted_v14-SL.docx` (main text, figure legends, Methods)
- `NE_NRVM_CDDIS_Supplementary_Information_v14.docx`
- `Final_Supplementary_Tables_S1-S8_v14_2.xlsx` (19 sheets)
- Figures 1-6 and S1-S5 as rendered PNGs

Checked on 2026-10-08 against commit `0934e01`. Every item below cites the
manuscript line and the repository file it was compared with.

---

## Summary

| | Count |
|---|---|
| Figures whose producing stage is correctly mapped | 11 of 11 |
| Panels remapped because the previous map was wrong | 13 |
| Published numbers reproduced exactly by the pipeline | Figure 2B (6 of 6 checked) |
| Conflicts between the Methods text and the code | 6 |
| Internal inconsistencies within the manuscript | 4 |
| Repository defects found and fixed by this check | 5 |

The pipeline's **scope** matches the manuscript. The problems are in the
details: a stale figure panel, several Methods statements that do not describe
what the code does, and cross-reference numbering that drifted.

---

## 1. The strongest result: Figure 2B reproduces exactly

Running `scripts/R/06_transcriptome_overview.R` against the frozen
Ensembl-114 metadata tables produces, with no tuning:

| Biotype | Pipeline | Manuscript (Fig 2B and Results) |
|---|---:|---:|
| protein coding | 21,981 | 21,981 |
| annotated lncRNA | 16,149 | 16,149 |
| snoRNA | 1,778 | 1,778 |
| snRNA | 1,517 | 1,517 |
| novel lncRNA | 840 | 840 |
| miRNA | 204 | 204 |

All six match. Every other biotype in the published panel (pseudogene 361,
ig_v_gene 324, rRNA 263, tr_v_gene 180, …) also matches. Only the circRNA bar
(8,736) is absent, because that branch needs resources this repository cannot
ship.

This is strong evidence that the frozen Ensembl-114 tables in `data/frozen/`
are the correct inputs, and that the metadata chain reconstructed as bash
stages 09b–09e is faithful.

---

## 2. Figure map corrections

The previous map came from the August audit
(`04_METADATA/final_figure_script_map.tsv`) and predates the final figures. The
published legends disagree with it in three places. `metadata/figure_script_map.tsv`
has been rewritten from the v14 legends.

### 2.1 Figure 6 and Figure S5 were transposed

The audit said Figure 6A was the phase-aware coherence plot and Figure 6B-F
were five networks. The published figures say otherwise:

| | Published v14 |
|---|---|
| **Figure 6A** | Integrated **UP-DOWN-UP** network (Sustained + Late), 197 triplets |
| **Figure 6B** | Integrated **DOWN-UP-DOWN** network (Early + Sustained + Late), 2,396 triplets |
| **Figure S5A** | Phase-aware median log2FC trajectories with IQR |
| **Figure S5B-F** | The five phase-specific networks |

Consequence: `scripts/R/23_ceRNA_phase_aware_plots.R` produced a file named
`Figure6A_phase_timepoint_ceRNA_coherence_RECOMMENDED.pdf`. That is Figure S5A.
Renamed to `FigureS5A_phase_timepoint_ceRNA_coherence`.

### 2.2 Figure S4 is a single panel

The audit listed S4A (top TF rankings), S4B (lncRNA-mRNA pair counts by KEGG)
and S4C (top lncRNA hubs). The published Figure S4 has **only one panel**:
"Curated KEGG representation among significant lncRNA-mRNA correlations" —
what the audit called S4B. The TF-ranking and lncRNA-hub panels are not in the
final supplement; their underlying tables are in Supplementary Tables S3B and
S5D.

### 2.3 Figure S2 panel letters

Published: **A** exon-count distribution, **B** isoform complexity,
**C** transcript length. `scripts/R/07_lncRNA_characterization_gtf.R` named all
three `SuppFig2A_*`. Renamed to `FigureS2A_`, `FigureS2B_`, `FigureS2C_`.

---

## 3. Methods statements that do not match the code

These need a decision: change the manuscript, or change the pipeline.

### 3.1 The reference GTF named in the Methods does not exist — BLOCKING

> "Transcript assembly was performed per sample using StringTie v3.0.0 guided
> by **Rattus_norvegicus.mRatBN7.2.115.gtf**"

and the same filename again for FEELnc Module 1. **Ensembl release 115 has no
mRatBN7.2 build.** Releases 114 and 115 ship GRCr8; 113 and earlier ship
mRatBN7.2. The file is `Rattus_norvegicus.GRCr8.115.gtf`.

The Methods also pair that GTF with
`Rattus_norvegicus.GRCr8.dna.primary_assembly` for HISAT2 — a GRCr8 FASTA with
an "mRatBN7.2" GTF, which could not have been run as written.

Two further points on the same sentence:

- Ensembl publishes no separate `dna.primary_assembly` for GRCr8; the toplevel
  file is equivalent and is what the repository downloads.
- The circRNA Methods say both callers referenced "the Ensembl rat **rn7**
  genome". Every path in the recovered circRNA scripts is rn8/GRCr8
  (`rat_primary_assembly_orderr8.fa`, `rn8.refflat.nochr.txt`,
  `Rattus_norvegicus.GRCr8.114.gtf`, `STAR_rn8`, index `Rn8`). The rn7 claim is
  wrong.

**Recommended Methods text:** genome `Rattus_norvegicus.GRCr8.dna.toplevel.fa`
and annotation `Rattus_norvegicus.GRCr8.115.gtf` (Ensembl release 115), noting
that the circRNA branch used GRCr8 release 114. See `docs/REFERENCE_AUDIT.md`.

### 3.2 BSJ threshold: Methods say ≥2, code says ≥3 — BLOCKING

> "only circRNAs with **at least two BSJ counts** were carried forward"

The code is `readNumber > 2` in both notebooks
(`01_circexp3.Rmd:91`, `02_ciriquant.Rmd:85`), i.e. **at least three**.

This changes the circRNA inventory and therefore the 8,736 locus count in
Figure 2B and every circRNA count in Figure 3. `scripts/R/04_circRNA_merge_export.R`
exposes `--bsj-min` and logs the choice; it currently defaults to 3, matching
the code rather than the prose.

### 3.3 FEELnc classifier window: Methods say 200 bp, the data says 1 kb

> "assessed cis-context using a 200 bp to 100 kb window
> (--window=200 --maxwindow=100000)"

The classifier output every downstream stage consumes is named
`lncRNA_classes_CNE_6_24_HISAT2_1kb_100kb.txt`. A filename is evidence of what
ran; the Methods are not. Stage 09 defaults to 1000 and now carries the
conflict as a comment.

### 3.4 Salmon index — FIXED in the repo

> "using a genome-aware **decoy/gentrome** index"

Stage 10 built a transcriptome-only index, which gives different
quantifications. It now builds a proper gentrome (transcripts + genome) with a
decoy list, matching the Methods. `BUILD_DECOY=0` restores the old behaviour
and warns.

### 3.5 lncRNA volcano fold-change cutoff

Methods: lncRNAs require adjusted p < 0.05 **and |log2FC| ≥ 0.5**.
`07_Figure4_clean_volcanoes.R` (ported as R16) sets `absfc_cutoff <- 0` — padj
alone, no fold-change requirement. Figure S3's legend says the dashed lines
mark "the thresholds used for classification", so the panel and the Methods
should agree.

`scripts/R/17_volcanoes_mRNA_miRNA.R` now applies the class-specific Methods
thresholds (`--lfc-cut-mrna 1`, `--lfc-cut-mirna 0.5`). R16 is left as ported;
changing it would change a published panel.

### 3.6 "No additional interaction tiers were used"

Methods: "Only interactions with MFE ≤ -83.7 kJ/mol (-20 kcal/mol) were
retained; **no additional interaction tiers were used**."

The ported `12_..._ALL_RECORDS_FINAL.R` (R21) defines three tiers
(`energy_high -20`, `energy_moderate -15`, `energy_exploratory -10`) and sets
`minimum_pair_energy <- energy_exploratory`, so intermediate outputs retain
pairs down to -10 kcal/mol. The **final** networks do use -20
(`HIGH_CONFIDENCE_ceRNA_triplets_E_le_minus20.csv`, and the MFE20 naming
throughout R24/R25), so the figures are consistent with the Methods — but the
statement is not literally true of the script. Either soften the sentence or
state that the tiers exist only as intermediate audit outputs.

### 3.7 Confirmed as matching

| Methods | Code | |
|---|---|---|
| fastp `-f 10 -F 10`, `--cut_tail_mean_quality 28`, phred 30 | stage 02 and 13-15 | OK |
| HISAT2 `--dta-cufflinks`, `-S` | stage 04 | OK |
| Salmon `-l A --validateMappings` | stage 10 | OK |
| StringTie guided by `-G`, then merge, then gffcompare | stages 06-07 | OK |
| FEELnc `--monoex=1`, three modules | stage 09 | OK |
| miRge3.0 rat miRBase, adapters `AGATCGGAAGAGCACACGTCT` / `GTTCAGAGTTCTACAGTCCGACGATC` | stage 11 | OK |
| miRWalk score ≥ 0.95, 3'UTR/5'UTR/CDS, miRDB support | 12 input files in `metadata/external_inputs.tsv` | OK |
| circRNA DE: nominal p < 0.05, \|log2FC\| ≥ 0.5 | R20 (`alpha_nominal`, `lfc_cut`) | OK |
| mRNA DE: padj < 0.05, \|log2FC\| ≥ 1 | R02 | OK |
| Correlation: Pearson, \|rho\| ≥ 0.8, p < 0.05, VST, Hmisc `rcorr` | R09 (`rho_cutoff 0.80`) | OK |
| GRN: TFLink, expressed > 0 in ≥ 1 replicate, mean normalised count > 10 | R05 | OK |
| IntaRNA v3.4.1, MFE ≤ -20 kcal/mol | stages 12/16, R21 | OK |
| **N = 3 biological replicates per condition** | `metadata/longRNA_samples.tsv` marks replicate 3 excluded | OK |

---

## 4. Internal inconsistencies within the manuscript

### 4.1 Figure 3B's lncRNA panel is stale — the figure needs regenerating

The Results text and Figure S3C agree with each other and **disagree with
Figure 3B** for lncRNA only:

| lncRNA program | Results text | Figure S3C (derived) | Figure 3B (published) |
|---|---|---|---|
| Early | 148 (62 up, 86 down) | 148 (62, 86) | **145 (61, 84)** |
| Sustained | 258 (86 up, 172 down) | 258 (86, 172) | **275 (95, 180)** |
| Late | 336 (109 up, 227 down) | 336 (109, 227) | **328 (104, 224)** |

Derivation from Figure S3C, summing the annotated and novel matrices:
Early = (47+15) up and (70+16) down; Sustained = (66+20) up and (128+44) down;
Late = (92+17) up and (163+64) down. Every one matches the Results text exactly.

mRNA, miRNA and circRNA reconcile perfectly across Figure 3A, 3B, S3C and the
text. **Only the lncRNA panel of Figure 3B is wrong**, and it is the figure
that is wrong, not the text. Regenerate it with
`scripts/R/10_figure3B_temporal_programs.R`.

### 4.2 Supplementary table cross-references disagree between legends and text

| | Results text cites | Figure legend cites |
|---|---|---|
| GSEA | Table S2 | Figure 4: "Supplementary Tables 4-5" |
| TF rankings | Table S3 | Figure 4: as above |
| lncRNA GBA | Table S5 | Figure 5: "Supplementary Tables 7-8" |
| miRNA-mRNA | Table S6 | Figure 5: as above |

The Results text is correct against the workbook's actual sheet names
(`S2_KEGG_GSEA`, `S3A_TF_pathway`, `S3B_TF_rankings`, `S5B_exact_GBA`,
`S6A_unique_pairs`). The figure legends need fixing. The author already noted
this at the top of the draft: *"Tablas y Figuras no siguen orden correlativo"*.

`metadata/supplementary_tables_map.tsv` maps all 19 sheets to producing stages.

### 4.3 A Supplementary Table 10 is cited but does not exist

> "the synchronized strict set comprised 2,593 triplets … (**Supplementary Table 10**)"

The workbook contains S1-S8. This should be Table S8 (`S8A_strict_triplets`).

### 4.4 Leftover review notes in the draft

Three Spanish editorial notes are still in the text and must come out before
submission: a comment after the Abstract ("Se ve muy general: cuál es el
descubrimiento de pathway más importante y novedoso"), "REFERENCIAS CON
NUMEROS", and "Tablas y Figuras no siguen orden correlativo". There is also a
typo at the start of the Results: "Differential.expression".

### 4.5 Placeholders still to fill

`[FINAL GITHUB REPOSITORY URL]` appears twice (Methods and Data availability)
and `[SEQUENCING ACCESSION(S)]` once. The repository URL is
`https://github.com/integrativebioinformatics/suz-nrvm-hypertrophy2026`.

---

## 5. Repository defects this check found, now fixed

| Defect | Impact |
|---|---|
| 16 ported stages rebound `out_dir` to a bare relative path **after** the CLI set it | Every one of those stages wrote into the current working directory instead of `--out-dir`. All now nest under `.suz_out_root`. The porter was fixed so future ports cannot reintroduce it. |
| Salmon built a transcriptome-only index | Did not match the Methods; different quantifications. Now builds a decoy-aware gentrome index. |
| `REPO_ROOT` was used but never defined in stage 10 | Would have aborted under `set -u`. |
| R07 named all three S2 panels `SuppFig2A_*` | Panel letters now follow the published figure. |
| R23 named its output `Figure6A_*` | It is Figure S5A. |

---

## 6. Scope: nothing in the manuscript is unaccounted for

Every panel in the published figures maps to a stage, an experimental source,
or a documented gap:

- **Produced by the pipeline:** Figure 2B-C, 3A-C, 4A-B, 5A-C, 6A-B, S2A-C,
  S3A-C, S4A, S5A-F.
- **Experimental, not transcriptomic:** Figure 1B-C, S1A-C (microscopy, ImageJ,
  RT-qPCR, GraphPad Prism).
- **Manual diagrams:** Figure 1A, 2A.
- **Manual layout on top of pipeline output:** Figure 6A-B and S5B-F
  (Cytoscape v3.10.3), and all multi-panel assembly (Affinity). See
  `docs/MANUAL_ARTIFACTS.md`.

Conversely, no pipeline stage is orphaned: all 26 R stages and 21 bash stages
feed a published panel, a supplementary table, or a documented prerequisite.

**The one scope gap is the circRNA branch**, which gates Figure 2B's circRNA
bar, all circRNA panels in Figures 3 and S3, Supplementary Tables S7A-B, and
the 9 circRNA triplets in Figure 6B. It depends on `rn8.yml`, `rn8.refflat`, a
STAR index, `CIRI2.pl` and `CIRI_AS.pl`, none of which were supplied, plus the
RECONSTRUCTED stage R04. See `docs/OPEN_ITEMS.md` sections A7 and A8.
