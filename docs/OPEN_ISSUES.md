# Open issues

Everything unresolved, missing, or waiting on a decision. **Ordered the way the
pipeline runs**, so each item sits next to the stage it affects.

Nothing here is a guess: each item cites the file and line it comes from. What
was changed and why is in
[`CHANGES_FROM_ORIGINALS.md`](CHANGES_FROM_ORIGINALS.md); this file is only
what is *not* settled.

---

## Index by owner

A round of decisions on 14 Oct 2026 (B. Borges) closed eleven items. What is
left is blocked on someone else or on executing a run — with two exceptions,
noted below, deliberately left open.

**Unrecoverable — nobody can supply these**
§4.1 the `exons/` split command · §4.2 which codpot run was final · §4.3 the
two server roots · §5.1 the 6 h miRge command and `DESmetadata_24h.csv` ·
§8.1 which library "rep3" was · §8.2 the notebook variant behind the frozen
tables

**S. Urquiza — still open by choice**
§3.4 StringTie flags · §7.2 ComBat_seq

**S. Urquiza — inputs to supply**
§5.2 small-RNA FASTQ checksums · §6.1 `TFLink_All.tsv` and the three Cytoscape
GRN exports · §6.2 the 12 miRWalk files and miRBase `mature.fa` · §12.3 the
Affinity and Cytoscape sources

**Needs a run, not a decision**
§7.5 the circRNA reconstruction · §11.1 the two new panels against the
published figures · §15 loose ends · "Closes on the re-run" below

**S. Urquiza — manuscript text**
§14, now twenty items

### What blocks what

| Blocked | By |
|---|---|
| Figure 4A-B and S4A | §6.1 |
| Figure 5B-C | §6.2 |
| Confirming the two new panels (Fig 3C circRNA, S3A-B mRNA/miRNA) | §11.1 — needs a run |
| Confirming the circRNA reconstruction | §7.5 — needs a run |

**The circRNA module is no longer externally blocked.** Every resource §7.4
listed is now either built by the pipeline (STAR index, refFlat, BWA index,
`.fai`, the CIRIquant config) or installed from bioconda (`ciri2`, `ciri-as`).
It is gated only on being run.

**Not on this list: the reference release.** Release 115 with the GRCr8 primary
assembly is a settled forward decision (B. Borges, 8 Oct 2026). Everything,
including the circRNA analysis, will be re-run at 115 and the results compared.

**Also not on this list: the three accepted analysis changes.** §7.1 (BSJ ≥ 2),
§7.3 (CIRI-AS on all samples) and §12.1 (lncRNA fold-change cutoff) are
decided, but they change published numbers — see each one for what has to be
re-read off the new run.

---

## Stage 00 — references

### 3.1 and 3.3 — CLOSED by decision: one assembly for every branch

Previously two items: the concatenation order behind
`rat_primary_assembly_orderr8.fa` was unrecorded (3.1), and that file was never
supplied so could not be checksum-verified (3.3).

**Decision (B. Borges, 14 Oct 2026): the circRNA branch uses the stage 00
primary assembly, the same file as the long-RNA branch.** Stages 12b, 14 and 15
all default to it, so there is nothing left to reconcile or verify.

What this costs: the original concatenation order is not reproduced, so BAM
header order may differ from the 2023 run. That affects header order and some
tools' output listings, **not alignment and not coordinates** — both are GRCr8
either way. In exchange, every branch is provably on one genome, which is what
made §3.3 a question in the first place.

### 3.2 Reference checksums — CLOSED, verified at download

Previously: every `checksum` cell in `metadata/reference_manifest.tsv` was
`TO_FILL`, to be filled by hand after a run.

Stage 00 now does both halves automatically:

- **Verifies each download against Ensembl's own `CHECKSUMS` file**, which
  every FTP directory carries in BSD `sum` format. A mismatch is fatal with the
  expected and actual values printed; a filename Ensembl does not list only
  warns, since not every file appears in every release.
- **Records the md5** of the finished FASTA and GTF into
  `reference_manifest.tsv`, matched on `final_file_name`.

So the cells populate on the first run, and a truncated download now fails
loudly instead of producing a silently short genome.

## Stage 06 — StringTie

### 3.4 The GRCr8 StringTie invocation was never recorded

Every StringTie command in `bash_scripts.txt` is rn7-era. `--conservative`
appears in five of them (`:107`, `:109`, `:152`, `:154`, `:156`) and two early
single-sample tests add `-f 0.05 -m 200 -t -s 3.0`. Which, if any, carried
into the final run is unknown.

The stage defaults to plain guided assembly, which is what the Methods
describe. `STRINGTIE_EXTRA="--conservative"` adds the flags back.

**Needed:** the flags used for the final run, from S. Urquiza.

## Stage 09 — FEELnc

### 4.2 Which `FEELnc_codpot` run was final

`bash_scripts.txt` records four successive runs on four different candidate
GTFs:

```
candidate_lncRNA_NE6_NE24.gtf
  → candidate_lncRNAs_NE6_NE24_oct.gtf
    → cand_lncRNAs_NE6_NE24_annotated.gtf
      → cand_lncRNAs_ALLPE_annotated.gtf
```

Resolved for the GRCr8 run — `:293` uses
`merged_CNE624_8_rn8_H.annotated.gtf`, which supersedes all four — but the
**Module 2 command for that run is missing entirely**: the heading `M2` at
`:295-297` has nothing under it. Only Module 1 and Module 3 were recorded.

The stage carries the parameters that were identical across all four earlier
runs. That is reasonable but unverified for the final one.

**Needed:** the final codpot command, from S. Urquiza.

## Stage 09b — the reconstructed `exons/` split

### 4.1 The original command is unrecovered — the single largest gap upstream

The metadata chain begins by reading an `exons/` directory whose generating
command is in neither the metadata folder nor `bash_scripts.txt`.

The reconstruction's evidence is strong and is set out in
`CHANGES_FROM_ORIGINALS.md` (stage 09b). `check_feelnc_split_counts.sh`
reproduces 46,030 / 42,127 / 3,822 with both invariants and **currently passes
exactly** against the historical directory. **Confidence: high** — but a
passing count check is not the command.

**Needed:** confirmation of what was actually run, before the repository is
frozen.

## Stages 09c-09e — the transcript-metadata chain

### 4.3 Two different server roots in one chain

Steps C and D ran under `/home/vmaracaja/shared-projects/SUZ_sr_2025`; step E
under `/home/vmaracaja/shared-projects/sebauz_sr`. The chain as received never
ran end to end in one place.

**Needed:** confirmation both pointed at the same data generation.

### 4.4 `isBest == 1` is a no-op — CLOSED, it is a Methods fix

Every row of `lncRNA_classes_CNE_6_24_HISAT2_1kb_100kb.txt` has
`isBest == 1` — 42,040 of 42,040. The filter in R08/R14 is harmless but is
**not** evidence that filtering occurred.

This was never an open question; it has a definite answer and only needed
writing down. Tracked as §14.7: do not describe it as a filtering step. The
filter stays in the code, because on a different FEELnc run it would do
something.

## Stage 11 — miRge3.0

### 5.1 Only the 24 h command survives, and its metadata file is missing

`archive/raw_scripts_original/miRge_scripts/miRge_3.txt` has three variants of
one 24 h invocation and **no 6 h command**, although
`miR.Counts_C6_NE6_3v3_new.csv` is a required input to R03.

`DESmetadata_24h.csv`, referenced as `-mdt`, is in no package received.

The exact **miRBase release** is also unrecorded. That `-db mirbase` was used
rather than miRGeneDB is confirmed (S. Urquiza, 2 Jul 2026); the release is
not.

**Needed from S. Urquiza:** the 6 h command, `DESmetadata_24h.csv`, and the
miRBase release.

### 5.2 Small-RNA FASTQ checksums

16 rows in `metadata/fastq_manifest.tsv` are still `TO_FILL` — paths, sizes
and MD5s. The 32 long-RNA rows are complete.

## Stages 13-16 — the circRNA branch

### 7.1 BSJ read-support threshold — BLOCKING

| | |
|---|---|
| Methods v14 | "only circRNAs with **at least two** BSJ counts were carried forward" |
| Code | `count_matrix[count_matrix$readNumber > 2,]` — `01_circexp3.Rmd:91`, `02_ciriquant.Rmd:85` |

`> 2` means **at least 3**. These are not the same rule and they give different
inventories, so Figure 2B's 8,736 loci and every circRNA count in Figures 3
and S3 move with the choice.

**CLOSED by decision (B. Borges, 14 Oct 2026): use >= 2, the Methods rule.**
`--bsj-min` now defaults to 2 and logs the choice on every run; `--bsj-min 3`
reproduces the historical notebook code for comparison.

Reasoning: the published text commits to "at least two", that is the number a
reader will check, >= 2 is the common convention, and `readNumber > 2` reads
like an off-by-one against the author's own stated intent rather than a
deliberate choice.

**Consequence to track:** this changes the circRNA inventory, so Figure 2B's
8,736 locus count and every circRNA count in Figures 3 and S3 must be re-read
off the new run, not carried over. No Methods text changes.

### 7.2 Was ComBat_seq used in the final circRNA DE?

| | |
|---|---|
| Historical notebooks | `ComBat_seq`, then `pvalue < 0.05 & abs(log2FoldChange) > 1` (`01_circexp3.Rmd:169-170, 253-254`; `02_ciriquant.Rmd:158-159, 273-274`) |
| Final figure script | nominal `p < 0.05` and `abs(log2FC) >= 0.5`, **no ComBat_seq** (`25_recreate_circRNA_volcanoes_with_all_NS.R:26-27`) |
| Methods | silent |

`--use-combat` exists so the historical route can be reproduced for
comparison; it is **off** by default.

Note that `batch` maps 1:1 to replicate index (`01_circexp3.Rmd:124`, `:212`),
so dropping replicate 3 removes **Batch3 entirely**. The exclusion is
batch-aligned, not random — which matters if ComBat_seq is retained anywhere.

**Needed:** confirmation that the final tables came from the no-ComBat route.

### 7.3 CIRI-AS covers SQ13-SQ16 only

`CIRI_2_run_all_seba.sh:2` is `for m in {13..16}`. There is no CIRI2/CIRI-AS
loop for H1-H4 or SQ1-SQ8 anywhere in any package.

**CLOSED by decision (B. Borges, 14 Oct 2026): run every sample.** Stage 15
now defaults to the whole sample sheet, and `--only
SQ13_ln,SQ14_ln,SQ15_ln,SQ16_ln` reproduces the original four-sample loop.
`--yes-widen` is gone.

Reasoning: restricting sequence reconstruction to four of sixteen libraries has
no methodological justification. The loop reads like it was written for one
batch and never extended, and the Methods describe reconstruction for the
analysis as a whole.

**Consequence to track:** the four circRNA IntaRNA target FASTAs, and through
them bash 16 and the circRNA ceRNA triplets, are now built from all samples.
The Methods sentence on CIRI-AS becomes true as written.

### 7.4 External resources the circRNA module needs

None of these was in any package received. Two are now built rather than
supplied; four are still outstanding.

| Resource | Needed by | Status |
|---|---|---|
| STAR index `STAR_rn8/` | bash 13 | **Resolved.** Stage 12b runs `STAR --runMode genomeGenerate` on the stage 00 primary assembly. No build command was ever recorded, so this is new code, not a recovery |
| `rn8.refflat.nochr.txt` | bash 13 | **Resolved.** Stage 12b derives the refFlat from the Ensembl GTF with `gtfToGenePred`. Ensembl names carry no `chr` prefix, which is what "nochr" recorded |
| `rn8.yml` — CIRIquant config | bash 14 | **Resolved.** With no `--config`, stage 14 generates an equivalent one from the active reference plus the BWA index (12b) and HISAT2 index (04), and prints it. Schema per the CIRIquant cookbook: `bwa_index` and `hisat_index` are prefixes, and the FASTA needs a `.fai` beside it, which 12b writes |
| `rat_primary_assembly_orderr8.fa` | bash 13, 15 | **Resolved by decision** — every branch uses the stage 00 primary assembly. See §3.1/§3.3 |
| `CIRI2.pl` v2.0.6 | bash 15 | **Resolved.** bioconda ships `ciri2`, whose current release is 2.0.6 — the exact version the Methods state. Pinned in `envs/environment_circrna.yml`; stage 15 finds it on `PATH` |
| `CIRI_AS_v1.2.pl` | bash 15 | **Resolved.** bioconda ships `ciri-as`, current release 1.2, again the exact Methods version |

So **nothing external is outstanding for the circRNA module.** It is now gated
only on running it — see §7.5.

`ciri-full` also bundles CIRI2 and CIRI-AS, but at its own versions and with a
Java dependency, and CIRI-full itself is commented out in the original script
(`CIRI_2_run_all_seba.sh:93,95-97`), as is CIRI-vis (`:105`). The two
standalone packages are both more faithful and smaller.

### 7.5 The circRNA reconstruction is unconfirmed

R04 + R04b rebuild the count-matrix assembly, tool merge and export from
`03_tool_comparation.Rmd`. They assert the recorded locus counts —
intersection 46,302, CIRCexplorer-only 3,419, CIRIquant-only 23,739 — and warn
rather than error on mismatch.

**Confidence: medium.** The FASTA export additionally depends on CIRI-AS
sequences, which §7.3 puts in doubt.

**Needed:** confirmation that the reconstruction matches what was actually run.

### 13.1 "CIRCexplorer3" — CLOSED, the Methods were right

Previously listed here as a naming discrepancy. It is not one.

**CLEAR and CIRCexplorer3 are the same tool.** The YangLab repository is
titled "CLEAR/CIRCexplorer3". There is no separate CIRCexplorer3 binary:
CLEAR declares CIRCexplorer2 >= 2.3.6 as a dependency and calls its `parse`
and `annotate` subcommands itself, adding `circ_quant` for FPBcirc, FPBlinear
and CIRCscore. The original script shows it directly — `CLEAR/src/run3.py` and
`CLEAR/src/circ_quant.py` wrapped around the two CIRCexplorer2 calls, under a
directory named `CIRCexp3`.

Two smaller points survive, folded into §14:

- **"v1.0.0" is invented.** CLEAR carries no version number. Cite the paper and
  the installed commit instead (§14.14).
- **The CLEAR paper is not cited** (§14.15): Ma XK, Wang MR, Liu CX, Dong R,
  Carmichael GG, Chen LL, Yang L. *A CLEAR pipeline for direct comparison of
  circular and linear RNA expression.* doi:10.1101/668657.

Worth knowing, though not a defect: `circ_quant` runs, so FPBcirc, FPBlinear
and CIRCscore are produced, but differential expression uses the raw
`readNumber` back-splice counts from `annotate`. CLEAR's distinguishing output
is computed and then not used.

### 13.2 The circRNA branch ran at GRCr8.114

`CIRI_2_run_all_seba.sh:6` and `vis_run_ciri2_seba.sh:2` use
`Rattus_norvegicus.GRCr8.114.gtf` — same assembly as the project lock, one
annotation release earlier.

**Not a decision — closes on the re-run.** Release 115 is settled, and the
circRNA branch is re-run with everything else, at which point this gap simply
does not exist. Listed under "Closes on the re-run" below.

Nothing to do but note it in the Methods: the published numbers came from 114.

## R01 — transcriptome metadata

### 8.2 The committed notebook is not quite the code that produced the shipped tables

The `.Rhistory` of the original analysis shows `run_full_analysis()` called
with four arguments (`config, ensembl, tx_attributes, exon_attributes`); the
committed notebook defines it with three. Exon coordinates came from a second
biomaRt round-trip in the earlier version and from the GTF in the committed
one. An intermediate variant produced the shipped files, and that variant was
not preserved.

This is one of the three reasons `data/frozen/` exists. See
`data/frozen/README.md`.

### 8.3 The frozen tables are Ensembl 114; the pipeline is pinned to 115

Both releases are GRCr8, so coordinates are unaffected — but gene names and
biotypes move between releases, and those feed Figure 2B's class counts and the
annotated-versus-novel lncRNA split.

**Not a decision — closes on the re-run.** Regenerate with R01, diff against
the frozen copies, and re-read Figure 2B's counts, Figure 2C and the Figure S2
lncRNA inventory off the new tables. Report the delta to S. Urquiza before
replacing `data/frozen/`. Listed under "Closes on the re-run" below.

## R02 — mRNA/lncRNA DESeq2

### 8.1 Which physical library is "rep3"?

`Integrated_Bulk_v3_6_24h.R` reads Salmon output from directories named
`Ctrl6_rep1_quant` … `NE24_rep4_quant`. **No file in any package maps those
aliases to `H1_ln`, `SQ13_ln`, …** So although replicate 3 is demonstrably
excluded from every final analysis, *which library* that was cannot be stated.

`metadata/longRNA_samples.tsv` has a `salmon_quant_alias` column, currently
`TO_CONFIRM` for every row.

**Needed:** the alias → `sample_id` mapping, and confirmation that marking
those four rows `no_excluded_from_final_DE` is right. Ideally the Methods
should say which library was excluded and why.

## R05 — TF and KEGG integration

### 6.1 Four inputs missing

| Missing | Note |
|---|---|
| `TFLink_All.tsv` | plus its version or download date |
| `GRN_NE06_*` edge export | manual Cytoscape export, not reproducible from the pipeline |
| `GRN_NE24_*` edge export | as above |
| `Merged_Network_*` edge export | as above |

Blocks Figures 4A-B and S4A.

## R13, R17 — the two new panels

### 11.1 Neither has been compared against the manuscript figure

Neither panel had a standalone script in any package; both were flagged `OPEN`
in the v3 audit and were written for this repository:

- Figure 3C circRNA → `scripts/R/13_figure3C_circRNA_heatmap.R` (`FIG3C_CIRC`)
- Figure S3A-B mRNA/miRNA → `scripts/R/17_volcanoes_mRNA_miRNA.R`
  (`S3_MRNA_MIRNA`)

They follow the conventions of their sibling panels, but they are new code and
their output has never been checked against the published figures.

**Needed:** visual confirmation against the current drafts.

## R14 — miRNA/mRNA integration

### 6.2 Three inputs missing

| Missing | Note |
|---|---|
| the 12 `miRWalk_miRNA_Targets_*.csv` files | also read by R22 |
| miRBase `mature.fa` | plus the exact release — see §5.1 |
| `miRDB target prediction data_4miRNAs_correlation.xls.xlsx` | optional, historical only |

Blocks Figures 5B-C.

## R16 — lncRNA volcanoes

### 12.1 padj alone, or padj plus |log2FC| ≥ 0.5?

| | |
|---|---|
| Methods v14 | lncRNAs: adjusted p < 0.05 **and \|log2FC\| ≥ 0.5** |
| Code | `07_Figure4_clean_volcanoes.R` sets `absfc_cutoff <- 0` — padj alone |

**CLOSED by decision (B. Borges, 14 Oct 2026): apply the Methods rule.**
R16's `absfc_cutoff` now defaults to 0.5 and is exposed as `--lfc-cut`;
`--lfc-cut 0` reproduces the original padj-only behaviour.

Reasoning: R17 already applies class-specific thresholds for mRNA and miRNA, so
leaving R16 on padj alone meant Figure S3's three panels were classified by two
different rules, and the legend's claim about "the thresholds used for
classification" was false for one of them.

**A second defect surfaced while making the change.** R16 drew its dashed
vertical lines at a hard-coded **±1** while the cutoff was 0 — so the lines
never marked the classification rule at all, under either threshold. They now
track `absfc_cutoff`, and are omitted entirely when it is 0.

**Consequence to track:** the Figure S3A-B lncRNA panels must be regenerated,
and their significant-lncRNA counts will drop, since a fold-change requirement
now applies where none did.

## R21-R25 — ceRNA networks and export

### 12.2 "No additional interaction tiers were used"

The Methods say "Only interactions with MFE ≤ −83.7 kJ/mol (−20 kcal/mol) were
retained; no additional interaction tiers were used." The ceRNA script defines
three tiers and sets the floor at −10, so intermediate outputs retain weaker
pairs.

The **final** networks do use −20, so the figures match the Methods; the
sentence is just not literally true of the script.

**CLOSED by decision (B. Borges, 14 Oct 2026): keep the tiers, fix the
sentence.** The tiers produce audit output that is worth having, and nothing
published depends on the −10 floor. It is a wording fix, now §14.16 — state
that intermediate audit outputs retain weaker pairs and only the −20 tier
reaches the figures. No code changes.

### 12.3 Figures are not pixel-reproducible

Multi-panel assembly happened in Affinity and the Cytoscape v3.10.3 layouts
were adjusted by hand. The pipeline reproduces the analytical content, not the
page composition or the node coordinates. Full detail in
`CHANGES_FROM_ORIGINALS.md` Part 6.

**Needed, if exact reproduction is wanted:** the Affinity sources
(`.afdesign` or exported SVG/PDF) and the Cytoscape `.cys` session or the final
node x/y coordinates plus the visual style.

---

## 14. Manuscript fixes needed

| # | Fix | From |
|---|---|---|
| 14.1 | `Rattus_norvegicus.mRatBN7.2.115.gtf` → `Rattus_norvegicus.GRCr8.115.gtf` (×2) — the named file has never existed | stage 00 |
| 14.2 | circRNA reference "Ensembl rat rn7 genome" → GRCr8 | stages 13-15 |
| 14.3 | FEELnc window "200 bp to 100 kb" → "1 kb to 100 kb" | stage 09 |
| 14.4 | miRge flags "(`-lib` and `-rat`)" → `-lib`, `-on rat`, `-db mirbase` | stage 11 |
| 14.5 | Say n = 3 per group, and ideally which library was excluded | §8.1 |
| 14.6 | Describe circRNA DE as exploratory — nominal p, not FDR | R20 |
| 14.7 | Do not describe the FEELnc `isBest` filter as a filtering step | §4.4 |
| 14.8 | **Regenerate Figure 3B's lncRNA panel** — published 145/275/328, correct 148/258/336 | R10 |
| 14.9 | Figure 4 and 5 legends cite the wrong supplementary tables | — |
| 14.10 | "Supplementary Table 10" → S8 | — |
| 14.11 | Remove the three Spanish editorial notes; fix "Differential.expression" | — |
| 14.12 | Fill `[FINAL GITHUB REPOSITORY URL]` ×2 and `[SEQUENCING ACCESSION(S)]` | — |
| 14.13 | Soften "no additional interaction tiers were used" | §12.2 |
| 14.14 | **"CIRCexplorer3 v1.0.0" -> drop the version.** CLEAR publishes none; give the installed commit instead | §13.1 |
| 14.15 | **Cite the CLEAR paper**: Ma XK et al., doi:10.1101/668657 | §13.1 |
| 14.16 | Rewrite "no additional interaction tiers were used" — say intermediate audit outputs retain weaker pairs and only the −20 tier reaches the figures | §12.2 |
| 14.17 | **Add CIRI2 to the Methods.** The text names CIRI-AS but not CIRI2, and CIRI-AS cannot run without it: `CIRI_2_run_all_seba.sh:90` passes CIRI2's `.ciri` call set via `-C`. Cite Gao et al. 2018 for CIRI2 alongside the existing CIRI-AS citation | new |
| 14.18 | **`bwa mem -T 19` is attributed to the wrong tool.** The text puts it under CIRIquant; the only recorded `-T 19` is `CIRI_2_run_all_seba.sh:84`, feeding CIRI2/CIRI-AS. CIRIquant runs BWA internally with its own settings | new |
| 14.19 | Say that circRNA DE uses the raw `readNumber` BSJ counts from CIRCexplorer2 `annotate`. `circ_quant` runs, so CIRCscore is computed, but nothing downstream uses it — a reader told the pipeline is CLEAR will look for it | §13.1 |
| 14.20 | The circRNA reference was GRCr8.114 in 2023 and is GRCr8.115 on the re-run. One sentence, since both are GRCr8 | §13.2 |

---

## Closes on the re-run

Not decisions and not blockers — they resolve themselves once the pipeline is
executed at release 115, and are listed only so nobody re-opens them as
questions.

| # | Item |
|---|---|
| §8.3 | The frozen Ensembl-114 metadata tables get regenerated at 115; diff them and re-read Figure 2B / 2C / S2 off the new output |
| §13.2 | The circRNA branch moves from GRCr8.114 to GRCr8.115 with everything else |
| §3.2 | `reference_manifest.tsv` checksums populate automatically on the first stage 00 run |
| §7.1, §7.3, §12.1 | The three accepted analysis changes take effect; the affected figures and counts must be re-read, not carried over |
| §15.2 | `envs/*.lock.yml` are written by `create_envs.sh`, fixing the tool versions for the Methods |
| §15.8 | The Docker image gets its first real build |

---

## 15. Repository loose ends

| # | Item |
|---|---|
| 15.1 | All `checksum` cells in `metadata/reference_manifest.tsv` are `TO_FILL` (§3.2) |
| 15.2 | **The environments pin no versions**, so `bash scripts/commands/create_envs.sh` installs current releases rather than the versions the Methods record (those are trailing comments in the `.yml` files). This is deliberate — the pipeline is being re-run at Ensembl 115 — but it means the Methods' tool versions must be rewritten from `envs/*.lock.yml` once the environments are built, and `metadata/software_tool_inventory.tsv` updated with them. Several rows there are still `TO_FILL` |
| 15.3 | `sessionInfo()` from the machine that produces the deposited figures is not recorded |
| 15.4 | **No full pipeline run has been executed.** The repository is verified structurally by the eight self-checks, and numerically only for Figure 2B |
| 15.5 | 16 small-RNA rows in `metadata/fastq_manifest.tsv` are `TO_FILL` (§5.2) |
| 15.6 | **CLEAR has no conda package** and installs from source into `env_circrna` (`git clone` + `python ./setup.py install`). Nothing pins or checksums it, so the version used is whatever `master` is on the day it is installed. Record the commit hash |
| 15.7 | **The environments have never been solved.** `clear` was listed as a bioconda dependency and does not exist there, and `mirge3=3.0` pinned a version that does not exist — both would have failed `conda env create`. Fixed, but no solve has been run, so more may surface. `bash scripts/commands/create_envs.sh --dry-run` solves without installing |
| 15.8 | **The image has never been built.** `Dockerfile` is statically validated only — every `COPY` source and every script it invokes exists, and all 8 `RUN` blocks parse as shell — but there is no Docker daemon on the machine it was written on. The first `docker build` is the real test |
| 15.9 | `Dockerfile` uses `condaforge/miniforge3:latest` and CLEAR's `master`, so neither the base image nor CLEAR is pinned. After the first good build, pin both: `docker buildx imagetools inspect` gives the base digest for `--build-arg BASE_IMAGE`, and the image records CLEAR's at `/opt/CLEAR/INSTALLED_COMMIT` for `--build-arg CLEAR_REF` |

---

## 16. Explicitly out of scope

| Item | Why |
|---|---|
| The STAR branch of the metadata analysis | HISAT2 only, confirmed by S. Urquiza |
| The Quarto notebooks of the metadata analysis | Logic ported to R01; the originals stay outside the repository by decision. `archive/ARCHIVE_NOTES.md` records where |
| Pixel-identical figures | Affinity assembly and manual Cytoscape layout sit between pipeline output and the manuscript (§12.3) |
| Multi-GB intermediates | Merged GTFs, `*_exon_lengths.csv` (64 MB), `*_feelnc_validated.csv` (26 MB) and `lncRNA_classes_*.txt` (27 MB) are regenerable by stages 06, 07, 09 and R01, and are manifest-only |
| Figure content validation against the drafts | Needs the final figures; only Figure 2B has been numerically verified |

---

## 17. Resolved during the rebuild

Kept so the same ground is not re-covered.

| Was | Now |
|---|---|
| `12_STRICT_DE_coherent_circRNA_triplets.csv` listed as an orphan external input | It is output `12_` of the cross-reference stage; R22 now runs before R23/R24 |
| Browser-download duplicates `name(1).csv` as first-preference inputs in 7 scripts | Removed; `find_input()` rejects the pattern outright |
| `base_dir` hard-coded in 17 scripts, three times in one of them | optparse `--in-dir` / `--out-dir` everywhere, enforced by `check_hardcoded_paths.sh` |
| Unpinned biomaRt in 4 call sites | One pinned `connect_ensembl()`, enforced by `check_biomart_pinning.sh` |
| Raw `12_*` and three `Figure8_*` writing identical filenames into one directory | Every stage writes to its own `results/R<NN>_*/` |
| The curated KEGG allowlist duplicated between an inline tribble and a CSV | `metadata/curated_kegg_allowlist.tsv`, read by both stages |
| 4 orphan input families | `check_io_graph.R` reports 0 orphans and 0 ordering violations |
| 16 stages rebinding `out_dir` after the CLI set it, writing to the working directory | All nest under `.suz_out_root` |
| Salmon building a transcriptome-only index | Decoy-aware gentrome index, matching the recorded command |
| Stage 00 downloading `dna.toplevel.fa` | Primary assembly, concatenated from the per-chromosome parts |
| `REPO_ROOT` used but never defined in stage 10 | Defined; the stage would have aborted under `set -u` |
| R07 naming all three Figure S2 panels `SuppFig2A_*` | Panel letters follow the published figure |
| R23 naming its output `Figure6A_*` | It is Figure S5A; renamed |
| Replicates framed as a circRNA-only inconsistency (4v4 vs 3v3) | n = 3 throughout; the historical notebooks are the outlier |
| `Rattus_norvegicus.mRatBN7.2.115.gtf` in the manifests | Corrected to `GRCr8.115.gtf` |
