# Changes from the originals

What this repository does differently from the scripts and notebooks received,
and why. **Ordered the way the pipeline runs**, so a stage can be reviewed
against its original without hunting.

Nothing here was changed silently. Every entry names the original file under
`archive/`, which holds every package received byte-for-byte and is never
executed.

### How to read an entry

| | |
|---|---|
| **Original** | the file under `archive/` this stage was built from |
| **Changed** | what differs, and why |
| **Manuscript** | where the v14 Methods disagree with the original script, and which one this repository follows |

Manuscript sources compared: `NE_NRVM_CDDIS_adapted_v14-SL.docx`,
`NE_NRVM_CDDIS_Supplementary_Information_v14.docx`,
`Final_Supplementary_Tables_S1-S8_v14_2.xlsx`, and figures 1-6 / S1-S5.
`bash_scripts.txt` below means
`archive/raw_scripts_original/bash_scripts/bash_scripts.txt`, the raw command
log; line numbers are to that file.

Anything still unresolved is in [`OPEN_ISSUES.md`](OPEN_ISSUES.md), not here.

---

# Part 1 — changes that apply to every stage

These were made once, in a shared place, rather than 27 times.

## 1.1 No path is baked in

The originals carried five different machine roots between them:
`C:/Users/sebau/...`, `C:/Users/sebau/Downloads`, `/mnt/data`,
`/home/vmaracaja/shared-projects/SUZ_sr_2025` and
`/home/vmaracaja/shared-projects/sebauz_sr`, plus `D:` and `E:` drive letters
in the circRNA notebooks. 17 scripts had one hard-coded;
`12_ceRNA_coherent_pairs_and_triplets_ALL_RECORDS_FINAL.R` had the same one
three times (`:8`, `:21`, `:28`).

Every R stage now takes `--in-dir` / `--out-dir` through optparse, and every
bash stage reads its directories from environment variables named in its
header. Inputs resolve by **basename** through `find_input()` in
`scripts/R/lib/io.R`, which searches `--in-dir` recursively and then the
working directory, and passes a full path through unchanged.

`scripts/checks/check_hardcoded_paths.sh` fails on any reintroduction.

## 1.2 Browser-download duplicates removed

Seven scripts listed `name(1).csv`-style filenames as **first-preference**
input candidates: `31_`, `32f_`, `35_`, `36b_`, `37c_`, `39e_` and
`24_FINAL_circRNA_lncRNA_miRNA_mRNA_cross_reference.R`. A stale copy left in a
download folder would have been picked over the real file.

All removed. `find_input()` rejects the `name(N).ext` pattern outright rather
than deprioritising it.

## 1.3 Every stage owns its output directory

The raw `12_*` script and all three `Figure8_*` variants wrote **identical
filenames** into one shared `Paper_ceRNA_candidates/` directory, so whichever
ran last won.

Each stage now writes to `results/R<NN>_<stage>/` and nothing writes into
another stage's directory. Each output directory also gets a
`00_INPUT_manifest.csv` recording which file was resolved for each input role,
with size and mtime.

## 1.4 One pinned Ensembl connection

Four unpinned biomaRt call sites existed: `useEnsembl(..., version = 114)` in
the metadata notebook, `useEnsembl(mirror = "uswest")`, and two bare
`useMart("ensembl")` calls in the circRNA annotation notebook.

All replaced by `connect_ensembl()` in `scripts/R/lib/ensembl.R`, pinned to the
archive host for the active release. `version=` is deliberately **not** passed
alongside an archive host — biomaRt's resolution between the two has changed
across package versions.

Two defects fixed in the same place:

- **`https://sep2025.archive.ensembl.org/` with a trailing slash** fails with
  `Error: Can't retrieve empty body`. `suz_normalize_host()` strips it.
- **`getBM(values = ...)` with no `filters = ...`**
  (`00_annotation_biomart.Rmd:47,129`) is silently ignored by biomaRt, so those
  chunks never retrieved what they appear to. `suz_getBM()` refuses the
  pattern instead of returning the wrong table.

`scripts/checks/check_biomart_pinning.sh` fails on any other call.

## 1.5 One CLI convention

`MASTER_results_report_builder_v3_8_1.R` was the only already-parameterised
script in the corpus. Its `commandArgs` pattern became the house style, lifted
to optparse: every R stage takes `--in-dir`, `--out-dir`, `--dpi` and `--help`.

While porting, 16 stages were found to rebind `out_dir` to a bare relative path
**after** the CLI had set it — every one of them would have written into the
current working directory instead of `--out-dir`. All now nest under
`.suz_out_root`.

## 1.6 Three shared libraries instead of copy-paste

| `scripts/R/lib/` | Replaces |
|---|---|
| `io.R` | per-script delimiter sniffing (`12_*:40-44` tab/comma, `24_FINAL_*` semicolons), the column resolver from `MASTER_results_report_builder:157`, and the Windows MAX_PATH guard from `32f_*:28-31` |
| `ensembl.R` | the four biomaRt call sites |
| `plotting.R` | per-script `ggsave` calls and per-script colour choices; one direction, RNA-class and temporal-phase palette so panels assembled side by side agree |

## 1.7 Figures are always a PDF/PNG pair

`save_figure()` writes `<stem>.pdf` (cairo where available, with a fallback)
and `<stem>_<dpi>dpi.png`, default 600 dpi. The originals mixed formats, dpi
and devices.

## 1.8 The conda environments

The original material records only one environment name, `IntaRNA_lncRNA`, in
`archive/raw_scripts_original/intarna_related/command_lines.txt`. Everything
else about the software stack had to be reconstructed from the commands and
version strings in the scripts, which is what
`metadata/software_tool_inventory.tsv` holds.

Five environments, matching the names the stages activate:

| File | Environment | Stages |
|---|---|---|
| `envs/environment_qc_align.yml` | `env_qc_align` | bash 01-10, 09b-09e |
| `envs/environment_mirge.yml` | `env_mirge` | bash 11 |
| `envs/environment_intarna.yml` | `env_intarna` | bash 12, 16 |
| `envs/environment_circrna.yml` | `env_circrna` | bash 13-15 |
| `envs/environment_r.yml` | `env_r` | all 27 R stages |

They are split rather than merged because the splits are forced: IntaRNA pulls
a Boost/ViennaRNA stack that conflicts with the aligners (and was separate in
the original analysis too), and the R/Bioconductor stack shares nothing with
the command-line tools.

**No version is pinned.** The solver takes the current release and the version
the Methods record sits in a trailing comment on each line, so nothing is lost.
`scripts/commands/create_envs.sh` writes `envs/*.lock.yml` from what actually
resolved, and that export — not the `.yml` — is the version record. This is a
deliberate consequence of re-running at Ensembl 115: see
[`OPEN_ISSUES.md`](OPEN_ISSUES.md) §15.2.

Four defects were fixed in the environment files themselves, none of which
would have survived a single `conda env create`:

| Defect | Fix |
|---|---|
| Every file declared `name: suz_<x>` while all 16 bash stages run `conda activate env_<x>`. **Nothing would have activated.** | Renamed to `env_<x>`, matching both the activation call and the error message `_activate_env.sh` prints. `create_envs.sh` now refuses to build a file whose `name:` does not match |
| `clear` was listed as a bioconda dependency. **It does not exist there** (404). CLEAR installs from source | Removed from `dependencies`, and its build requirements added instead — CIRCexplorer2, HISAT2, **StringTie** (which was missing), pysam, pybedtools. The `git clone` + `setup.py install` command is in the file header and `create_envs.sh` prints it |
| `mirge3=3.0` pinned a version that does not exist. The tool is branded miRge3.0 but the bioconda package is versioned independently, currently 0.1.4 | Unpinned, with the discrepancy explained in the header so nobody re-pins it |
| All five listed the `defaults` channel alongside conda-forge and bioconda, which breaks bioconda's required channel order and can fail outright on Anaconda's terms of service | `conda-forge`, `bioconda`, `nodefaults`. `create_envs.sh` also sets `channel_priority strict`, which bioconda requires and a YAML file cannot express |

Also: `openjdk` was dropped from the circRNA environment — it was there only
for CIRI-vis, which the same comment called "not part of the pipeline".

**Stages 09b-09e did not activate any environment**, although the chain's whole
reason for pinning `gawk` is that stage 09e needs GNU awk's three-argument
`match()`. All four now activate `env_qc_align`, after argument parsing so
`--help` still works without conda.

### The container

The README has always said the pipeline is Linux-first and that Windows should
go through WSL2 or Docker, but no `Dockerfile` existed. There is one now: a
single image carrying all five environments, with CLEAR built in.

Three choices in it are load-bearing:

- **The base is `condaforge/miniforge3`, not micromamba.**
  `scripts/bash/_activate_env.sh` calls `conda info --base` and then
  `conda activate` directly. micromamba provides neither, so every bash stage
  would fail to activate. miniforge3 also defaults to conda-forge and ships no
  `defaults` channel, which matches `envs/*.yml`.
- **`ENV BASH_ENV=/opt/conda/etc/profile.d/conda.sh`** rather than an entrypoint
  wrapper. `conda activate` does not work in a non-interactive shell unless
  `conda.sh` has been sourced, and `BASH_ENV` makes every `bash -c` source it.
  `/etc/bash.bashrc` covers the interactive case.
- **`--platform linux/amd64`.** bioconda publishes linux-64 only; without it the
  build silently fails to solve on Apple silicon.

The environments are built from `envs/` *before* the rest of the repository is
copied, so editing a stage script does not invalidate the expensive
environment layers. The build then runs the repository's own checks and
`create_envs.sh --verify-only`, so a broken image fails instead of shipping,
and finally writes `envs/*.lock.yml` inside the image — the container carries
its own version manifest.

`archive/` is deliberately **not** in `.dockerignore`: it is 1.5 MB, and 90 of
the 104 rows in `metadata/frozen_inputs_checksums.tsv` point into it, so
excluding it would make `make_checksums.sh --verify` fail inside the container.

## 1.9 Line endings

`.gitattributes` forces LF on every script. The pipeline is Linux-first and
CRLF breaks shebang lines.

---

# Part 2 — the bash stages, in run order

## Stage 00 — `00_download_references.sh`

**Original:** `archive/starter_v1_original/metadata/reference_manifest.tsv` and
`commands/example_future_pipeline_call.sh`.

**Changed — the annotation filename.** The starter package, the circRNA
addendum and the manuscript Methods all name
`Rattus_norvegicus.mRatBN7.2.115.gtf`. **That file has never existed.** The rat
assembly switched from mRatBN7.2 to GRCr8 at release 114, so release 115 ships
`Rattus_norvegicus.GRCr8.115.gtf`. Wherever the impossible name appears, a
release number and an assembly name have been combined from two different
sources.

| Release | Rat GTF | Assembly | biomaRt archive |
|---|---|---|---|
| 112, 113 | `mRatBN7.2.<rel>.gtf.gz` | mRatBN7.2 (rn7) | — |
| **114** | `GRCr8.114.gtf.gz` | **GRCr8 (rn8)** | `may2025.archive.ensembl.org` |
| **115** | `GRCr8.115.gtf.gz` | **GRCr8 (rn8)** | `sep2025.archive.ensembl.org` |

**Changed — primary assembly, not toplevel.** The stage now builds
`dna.primary_assembly`, preferring a single published archive and otherwise
downloading the per-chromosome parts and concatenating them. It **refuses to
fall back to toplevel silently.**

The manuscript Methods were right about this and an earlier version of this
repository was wrong: it had been downloading `dna.toplevel.fa` on the
incorrect claim that Ensembl publishes no primary assembly for GRCr8. Ensembl
publishes `Rattus_norvegicus.GRCr8.dna.primary_assembly.*`, split per
chromosome, and its own FASTA README settles what that means — toplevel
"includes ... N padded haplotype/patch regions", while the primary assembly
excludes them and "is best used for performing sequence similarity searches
where patch and haplotype sequences would confuse analysis". Only when no
primary-assembly file is published are the two equivalent.

The raw command log shows this was understood at the time: `:46` marks the
toplevel index **`(NO USAR)`** — do not use — and `:49-50` points at the
primary-assembly index. The circRNA branch agrees:
`rat_primary_assembly_orderr8.fa`.

`PRIMARY_PARTS` defaults to `1..20 X Y MT`. **`nonchromosomal` is excluded**
(B. Borges, 8 Oct 2026).

**Changed — one release switch.** `metadata/reference_release.tsv` drives the
GTF path, the FASTA URLs and the biomaRt archive together. `ENSEMBL_RELEASE`
overrides for one run. A `url_exists()` probe now fails with the URL it tried
rather than writing a zero-byte file.

**Manuscript:** genome FASTA correct; GTF filename impossible and the release
number contested — see [`OPEN_ISSUES.md`](OPEN_ISSUES.md) §3.1.

## Stages 01-05 — QC, trimming, alignment, BAM

**Original:** carried over from the starter package's stage scripts, which were
already sound.

**Changed:** headers normalised to one house style (purpose, settings, `Usage`,
`Env`, `Next`). No parameter changed.

**Manuscript:** agrees throughout — fastp
`-f 10 -F 10 -3 --cut_tail_mean_quality 28 -q 30` (`:66`), HISAT2
`--dta-cufflinks -S` (`:78-90`). The Methods omit fastp's `-3` and the thread
count; both are carried.

## Stage 06 — `06_stringtie_assemble_merge.sh`

**Changed:** `STRINGTIE_EXTRA` added, defaulting to **empty** — plain
reference-guided assembly, which is what the Methods describe.

**The GRCr8 StringTie invocation was never recorded.** Every StringTie command
in `bash_scripts.txt` is rn7-era (`mRatBN7.2.ensembl.gtf` or `.112.gtf`). Among
those, `--conservative` appears in five (`:107`, `:109`, `:152`, `:154`,
`:156`) and two early single-sample tests add `-f 0.05 -m 200 -t -s 3.0`
(`:107`, `:109`). Which, if any, carried into the final run is unknown, so the
flags are a variable rather than a default:

```bash
STRINGTIE_EXTRA="--conservative" bash scripts/bash/06_stringtie_assemble_merge.sh
```

**Manuscript:** describes guided assembly with no extra parameters, which is
the default here. See [`OPEN_ISSUES.md`](OPEN_ISSUES.md) §3.4.

## Stages 07-08 — gffcompare, gffread

**Changed:** headers only. gffcompare `-r <reference>` then the `.tmap`;
gffread extracts all transcripts plus the novel intergenic subset
(`class_code "u"`). Both match `:149`, `:182`, `:185`, `:188` and the Methods.

## Stage 09 — `09_feelnc_classify.sh`

**Changed — the classifier window is 1 kb, not 200 bp.**

The Methods say "a 200 bp to 100 kb window (`--window=200 --maxwindow=100000`)".
The final 6+24 h GRCr8 run, `:300`, is:

```
FEELnc_classifier.pl ... --window=1000 --maxwindow=100000
  > lncRNA_classes_CNE624_rn8_8_annotated_1kb_100kb.txt
```

That output filename is exactly what the downstream stages consume.
`--window=200` appears only in six earlier exploratory runs (`:199`, `:207`,
`:209`, `:219`, `:227`), every one of which used maxwindow 10 000 or 50 000 and
produced a differently named file that nothing reads. The Methods sentence
pairs the window from the discarded runs with the maxwindow from the final one.

This repository follows the data: **1 kb to 100 kb**, with the evidence quoted
in the script header. The Methods text needs correcting.

**Changed — codpot parameters carried forward.** `bash_scripts.txt` has the
heading `M2` with **nothing under it** (`:295-297`): Module 1 and Module 3 of
the final run were recorded, Module 2 was not. The stage uses the parameters
that were identical across all four earlier recorded codpot runs —
`--mode=intergenic -b transcript_biotype=protein_coding -b transcript_status=KNOWN`
— plus `--monoex=1`, which the Methods confirm.

**Manuscript:** three modules and `--monoex=1` agree; the window does not.

## Stage 09b — `09b_feelnc_split_exon_gtfs.sh` — RECONSTRUCTED

**Original: none.** The metadata chain begins by reading an `exons/` directory
whose generating command is in neither the metadata folder nor
`bash_scripts.txt`.

It was identified from the data, and the evidence is strong:

- FEELnc's 42,040 classified lncRNA transcripts are **100 % contained** in
  `exons/*.lncRNA.gtf` — the set difference in that direction is exactly 0;
- `.lncRNA.gtf` and `.otherRNAs.gtf` are **disjoint**;
- all 46,030 transcript IDs are gffcompare `.tmap` `qry_id`s;
- the subset spans all 13 class codes, so it is not class-code driven;
- the files are exon lines only, with the minimal attribute set.

That is the FEELnc_codpot signature: `.annotated.gtf` is the candidate set,
`.lncRNA.gtf` the coding-potential lncRNA calls, `.otherRNAs.gtf` the
complement. Two residuals confirm it — 42,127 − 42,040 = **87** lncRNA
transcripts the classifier drops as orphans with no partner in the
1 kb/100 kb window, and 46,030 − 42,127 − 3,822 = **81** in codpot's `noORF`
bucket.

`scripts/checks/check_feelnc_split_counts.sh` reproduces 46,030 / 42,127 /
3,822 with both invariants, and **currently passes exactly** against the
historical directory. That check is the only evidence the reconstruction is
right; the command itself still needs confirming.

## Stages 09c-09e — the transcript-metadata chain

**Originals:** `archive/transcriptome_metadata_original/` —
`extract_transcript_ids.sh`, `subset_merged_gtfs.sh`,
`filter_out_lncRNAs.sh`.

These three shell scripts run **before** the metadata notebook that became R01,
and R01 is unreproducible without them. Previous versions of this repository
omitted the chain entirely.

**Changed — the hard-coded `case` block.** `subset_merged_gtfs.sh` dispatched on
a 12-arm `case` with no arm for `CNE_24_STAR.annotated`, only the misspelled
`CNE24_STAR.annotated` that happened to match one odd filename. A
correctly-named file fell through to `*)` and was skipped with a warning only.
Stage 09d derives the dataset list from the sample sheet instead.

**Changed — one working directory.** The three originals carry two different
server roots: `SUZ_sr_2025` in the first two, `sebauz_sr` in the third. That
looked like a sign the chain had been run in pieces against different data, but
it was not — the directory was renamed by accident partway through, and it is
one data generation (B. Borges, 15 Oct 2026). All three stages now take one
`--work-dir`, which is what the originals meant.

**Changed — intermediates go under `results/`.** The `<id>_transcript_ids.txt`
files and `transcript_extraction_*.log` were written beside their inputs. They
now land in the stage's own output directory.

**Changed — gawk is now required, loudly.** `filter_out_lncRNAs.sh:26` uses GNU
awk's three-argument `match()`. Under mawk or BSD awk the capture group is
never populated, **every line passes the filter**, and the output is an
unfiltered GTF that looks correct. Stage 09e checks for gawk and aborts;
`envs/environment_qc_align.yml` pins it as a hard dependency, not a nicety.

**Changed — HISAT2 only.** The originals covered both the HISAT2 and STAR
branches. Confirmed with S. Urquiza that only HISAT2 feeds the manuscript, so
the STAR half is dropped.

## Stage 10 — `10_salmon_index_quant.sh`

**Changed — a decoy-aware gentrome index.** The Methods say "a genome-aware
decoy/gentrome index", and they are right. `:232-235`:

```
salmon index -t gentrome_CNE624_8_H.fa -i transcripts_CNE624_H_index \
             --decoys decoys.txt -p 6 -k 31
```

This repository had been building a transcriptome-only index, which gives
different quantifications. Stage 10 now builds a gentrome (transcripts +
genome) with a decoy list, same `--decoys` and `-k 31`. `BUILD_DECOY=0`
restores the old behaviour and warns.

Quantification at `:286-288` uses `-l A -p 8 --validateMappings`, and Salmon
v1.10.3 is named at `:280`; both agree with the Methods.

**Changed:** `REPO_ROOT` was used but never defined — the stage would have
aborted under `set -u`.

## Stage 11 — `11_mirge3_run.sh`

**Original:** `archive/raw_scripts_original/miRge_scripts/miRge_3.txt`.

**Changed — the flags the Methods name do not exist.** The Methods say "the rat
miRBase library (`-lib` and `-rat`)". There is no `-rat` flag. The recorded
invocation is:

```
-lib <miRge3_Lib> -on rat -db mirbase
-a AGATCGGAAGAGCACACGTCT          3' adapter
-g GTTCAGAGTTCTACAGTCCGACGATC     5' adapter
-dex -mdt DESmetadata_24h.csv
```

`-on rat` is the organism flag and `-db mirbase` is what makes it miRBase.
Confirmed by S. Urquiza, 2 Jul 2026: the Methods text will be corrected.

**Changed:** parameterised. Only the 24 h command survives — no 6 h command was
ever recorded, although `miR.Counts_C6_NE6_3v3_new.csv` is a required input.

## Stage 12 — `12_intarna_run.sh`

**Original:**
`archive/final_package_v3_original/03_INTARNA_COMMANDS/IntaRNA_commands_as_used.sh`
and `archive/raw_scripts_original/intarna_related/command_lines.txt`, the only
record of `conda activate IntaRNA_lncRNA`.

**Changed:** the four lncRNA–miRNA commands are parameterised, keeping
`--mode=M --outMode=C` and, critically, **the output filenames**, which must
keep matching the `Interact_*.csv` glob that R21 reads. Verified against
`IntaRNA_manifest.tsv`.

All eight runs (four here, four in stage 16) deliberately pair UP ncRNAs
against DOWN miRNAs and vice versa, because a sponging relationship predicts
opposite directions. The Methods' MFE ≤ −20 kcal/mol is applied downstream in
R21.

## Stages 13-16 — the circRNA branch

The 31 Jul 2026 addendum reintroduced circRNA as an independent module. The
previous repository still described it as archived and "not included in the
final 7-figure manuscript", which is no longer true.

### Stage 13 — `13_circrna_circexplorer.sh`

**Originals:** `CIRCexp3_run_all_2023_SebaNECardio.sh`, `_2.sh`, `_3.sh`.

**Changed:** three near-identical scripts consolidated into one sample-sheet
loop over H1-H4, SQ1-SQ8 and SQ13-SQ16. `set -euo pipefail` added. **Each
sample gets its own scratch directory** — the originals shared one `*.fq` glob
and were unsafe to run concurrently.

The commands themselves are kept verbatim: `CIRCexplorer2 parse`,
`CIRCexplorer2 annotate`, and CLEAR's `run3.py` and `circ_quant.py`.

**Manuscript: the Methods are right, and an earlier version of this document
was wrong to question them.** It said the Methods should "name what ran",
implying CIRCexplorer3 was a misnomer for CIRCexplorer2. It is not.

**CLEAR and CIRCexplorer3 are the same tool.** The YangLab repository is
titled "CLEAR/CIRCexplorer3". There is no separate CIRCexplorer3 binary to
invoke: CLEAR declares CIRCexplorer2 >= 2.3.6 as a dependency and calls its
`parse` and `annotate` subcommands itself, then adds `circ_quant` for FPBcirc,
FPBlinear and CIRCscore. The original script makes this plain — it runs
`.../CIRCexp3/CLEAR/src/run3.py` and `.../CIRCexp3/CLEAR/src/circ_quant.py`
around the two CIRCexplorer2 calls.

So seeing `CIRCexplorer2` in the commands is what using CIRCexplorer3 looks
like, and the pipeline is already doing it.

Two real defects remain in that Methods sentence:

- **"v1.0.0" is invented.** CLEAR publishes no version number at all. Replace
  it with the citation and the installed commit, which the image records at
  `/opt/CLEAR/INSTALLED_COMMIT`.
- **The CLEAR paper is not cited.** Ma XK, Wang MR, Liu CX, Dong R, Carmichael
  GG, Chen LL, Yang L. *A CLEAR pipeline for direct comparison of circular and
  linear RNA expression.* doi:10.1101/668657.

One note for the Methods, not a defect: `circ_quant` is run, so FPBcirc,
FPBlinear and CIRCscore are computed, but differential expression uses the raw
`readNumber` back-splice counts from `annotate`. A reader told the pipeline is
CLEAR will expect CIRCscore to appear somewhere; it does not.

### Stage 14 — `14_circrna_ciriquant.sh`

**Originals:** `CIRIquant_all_run_2023_Seba_NE_Cardio_1.sh`, `_2.sh`, `_3.sh`.

**Changed:** same consolidation. The `rn8.yml` config is not in any package, so
the stage **fails loudly with a pointer** rather than guessing a configuration.

### Stage 15 — `15_circrna_ciri_as.sh`

**Original:** `CIRI_2_run_all_seba.sh`.

**Changed:** parameterised BWA-MEM `-T 19 -t 20` → CIRI2 → CIRI-AS. The sample
loop is **deliberately not widened**: the original is `for m in {13..16}`, so
SQ13-SQ16 only, and no CIRI2/CIRI-AS loop exists for H1-H4 or SQ1-SQ8 anywhere.
The stage defaults to those four and refuses `--all` without `--yes-widen`.

**Not carried:** two defects in the sibling `vis_run_ciri2_seba.sh` — its sample
list reads `Q13 Q14 Q15 Q16` instead of `SQ13…SQ16`, and it `rm`s the PDFs it
has just generated. Both are preserved in `archive/` and not reproduced.

### Stage 16 — `16_circrna_intarna.sh` — RECONSTRUCTED

**Original:** the four circRNA–miRNA commands in
`IntaRNA_commands_as_used.sh`, with no runner script.

**Changed:** rebuilt as a stage, same `--mode=M --outMode=C` and the same output
names. It depends on CIRI-AS sequences, so stage 15's coverage question gates
it.

---

# Part 3 — the R stages, in run order

## R01 — `01_transcriptome_metadata_annotation.R`

**Original:** the transcriptome-metadata Quarto notebook
(`2_nrvm_transcriptome_metadata_handling.qmd`, B. Borges). Per the decision
taken for this rebuild, no `.qmd` enters the repository; the notebook originals
stay outside it and `archive/ARCHIVE_NOTES.md` records where. The three shell
scripts of the same chain **do** come in, as stages 09c-09e.

**Changed:** ported to a clean `Rscript` with the standard CLI, HISAT2 config
only, carrying steps 1-8 of the notebook — tmap read, `getBM`, compatibility
filter, exon lengths from the GTF, `_meta_ann.csv`, `_meta_novel_lncRNA.csv`,
FEELnc validation, `_final_annotated.gtf`.

**Changed:** biomaRt pinned. The notebook called
`useEnsembl(..., version = 114)`; the stage uses `connect_ensembl()` and adds
`--dry-run-biomart`, which reports the host and release without running the
analysis.

**Changed:** depends on bash 09b-09e, which previous versions of the repository
did not contain at all.

**Noted, not changed:** every row of the shipped classifier output
`lncRNA_classes_CNE_6_24_HISAT2_1kb_100kb.txt` already has `isBest == 1`
(42,040 of 42,040), so the `isBest` filter in R08/R14 is a **no-op on this
input**. The filter is kept, but it is not evidence that filtering occurred and
should not be described as a filtering step.

**The shipped tables are frozen rather than regenerated.** `data/frozen/` holds
`hisat2_6_24_meta_ann.csv.gz` and `hisat2_6_24_meta_novel_lncRNA.csv.gz`,
produced by this notebook at Ensembl 114, because Figure 2B was built from
them and they cannot be regenerated byte-identically — the committed notebook
is not quite the code that produced them (`.Rhistory` shows
`run_full_analysis()` called with a fourth `exon_attributes` argument the
committed version does not accept). `data/frozen/README.md` has the full
reasoning, the checksums and the column lists. Pass `--meta-ann` /
`--meta-novel` to use freshly generated output instead.

**Figure 2B reproduces exactly from them** — 21,981 protein-coding / 16,149
annotated lncRNA / 840 novel lncRNA / 1,778 snoRNA / 1,517 snRNA / 204 miRNA,
every number in the published panel, with no tuning. That is the strongest
single piece of evidence that the frozen tables are the right inputs and that
the reconstructed chain 09b-09e is faithful.

## R02 — `02_deseq2_mRNA_lncRNA.R`

**Original:**
`archive/raw_scripts_original/R_scripts/Integrated_Bulk_v3_6_24h.R` (2,564
lines).

**Changed:** CLI and path handling per Part 1. The multifactorial design
(`~ time + treatment + time:treatment`), the thresholds (padj < 0.05,
|log2FC| ≥ 1) and the correlation settings (**Pearson**, |rho| ≥ 0.8, p < 0.05,
on VST counts) are carried unchanged and all agree with the Methods.

**Noted:** `:30-33` reads `Ctrl6_rep1, Ctrl6_rep2, Ctrl6_rep4, …` — **three
replicates per group, replicate 3 excluded.** That holds across the whole
corpus: the miRNA matrices are `*_3v3_new.csv` and the final circRNA volcano
script drops `_3`. The manuscript legends say N = 3, so the manuscript is
internally consistent; the 4v4 historical circRNA notebooks are the outlier.
`metadata/longRNA_samples.tsv` marks the four replicate-3 rows
`include_final_analysis = no_excluded_from_final_DE`.

The sample sheets describe 16 long-RNA libraries; 12 were analysed.

**Changed — one replicate numbering across every sheet.** Two incompatible
schemes were in `metadata/` at once, and `fastq_manifest.tsv` disagreed with
itself: its long-RNA rows used one, its small-RNA rows the other, which is also
the one in `longRNA_samples.tsv` and Seba's starter sheet. Because R04 keys its
sample labels on `condition + time + replicate` and R04b filters with
`--replicates`, "drop rep3" resolved to two disjoint sets of libraries.

The long-RNA block of `fastq_manifest.tsv` was confirmed authoritative
(B. Borges, 19 Oct 2026) and every other sheet renumbered to it. It is the only
scheme that pairs each control with its norepinephrine partner by consecutive
sample ID at both time points — 6 h: H1/H2, SQ1/SQ2, SQ5/SQ6, SQ13/SQ14;
24 h: SQ15/SQ16, H3/H4, SQ3/SQ4, SQ7/SQ8. The other breaks that pattern at
24 h.

So the excluded libraries are **`SQ5_ln`, `SQ6_ln`, `SQ3_ln` and `SQ4_ln`**,
and the Methods can name them. Full evidence in
[`OPEN_ISSUES.md`](OPEN_ISSUES.md) §8.1. `salmon_quant_alias` remains
`TO_CONFIRM`: nothing records whether the person who made the quant
directories used this numbering, so the mapping is a strong inference rather
than a document.

## R03 — `03_mirna_deseq2.R`

**Original:** `archive/raw_scripts_original/R_scripts/miRNAs_6_24.R` (1,560
lines).

**Changed:** CLI and path handling only. Input is `miR.Counts_C6_NE6_3v3_new.csv`
and `miR.Counts_C24_NE24_3v3_new.csv`.

## R04 and R04b — the circRNA reconstruction

**Original: none.** Seven final scripts read four file families that **nothing
in any package produces**:

| File family | Read by |
|---|---|
| `raw_count_matrix_all_{6,24}.txt` | R06, R20 |
| `0{6,7}_circRNA_sequences_de{6,24}_DE_complete_table.txt` | R13, R20, R21, R22 |
| `circRNA_{UP,DOWN}_{6,24}h_all_records.fa` | bash 16 |
| `12_STRICT_DE_coherent_circRNA_triplets.csv` | R23, R24 |

The last one turned out **not** to be an orphan: it is output `12_` of the
cross-reference stage. Reordering so R22 runs before R23/R24 removed that gap
outright. The other three are rebuilt here from
`archive/final_package_v3_original/02_CIRCRNA_UPSTREAM_RAW_PROVENANCE/03_tool_comparation.Rmd`.

**Changed — split into two stages**, following the v3 merge plan's
recommendation rather than collapsing them, and the split is load-bearing:

- **R04 is deterministic** — parse both callers, reconcile coordinates, merge
  loci, emit an **unfiltered** BSJ count matrix and a locus annotation table.
- **R04b holds every contested parameter** as a flag: `--bsj-min` (3),
  `--replicates` (1,2,4), `--use-combat` (off), `--alpha` (0.05), `--lfc-cut`
  (0.5).

Because the BSJ rule, the replicate set and ComBat_seq are all still open, DE
has to be re-runnable under either rule without re-parsing every caller output.

**Carried verbatim from the notebook:**

- the **asymmetric `+1` start-coordinate correction applied to CIRCexplorer
  only, not CIRIquant** (`03_tool_comparation.Rmd:16` vs `:23`) — a deliberate
  0-based/1-based reconciliation, not a bug;
- merge key `c("chr","Start","End")` and the `software` tags
  `circexp3_and_ciriquant` / `circexp3` / `ciriquant`;
- the recorded locus counts — intersection 46,302, CIRCexplorer-only 3,419,
  CIRIquant-only 23,739 — asserted on every run, **warning not erroring** on
  mismatch, so a divergent rerun says so immediately.

**Changed:** R04b emits **raw** counts, not the normalised matrix the notebook
writes, because the consumers want `raw_count_matrix_all_*`.

**Changed:** `circRNA_sequences_de{6,24}_2.fa` are written by the notebook with
`write.table(sep="\t", col.names=TRUE)` and read back with
`read.table(header=TRUE)` — they are **TSVs with a `.fa` extension**. The
behaviour is kept; the outputs are named `.tsv`.

**Manuscript — two contradictions this stage exposes.**

- **BSJ threshold.** The Methods say "only circRNAs with **at least two** BSJ
  counts were carried forward". The code is `readNumber > 2`
  (`01_circexp3.Rmd:91`, `02_ciriquant.Rmd:85`) — i.e. **at least three**.
  These give different inventories, so Figure 2B's 8,736 loci and every
  circRNA count in Figures 3 and S3 move with it. `--bsj-min` defaults to 3,
  matching the code rather than the prose, and logs the choice on every run.
  **Unresolved and blocking.**
- **ComBat_seq.** The Methods do not mention it. The historical notebooks apply
  it before DE (`01_circexp3.Rmd:132`, `:219`); the final volcano script does
  not. `--use-combat` is off by default. Note that `batch` maps 1:1 to
  replicate index (`01_circexp3.Rmd:124`, `:212`), so excluding replicate 3
  removes **Batch3 entirely** — the exclusion is batch-aligned, not random.
  **Unresolved.**

## R05 — `05_tf_kegg_integration.R`

**Original:**
`archive/final_package_v3_original/01_FINAL_FIGURE_AND_ANALYSIS_SCRIPTS/script_final_depurado_TF_KEGG.R`
(1,300 lines, three stages in one file).

**Changed:** CLI and path handling; the three stages take explicit inputs
instead of relying on objects left in the session. The GRN rules (TFLink,
expressed > 0 in at least one replicate, mean normalised count > 10) are
carried unchanged and agree with the Methods.

Needs `TFLink_All.tsv` and three Cytoscape GRN edge exports, none of which were
supplied.

## R06 — `06_transcriptome_overview.R`

**Original:** `.../26_transcriptome_overview_with_circRNA.R`.

**Changed:** reads the frozen metadata tables through `read_frozen()` by
default. Produces **Figure 2B**, which reproduces the published panel exactly
(see R01).

## R07 — `07_lncRNA_characterization_gtf.R`

**Original:** `.../15_lncRNA_characterization_from_GTF.R`.

**Changed — panel letters.** The script named all three outputs `SuppFig2A_*`.
The published Figure S2 is **A** exon-count distribution, **B** isoform
complexity, **C** transcript length, so they are now `FigureS2A_`,
`FigureS2B_`, `FigureS2C_`.

**Noted:** the script infers novelty from the `^MSTRG` / `^ENSRNOG` prefixes,
and its own summary carries a caveat worth reusing in the manuscript — because
this is a transcriptome reconstruction framework, novel MSTRG entries are best
described conservatively as novel reconstructed lncRNA loci rather than
definitively novel lncRNA genes, unless orthogonal evidence supports
locus-level novelty.

## R08 — `08_feelnc_panels.R`

**Original:** `archive/raw_scripts_original/R_scripts/14_FEELnc_supplement_panels.R`,
with the equivalent block at `script_miRNA_mRNA_6_24.R:438+`.

**Changed:** CLI and path handling. The `isBest == 1` filter is kept and
documented as a no-op on this input (see R01).

## R09 — `09_de_counts_heatmaps_gba.R`

**Original:**
`.../35_FINAL_missing_figures_miRNA_heatmap_DE_GBA.R`.

**Changed — the curated KEGG allowlist moved out of the code.** It was an
inline `tribble` at `35_:30+`, silently duplicated against
`01C3_selected_KEGG_allowlist_used.csv` written by `32f_`. It is now
`metadata/curated_kegg_allowlist.tsv`, and both this stage and R24 read the
same copy.

Produces **Figures 3A, 3C miRNA, 5A and S4A**.

## R10 — `10_figure3B_temporal_programs.R`

**Original:** `.../36b_Figure3_panelB_DOWN_bottom_UP_top.R`.

**Changed:** the `/mnt/data` input path at `36b_:30` removed. Stacked-bar
ordering (DOWN at the base, UP on top, white in-bar counts) is kept — both
`36b_` and `39e_` insist on it, so it lives in `lib/plotting.R` now.

**Manuscript — the published panel is stale.** For lncRNA only, Figure 3B
disagrees with both the Results text and Figure S3C, which agree with each
other:

| lncRNA program | Results text | Figure S3C (derived) | Figure 3B (published) |
|---|---|---|---|
| Early | 148 (62 up, 86 down) | 148 (62, 86) | **145 (61, 84)** |
| Sustained | 258 (86 up, 172 down) | 258 (86, 172) | **275 (95, 180)** |
| Late | 336 (109 up, 227 down) | 336 (109, 227) | **328 (104, 224)** |

Summing the annotated and novel matrices of Figure S3C gives
Early = (47+15) up, (70+16) down; Sustained = (66+20), (128+44);
Late = (92+17), (163+64) — every number matching the Results text. mRNA, miRNA
and circRNA reconcile across Figure 3A, 3B, S3C and the text. **It is the
figure that is wrong**, and regenerating it with this stage fixes it.

## R11, R12 — Figure 3C mRNA and lncRNA heatmaps

**Originals:** `.../05_Figure3_mRNA_programs.R` and
`.../37c_lncRNA_temporal_heatmaps_annotated_novel_MSTRG_FIXED.R`.

**Changed:** CLI and path handling. `37c_`'s reclassification of MSTRG-prefixed
entries as novel is kept.

## R13 — `13_figure3C_circRNA_heatmap.R` — NEW

**Original: none.** No standalone script for this panel existed in any package;
it was flagged `FIG3C_CIRC` as open in the v3 audit.

Written to follow its sibling panels: classifies loci as Early, Sustained, Late
or Reversal from the two state-specific DE tables, draws the heatmap, and
writes the per-phase locus lists R21 reads. **Its output has never been
compared against the manuscript figure.**

## R14 — `14_mirna_mRNA_integration.R`

**Original:** `.../script_miRNA_mRNA_6_24.R`.

**Changed:** CLI and path handling; emits the miRNA FASTAs that bash 12 and 16
consume. Needs the 12 `miRWalk_miRNA_Targets_*.csv` files and miRBase
`mature.fa`, neither supplied. The miRWalk settings (score ≥ 0.95,
3'UTR/5'UTR/CDS, miRDB support) agree with the Methods.

Produces **Figures 5B-C**.

## R15 — `15_lncRNA_hubs.R`

**Original:** `.../06_Figure4_refine_heatmap_hubs.R`.

**Changed — this is no longer a figure.** The v3 audit map listed a top
lncRNA-hubs panel as Figure S4C. The published Figure S4 has **only one
panel**, the curated KEGG representation (R09). The hub rankings live in
Supplementary Table S5D instead, so this stage produces tables.

## R16 — `16_volcanoes_lncRNA.R`

**Original:** `.../07_Figure4_clean_volcanoes.R`.

**Changed:** CLI and path handling only.

**Manuscript:** the Methods require lncRNAs to pass adjusted p < 0.05 **and
|log2FC| ≥ 0.5**. The script sets `absfc_cutoff <- 0` — significance on padj
alone. This stage is **left as ported**, because it produces a published panel
and changing the threshold would change the figure. Figure S3's legend says the
dashed lines mark "the thresholds used for classification", so the panel and
the Methods should be made to agree one way or the other. **Unresolved.**

## R17 — `17_volcanoes_mRNA_miRNA.R` — NEW

**Original: none** as a standalone script. The plotting code was embedded in
`Integrated_Bulk_v3_6_24h.R` and `miRNAs_6_24.R` and was never isolated;
flagged `S3_MRNA_MIRNA` as open in the v3 audit.

Applies the class-specific Methods thresholds: `--lfc-cut-mrna 1`,
`--lfc-cut-mirna 0.5`, both with padj < 0.05. **Its output has never been
compared against the manuscript figure.**

## R18 — `18_directional_overlap.R`

**Original:**
`.../39e_lncRNA_directional_overlap_barplot_UP_top_DOWN_bottom_FINAL.R`.

**Changed:** extended from lncRNA only to mRNA, lncRNA and miRNA, because the
published Figure S3C has all of them. The circRNA panel of S3C comes from R22.

## R19 — `19_extract_lncRNA_for_intarna.R`

**Original:** `archive/raw_scripts_original/R_scripts/09_extract_lncRNA_candidates_for_IntaRNA.R`.

**Changed:** CLI and path handling. Writes the lncRNA target FASTAs bash 12
reads, so it must run before it.

## R20 — `20_circRNA_volcanoes.R`

**Original:** `.../25_recreate_circRNA_volcanoes_with_all_NS.R`.

**Changed:** CLI and path handling. Thresholds kept: nominal p < 0.05 and
|log2FC| ≥ 0.5 (`:26-27`), which agree with the Methods.

**Noted:** this is a **nominal** p value, not an FDR. circRNA differential
expression is exploratory throughout — the phase-aware script says so
explicitly — and Figure 3C circRNA and Figure S3A-B circRNA should be described
that way.

## R21 — `21_ceRNA_coherent_pairs_triplets.R`

**Original:**
`.../12_ceRNA_coherent_pairs_and_triplets_ALL_RECORDS_FINAL.R` (982 lines).

**Changed:** the same hard-coded `base_dir` appeared three times (`:8`, `:21`,
`:28`); all removed. Reads the eight `Interact_*.csv` files from bash 12 and
16.

**Manuscript — "no additional interaction tiers were used".** The Methods say
"Only interactions with MFE ≤ −83.7 kJ/mol (−20 kcal/mol) were retained; no
additional interaction tiers were used." The script defines three —
`energy_high -20`, `energy_moderate -15`, `energy_exploratory -10` — and sets
`minimum_pair_energy <- energy_exploratory`, so intermediate outputs retain
pairs down to −10.

The **final** networks do use −20
(`HIGH_CONFIDENCE_ceRNA_triplets_E_le_minus20.csv`, and the MFE20 naming
through R24/R25), so the figures match the Methods; the sentence is just not
literally true of the script. The tiers are retained because they produce audit
output. Either soften the sentence or state that the tiers exist only as
intermediate audit output.

## R22 — `22_circRNA_cross_reference.R`

**Original:**
`archive/circrna_addendum_v1_original/scripts/24_FINAL_circRNA_lncRNA_miRNA_mRNA_cross_reference.R`
(894 lines).

**Changed — moved earlier in the order.** Its output `12_` is
`12_STRICT_DE_coherent_circRNA_triplets.csv`, which R23 and R24 read. In the
original ordering those two ran first and the file looked like an orphan
external input. Running R22 before them removes the gap.

**Changed:** maps mRNA IDs to symbols through the frozen `meta_ann` table.

Produces the **Figure S3C circRNA** panel and the cross-ncRNA modules.

## R23 — `23_ceRNA_phase_aware_plots.R`

**Original:** `.../31_phase_aware_ceRNA_coherence_plots.R`.

**Changed — this is Figure S5A, not Figure 6A.** The v3 audit map had Figure 6
and Figure S5 transposed. The published figures are:

| | Published v14 |
|---|---|
| **Figure 6A** | integrated **UP-DOWN-UP** network (Sustained + Late), 197 triplets |
| **Figure 6B** | integrated **DOWN-UP-DOWN** network (Early + Sustained + Late), 2,396 triplets |
| **Figure S5A** | phase-aware median log2FC trajectories with IQR |
| **Figure S5B-F** | the five phase-specific networks |

The output was named `Figure6A_phase_timepoint_ceRNA_coherence_RECOMMENDED.pdf`
and is now `FigureS5A_phase_timepoint_ceRNA_coherence`.
`metadata/figure_script_map.tsv` was rewritten from the v14 legends rather than
the audit map, which predates the final figures — **13 panels were remapped.**

## R24 — `24_ceRNA_cytoscape_export.R`

**Original:** `.../32f_export_ceRNA_cytoscape_TF11_curatedKEGG_MFE20.R`
(1,163 lines).

**Changed — an explicit `ALL_MFE20` pair.** R25 looks for
`ALL_MFE20_{nodes,edges}.csv`, but this script wrote `<module>_{nodes,edges}.csv`
— the two only ever met if a module happened to be named `ALL`. The stage now
emits an explicit `ALL_MFE20_nodes.csv` / `ALL_MFE20_edges.csv` pair alongside
the per-module tables.

**Changed:** reads the curated KEGG allowlist from `metadata/` rather than
rebuilding it. The Windows MAX_PATH guard from `32f_:28-31` moved into
`lib/io.R`, where every stage gets it.

Produces the node/edge tables behind **Figures 6A-B and S5B-F**.

## R25 — `25_cytoscape_key_validation.R`

**Original:** `.../33_validate_fix_cytoscape_import_keys.R`.

**Changed:** takes the `ALL_MFE20` pair as a CLI argument instead of assuming
its filename. Writes `*_CYTOSCAPE_SAFE.tsv`,
`ALL_MFE20_edges_MINIMAL.tsv` and an ID-match report.

## R26 — `26_results_report.R`

**Original:** `archive/raw_scripts_original/R_scripts/MASTER_results_report_builder_v3_8_1.R`
(1,512 lines) — the only script in the corpus that already took command-line
arguments. Its `commandArgs` pattern and its column resolver became the house
style for all 27 stages.

**Changed:** CLI lifted to optparse; the column resolver generalised into
`choose_col()` in `lib/io.R`.

---

# Part 4 — the full Methods-versus-scripts list

Everything from Parts 2 and 3 in one table, so the Methods can be revised in
one pass. "Follows" is what this repository does.

| # | Stage | Methods say | Scripts do | Follows | State |
|---|---|---|---|---|---|
| 1 | 00 | `Rattus_norvegicus.mRatBN7.2.115.gtf` | `GRCr8.114.gtf` (`:293`, `:300`, `CIRI_2_run_all_seba.sh:6`) | `GRCr8.${RELEASE}.gtf` | the Methods filename is impossible; fix the text |
| 2 | 00 | release 115 implied | every executed command is 114 | 115, by decision | settled: everything will be re-run at 115 |
| 3 | 00 | `GRCr8.dna.primary_assembly` | agrees (`:46` marks toplevel `NO USAR`) | primary assembly, `1..20 X Y MT` | **the Methods were right**; this repository was wrong and is fixed |
| 4 | 06 | guided assembly, no parameters | `--conservative`, sometimes `-f 0.05 -m 200 -t -s 3.0` (rn7-era only) | plain guided; `STRINGTIE_EXTRA` to add | the GRCr8 invocation was never recorded |
| 5 | 09 | `--window=200 --maxwindow=100000` | `--window=1000 --maxwindow=100000` (`:300`) | 1 kb | the Methods are wrong; fix the text |
| 6 | 09 | three modules, `--monoex=1` | agrees | as recorded | OK |
| 7 | 10 | decoy/gentrome index | agrees (`:232-235`) | gentrome, `--decoys`, `-k 31` | **the Methods were right**; this repository was wrong and is fixed |
| 8 | 11 | "(`-lib` and `-rat`)" | `-lib … -on rat -db mirbase` | the real flags | minor; fix the flag names |
| 9 | 13 | "CIRCexplorer3 v1.0.0" | CLEAR (= CIRCexplorer3) driving `CIRCexplorer2 parse`/`annotate` | the same thing | **the Methods are right**; only "v1.0.0" is invented and the CLEAR citation is missing |
| 10 | 13-15 | circRNA reference "Ensembl rat **rn7**" | every path is rn8/GRCr8 | GRCr8 | the Methods are wrong; fix the text |
| 11 | 15 | reconstruction implied for the analysis as a whole | `for m in {13..16}` — SQ13-SQ16 only | SQ13-SQ16, `--all` needs `--yes-widen` | **unresolved** |
| 12 | 04b | "at least two BSJ counts" | `readNumber > 2` — at least three | 3, via `--bsj-min` | **unresolved, blocking** |
| 13 | 04b | ComBat_seq not mentioned | notebooks yes, final script no | off, via `--use-combat` | **unresolved** |
| 14 | 16, 21 | IntaRNA `--mode=M --outMode=C`, MFE ≤ −20 | agrees | as recorded | OK |
| 15 | 21 | "no additional interaction tiers were used" | three tiers, floor at −10 | tiers kept as audit output; finals at −20 | soften the sentence |
| 16 | 16 | lncRNA: padj < 0.05 **and \|log2FC\| ≥ 0.5** | `absfc_cutoff <- 0` | as ported (padj alone) | **unresolved** |
| 17 | 20 | circRNA: nominal p < 0.05, \|log2FC\| ≥ 0.5 | agrees (`:26-27`) | as recorded | OK |
| 18 | 02 | mRNA: padj < 0.05, \|log2FC\| ≥ 1 | agrees | as recorded | OK |
| 19 | 02 | Pearson, \|rho\| ≥ 0.8, p < 0.05, VST | agrees | as recorded | OK |
| 20 | all | N = 3 per condition | agrees — rep 1/2/4 everywhere | n = 3 | consistent; the excluded libraries are now named (SQ5_ln, SQ6_ln, SQ3_ln, SQ4_ln) |

### Parameters the Methods omit

Not contradictions, but the Methods alone are not sufficient to reproduce
these. All are carried in the corresponding stage.

| Step | In the scripts, absent from the Methods |
|---|---|
| fastp | `-3`, the thread count; the CIRI2 branch adds `-b 140 -B 140 -l 140` |
| StringTie | `--conservative` and, in two early runs, `-f 0.05 -m 200 -t -s 3.0` (`:107`, `:109`) |
| FEELnc codpot | `--mode=intergenic`, `-b transcript_status=KNOWN` (`:195`, `:205`, `:217`, `:225`) |
| Salmon index | `-k 31` |
| BWA-MEM | `-t 20` alongside the documented `-T 19` |
| STAR | `--chimSegmentMin 10 --runThreadN 20` |

### Verified as matching, no action

fastp `-f 10 -F 10 -3 --cut_tail_mean_quality 28 -q 30` (`:66`); HISAT2
`--dta-cufflinks -S` (`:78-90`); Salmon quant `-l A --validateMappings`
(`:241-288`) and v1.10.3 (`:280`); StringTie guided then merge then gffcompare;
FEELnc maxwindow 100 000 (`:300`); miRge adapters
`AGATCGGAAGAGCACACGTCT` / `GTTCAGAGTTCTACAGTCCGACGATC`; miRWalk score ≥ 0.95
with 3'UTR/5'UTR/CDS and miRDB support; BWA-MEM `-T 19`
(`CIRI_2_run_all_seba.sh:98`); the GRN rules; all eight IntaRNA output names
against the `Interact_*.csv` glob; and every tool version in
`metadata/software_tool_inventory.tsv`.

### Manuscript-internal issues, not script-related

| # | Issue |
|---|---|
| 1 | **Figure 3B's lncRNA panel is stale** — see R10 above. Published 145/275/328; correct 148/258/336. |
| 2 | Figure 4's legend cites "Supplementary Tables 4-5" and Figure 5's "7-8", where the Results correctly cite S2/S3 and S5/S6. The Results text matches the workbook's actual sheet names; the legends need fixing. |
| 3 | "Supplementary Table 10" is cited; the workbook contains S1-S8. Should be S8 (`S8A_strict_triplets`). |
| 4 | Three Spanish editorial notes are still in the draft, plus the typo "Differential.expression" at the start of the Results. |
| 5 | `[FINAL GITHUB REPOSITORY URL]` ×2 and `[SEQUENCING ACCESSION(S)]` ×1 are unfilled. The repository URL is `https://github.com/integrativebioinformatics/suz-nrvm-hypertrophy2026`. |
| 6 | The v3 merge plan's final-validation item cites "Table 4", which does not exist in v14; assumed to mean S4 (lncRNA summary), which matches Figure 2C. **Confirm.** |

---

# Part 5 — conformance with the v3 handoff instructions

The 18 Aug 2026 `NE_NRVM_FINAL_repository_update_v3` package came with five
instruction documents, preserved under
`archive/final_package_v3_original/`. Item by item:

## `README_FIRST.md`

| Instruction | Status |
|---|---|
| Split final manuscript-aligned scripts from raw circRNA provenance | Done. `scripts/` is runnable; `archive/` is verbatim and on no runner's path. |
| "Raw circRNA provenance scripts … should be consolidated and parameterized rather than run unchanged" | Done. Bash 13-16 and R04/R04b. |
| Read the provenance audit before changing references, BSJ filters or DE logic | Done. All seven items are in `OPEN_ISSUES.md` §2, and each parameter is a CLI flag rather than an edit. |
| Do not duplicate the existing cross-reference script | Done, by a different route. The repository was rebuilt from zero, so there was nothing to preserve; it exists once, as R22. |
| Figures need not be pixel-identical | Done. See Part 6. |

## The seven audit discrepancies

All seven were independently re-verified against the code during the rebuild.
All seven hold.

| # | Discrepancy | Where it is now |
|---|---|---|
| 1 | Reference files vs the stated lock; "do not silently replace the historical paths" | Part 2 stage 00 and `OPEN_ISSUES.md` §3. Provenance paths untouched; the runnable stages take the reference as an argument. Escalated: the lock as written names a GTF that does not exist. |
| 2 | BSJ filter: prose ≥2 vs code `> 2` | Part 4 #12. **Blocking** — v14 commits to ≥2. |
| 3 | Historical DE: ComBat_seq + \|LFC\|>1 vs final nominal p + \|LFC\|≥0.5 | Part 4 #13. |
| 4 | Replicates: notebooks 4+4 vs final 3+3 | Largely resolved — see R02. Every v14 legend says N = 3 and the mRNA/lncRNA and miRNA analyses are 3v3 too; the notebooks are the outlier. Which library "rep3" was is still open. |
| 5 | CIRI-AS covers SQ13-SQ16 only | Part 2 stage 15. |
| 6 | CIRCexplorer3 in prose vs CIRCexplorer2 + CLEAR in code | **Resolved, not a discrepancy.** CLEAR IS CIRCexplorer3; CIRCexplorer2 is its dependency. Part 2 stage 13. |
| 7 | Final table lineage not explicit | R04 + R04b, labelled RECONSTRUCTED, asserting the recorded counts. |

## `handoff_summary.txt` — the five priorities

| # | Priority | Status |
|---|---|---|
| 1 | GRCr8/release-115 vs historical rn8/GRCr8.114 | Resolved as far as the evidence allows: 114 and 115 are both GRCr8, so this is an annotation-release gap, not an assembly conflict. The forward decision is 115. |
| 2 | BSJ threshold | Open, blocking. |
| 3 | circRNA DE provenance | Open. |
| 4 | 3 vs 4 replicates | Largely resolved; the manuscript says N = 3. |
| 5 | Final circRNA count/sequence table lineage | Reconstructed in R04 + R04b. |

## `repository_merge_plan.md`

**A. Compare/update rather than duplicate.** Followed in spirit. The rebuild
started from zero, so there was nothing to compare against — but the sound bash
stages 00-12 were carried over rather than rewritten, and each ported R stage
exists exactly once, diffable against its archived original.

**B. The seven recommended consolidated stages**, with the one deliberate
departure:

| v3 recommendation | This repository |
|---|---|
| `circRNA_01_circexplorer.sh` | `scripts/bash/13_circrna_circexplorer.sh` |
| `circRNA_02_ciriquant.sh` | `scripts/bash/14_circrna_ciriquant.sh` |
| `circRNA_03_ciri_as.sh` | `scripts/bash/15_circrna_ciri_as.sh` |
| `circRNA_04_annotation_merge.R` | `scripts/R/04_circRNA_annotation_merge.R` |
| `circRNA_05_final_DE_temporal.R` | `scripts/R/04b_circRNA_de_temporal.R` (temporal classification in R13) |
| `circRNA_06_intarna.sh` | `scripts/bash/16_circrna_intarna.sh` |
| the existing cross-reference script | `scripts/R/22_circRNA_cross_reference.R` |

**The departure is numbering only.** The circRNA stages are numbered into the
single global pipeline sequence rather than a separate `circRNA_0N_` series, so
`check_io_graph.R` can enforce one dependency order across all 27 R stages and
21 bash stages. The content of each follows the plan; only the filename prefix
differs, and this table is the crosswalk.

**One correction made in response to this check.** Stages 4 and 5 were
initially collapsed into a single `04_circRNA_merge_export.R`. They are now
split as the plan specifies — see R04/R04b for why the split is load-bearing
rather than cosmetic.

**C. Do not require pixel-identical figures.** Done; see Part 6.

**D. Final validation before freezing.**

| Item | Status |
|---|---|
| Reference version and checksum | Version resolved. **Checksums still `TO_FILL`** in `reference_manifest.tsv`. |
| BSJ threshold | Open, blocking. |
| circRNA sample inclusion | Largely resolved (N = 3); which library "rep3" was is open. |
| Whether ComBat_seq belongs in the final branch | Open. |
| Final counts in Figures 2-3 and "Table 4" | Figure 2B **verified exactly**. Figure 3A/3C and S3C reconcile; **Figure 3B's lncRNA panel does not**. "Table 4" does not exist in v14 — see Part 4, manuscript issue 6. |
| IntaRNA filenames and MFE ≤ −20 filtering | Verified. One caveat: Part 4 #15. |

---

# Part 6 — what the pipeline does not reproduce

From `archive/final_package_v3_original/04_METADATA/manual_editing_provenance.tsv`.

**Affinity — multi-panel assembly.** Every multi-panel figure was composed and
cosmetically edited in Affinity. The pipeline guarantees the numerical content,
the selected features, and each source panel as its own PDF/PNG pair. It does
not guarantee the page composition, panel lettering, cropping, shared legends
or typography. The audit's own wording is the right standard: the analytical
scripts should be expected to reproduce the underlying source plots and data,
not necessarily the exact final page composition.

*To close this gap:* deposit the Affinity sources (`.afdesign`) or exported
SVG/PDF of the final compositions.

**Cytoscape v3.10.3 — network layouts.** Figures 6A-B and S5B-F are Cytoscape
renderings whose layout was adjusted by hand. The pipeline guarantees network
membership and every node and edge attribute, in import-safe tables: R24 writes
`ALL_MFE20_nodes.csv`, `ALL_MFE20_edges.csv`, the per-module tables, a
visual-style template and an export manifest; R25 writes the
`*_CYTOSCAPE_SAFE.tsv` pair and an ID-match report. It does not guarantee node
coordinates — re-importing reproduces the same graph with a different layout.

*To close this gap:* deposit the saved `.cys` session, or the final node x/y
coordinates plus the visual style.

**Not transcriptomic at all.** Figure 1A-D (workflow, immunofluorescence,
area/perimeter), Figure 2A (hand-drawn integrated workflow diagram) and Figure
S1A-E (phenotypic, qPCR and western validation). Preserve the editable source
for Figure 2A if it exists.

**One notebook's panels exist only as HTML.** The visualization notebook from
the metadata analysis contains more than forty characterization panels and
**zero `ggsave()` calls** — every one of them exists only inside the rendered
HTML. If any is wanted for the supplement it has to be re-emitted through
`save_figure()` so there is a file on disk. That notebook is not in the
repository.

---

# Part 7 — defects preserved in `archive/`, not reproduced

Real defects in the files under `archive/`. They are kept verbatim because that
directory is the historical record. The active stages do not repeat them.

| Where | Defect |
|---|---|
| `03_tool_comparation.Rmd:28-30` | `c3_cq <- length(c3_cq$chr)` overwrites the data frame with a scalar, then `c3_cq$circRNA_ID <- …` is attempted on an atomic vector. **The chunk cannot run past line 30 as committed.** |
| `03_tool_comparation.Rmd:169,173` | `circRNA_sequences_de{6,24}_2.fa` written with `write.table(sep="\t", col.names=TRUE)` — TSVs with a `.fa` extension. |
| `01_circexp3.Rmd:159` | Writes the CIRCexplorer 6 h normalised matrix into the **CIRIquant** results folder. |
| `01_circexp3.Rmd:244` | Bypasses `paste0(path, …)` with a hard-coded `D:` path. |
| `00_annotation_biomart.Rmd:12` | Carries `#NO VOLVER A EJECUTAR` — do not re-run. Honoured: the stage is not in the runnable path. |
| `00_annotation_biomart.Rmd:47,129` | `getBM(values = …)` with no `filters = …`; biomaRt ignores `values`, so these chunks never retrieved what they appear to. |
| `vis_run_ciri2_seba.sh:3` | Sample list reads `Q13 Q14 Q15 Q16` instead of `SQ13…SQ16`. |
| `vis_run_ciri2_seba.sh` | `rm ${sample}_vis/*.pdf` — deletes the PDFs it just generated. |
| `subset_merged_gtfs.sh` | No `case` arm for `CNE_24_STAR.annotated`; a correctly-named file falls through to `*)` and is skipped with a warning only. |
| `filter_out_lncRNAs.sh:26` | Requires GNU awk's three-argument `match()`. Under mawk the filter silently passes every line through. |
