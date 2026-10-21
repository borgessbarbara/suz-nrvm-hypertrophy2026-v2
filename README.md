# suz-nrvm-hypertrophy2026

Reproducibility pipeline for the NE-NRVM manuscript: **norepinephrine versus
control in neonatal rat ventricular myocytes at 6 h and 24 h**, integrating
mRNA, lncRNA, miRNA and circRNA into ceRNA networks.

6 main figures and 5 supplementary figures. 21 bash stages, 27 R stages, one
orchestrator.

### The three documents

| | |
|---|---|
| **This file** | what the pipeline is, how it is structured, how to run it |
| [`docs/CHANGES_FROM_ORIGINALS.md`](docs/CHANGES_FROM_ORIGINALS.md) | every change made to the original scripts, in pipeline order, each pointing at its original in `archive/` |
| [`docs/OPEN_ISSUES.md`](docs/OPEN_ISSUES.md) | what is missing, unresolved or waiting on a decision |

Read `OPEN_ISSUES.md` before treating any output as final. Two stages are
reconstructions of scripts that were never supplied, and several thresholds in
the manuscript Methods do not match the code.

---

## Reference lock

| | |
|---|---|
| Genome | `Rattus_norvegicus.GRCr8.dna.primary_assembly.fa` (chromosomes `1..20 X Y MT`, concatenated) |
| Annotation | `Rattus_norvegicus.GRCr8.116.gtf` |
| Ensembl release | **116**  (114 and 115 also selectable) |
| biomaRt | `https://jun2026.archive.ensembl.org` (= `ensembl_mart_116`) |

The release is set in **one place**, `metadata/reference_release.tsv`. Its
`active` column drives the GTF path, the FASTA URLs and the biomaRt archive
together, so the annotation and the genome can never drift apart. Override for
a single run with `ENSEMBL_RELEASE=114` or `=115`.

```bash
# which release will actually be used, and is the archive reachable?
Rscript scripts/R/01_transcriptome_metadata_annotation.R --dry-run-biomart
```

**Linux-first.** Windows via WSL2 or Docker, not native paths. `.gitattributes`
forces LF on every script; CRLF breaks the shebang lines.

---

## Before you run

You need three things this repository cannot contain.

**1. The FASTQ files.** On the server and an external drive. Paths, sizes and
MD5s for the 32 long-RNA files are in `metadata/fastq_manifest.tsv`; the 16
small-RNA rows are still `TO_FILL`.

**2. The references.** `scripts/bash/00_download_references.sh` fetches them, or
pass your own with `--genome-fasta` and `--gtf`.

**3. External resources** listed in `metadata/external_inputs.tsv` — 29 rows.
The circRNA branch needs `rn8.yml`, `rn8.refflat.nochr.txt`, a STAR index,
`CIRI2.pl` and `CIRI_AS.pl`; the miRNA and TF integration need
`TFLink_All.tsv`, the 12 miRWalk target files and miRBase `mature.fa`. Each
stage that needs one of these fails with a pointer rather than guessing.

Add `--skip-circrna` to run everything that does not depend on them.

---

## Install and verify

```bash
bash scripts/commands/create_envs.sh            # all five, then verify and lock
bash scripts/commands/create_envs.sh --dry-run  # solve only, change nothing
```

That creates `env_qc_align`, `env_mirge`, `env_intarna`, `env_circrna` and
`env_r`, sets the strict channel priority bioconda requires, checks that every
binary the stages call is on `PATH` inside each environment, and writes
`envs/*.lock.yml`. Or create one at a time:

```bash
conda env create -f envs/environment_qc_align.yml   # -> env_qc_align
```

**The environment files pin no versions**, so you get current releases; the
version each tool had in the manuscript Methods is a trailing comment on its
line. The `.lock.yml` files are the version record — update
`metadata/software_tool_inventory.tsv` from them once the environments are
built.

**Three tools are not conda packages.** CLEAR (`run3.py`, `circ_quant.py`)
installs from source into `env_circrna`; CIRI2 and CIRI-AS are standalone Perl
scripts passed to stage 15 with `--ciri2-path` and `--ciri-as-path`. See the
header of `envs/environment_circrna.yml`.

### Or use Docker

The `Dockerfile` builds one image carrying all five environments, with CLEAR
already installed. It runs the repository's own checks during the build, so a
broken image fails rather than ships.

```bash
docker build --platform linux/amd64 -t suz-nrvm:latest .
```

`--platform linux/amd64` matters: bioconda publishes linux-64 only, so an
Apple-silicon Mac has to emulate. Expect a long first build and several GB —
`--build-arg ENVS="qc_align r"` builds a subset.

```bash
docker run --rm -it --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  -v /path/to/fastq:/data/fastq:ro \
  -v /path/to/references:/data/references:ro \
  -v "$PWD/results:/opt/suz/results" \
  suz-nrvm:latest bash
```

Inside the container everything works exactly as on a host — the stages
activate their own environments. CIRI2 and CIRI-AS are not redistributable
through conda, so mount them and pass `--ciri2-path` / `--ciri-as-path`.

Then confirm the repository is internally consistent. All of these run in
seconds and need no data:

```bash
bash scripts/checks/check_biomart_pinning.sh      # one pinned Ensembl connection
bash scripts/checks/check_hardcoded_paths.sh      # no foreign paths, no name(1).csv
Rscript scripts/checks/check_io_graph.R           # every input has a declared source
bash scripts/checks/check_feelnc_split_counts.sh  # the stage 09b reconstruction
bash scripts/commands/make_checksums.sh --verify   # frozen inputs intact

for f in scripts/R/[0-9]*.R; do Rscript "$f" --help >/dev/null || echo "FAIL $f"; done
for f in scripts/bash/*.sh suz_pipeline.bash; do bash -n "$f" || echo "FAIL $f"; done
```

---

## Run

Whole pipeline:

```bash
./suz_pipeline.bash \
  --longrna-samples metadata/longRNA_samples.tsv \
  --smallrna-samples metadata/smallRNA_samples.tsv \
  --fastq-manifest metadata/fastq_manifest.tsv \
  --threads 16 --out results --run-r-analysis
```

Just the R analysis, or a slice of it:

The bash stages activate their own environment. The R stages do not, so
activate `env_r` first:

```bash
conda activate env_r
bash scripts/bash/_run_r_analysis.sh --results results --dpi 600
bash scripts/bash/_run_r_analysis.sh --stages 5-13          # a range
bash scripts/bash/_run_r_analysis.sh --stages 21,22,24      # a selection
bash scripts/bash/_run_r_analysis.sh --stages 14- --dry-run # from 14 to the end
```

One stage on its own:

```bash
Rscript scripts/R/06_transcriptome_overview.R --in-dir results --out-dir results --dpi 600
bash scripts/bash/09_feelnc_classify.sh
```

**Every stage has `--help`.** R stages take `--in-dir`, `--out-dir` and
`--dpi`; bash stages are configured by environment variables, named in each
script's header. Nothing has a path baked in.

---

## The pipeline, in run order

### Bash stages

| # | Script | Env | What it does |
|---|---|---|---|
| 00 | `00_download_references.sh` | qc_align | GRCr8 primary-assembly FASTA + GRCr8.116 GTF |
| 01 | `01_qc_raw.sh` | qc_align | FastQC on raw reads |
| 02 | `02_trim_fastp.sh` | qc_align | fastp trimming |
| 03 | `03_qc_trimmed.sh` | qc_align | FastQC on trimmed reads |
| 04 | `04_align_hisat2.sh` | qc_align | HISAT2 index + alignment (`--dta-cufflinks`) |
| 05 | `05_sort_index_bam.sh` | qc_align | samtools sort and index |
| 06 | `06_stringtie_assemble_merge.sh` | qc_align | StringTie per sample, then `--merge` |
| 07 | `07_gffcompare.sh` | qc_align | class codes and the `.tmap` |
| 08 | `08_gffread_extract_fasta.sh` | qc_align | transcript FASTAs |
| 09 | `09_feelnc_classify.sh` | qc_align | FEELnc filter → codpot → classifier |
| 09b | `09b_feelnc_split_exon_gtfs.sh` | qc_align | the `exons/` layout — **RECONSTRUCTED** |
| 09c | `09c_extract_transcript_ids.sh` | qc_align | transcript ID lists |
| 09d | `09d_subset_merged_gtfs.sh` | qc_align | full gffcompare records for those IDs |
| 09e | `09e_filter_out_lncRNAs.sh` | qc_align | lncRNA-depleted merged GTF — **needs gawk** |
| 10 | `10_salmon_index_quant.sh` | qc_align | decoy-aware gentrome index, then quant |
| 11 | `11_mirge3_run.sh` | mirge | miRge3.0, rat miRBase |
| 12 | `12_intarna_run.sh` | intarna | 4 lncRNA–miRNA duplex runs |
| 12b | `12b_circrna_references.sh` | circrna | STAR index, refFlat, BWA index and `.fai` for the circRNA branch |
| 13 | `13_circrna_circexplorer.sh` | circrna | CIRCexplorer2 + CLEAR |
| 14 | `14_circrna_ciriquant.sh` | circrna | CIRIquant |
| 15 | `15_circrna_ciri_as.sh` | circrna | BWA-MEM → CIRI2 → CIRI-AS |
| 16 | `16_circrna_intarna.sh` | intarna | 4 circRNA–miRNA duplex runs |

Stages `09b`–`09e` are the transcript-metadata chain. **`scripts/R/01` cannot
run without them**, and they run before any R stage.

Helpers: `_activate_env.sh` (non-interactive conda activation, sourced by every
stage) and `_run_r_analysis.sh` (the R stage runner).

### R stages

Numbered in dependency order. `scripts/checks/check_io_graph.R` enforces that
no stage reads a file a later stage produces.

| # | Stage | Produces |
|---|---|---|
| 01 | `transcriptome_metadata_annotation` | the `meta_*` tables — **the only biomaRt call** |
| 02 | `deseq2_mRNA_lncRNA` | mRNA/lncRNA DE, temporal phases, correlations |
| 03 | `mirna_deseq2` | miRNA DE and temporal tables |
| 04 | `circRNA_annotation_merge` | locus merge, unfiltered BSJ matrix — **RECONSTRUCTED** |
| 04b | `circRNA_de_temporal` | circRNA DE and exports — **RECONSTRUCTED** |
| 05 | `tf_kegg_integration` | **Fig 4A-B** |
| 06 | `transcriptome_overview` | **Fig 2B** |
| 07 | `lncRNA_characterization_gtf` | **Fig 2C, S2A-C** |
| 08 | `feelnc_panels` | Figure S2 support tables |
| 09 | `de_counts_heatmaps_gba` | **Fig 3A, 3C miRNA, 5A, S4A** |
| 10 | `figure3B_temporal_programs` | **Fig 3B** |
| 11 | `figure3C_mRNA_heatmap` | **Fig 3C mRNA** |
| 12 | `figure3C_lncRNA_heatmap` | **Fig 3C lncRNA** |
| 13 | `figure3C_circRNA_heatmap` | **Fig 3C circRNA** |
| 14 | `mirna_mRNA_integration` | **Fig 5B-C**, and the miRNA FASTAs |
| 15 | `lncRNA_hubs` | lncRNA hub tables |
| 16 | `volcanoes_lncRNA` | **Fig S3A-B lncRNA** |
| 17 | `volcanoes_mRNA_miRNA` | **Fig S3A-B mRNA/miRNA** |
| 18 | `directional_overlap` | **Fig S3C** mRNA/lncRNA/miRNA |
| 19 | `extract_lncRNA_for_intarna` | lncRNA target FASTAs for stage 12 |
| 20 | `circRNA_volcanoes` | **Fig S3A-B circRNA** |
| 21 | `ceRNA_coherent_pairs_triplets` | coherent pairs and triplets |
| 22 | `circRNA_cross_reference` | **Fig S3C circRNA**, cross-ncRNA modules |
| 23 | `ceRNA_phase_aware_plots` | **Fig S5A** |
| 24 | `ceRNA_cytoscape_export` | **Fig 6A-B, S5B-F** node/edge tables |
| 25 | `cytoscape_key_validation` | import-key QC for stage 24 |
| 26 | `results_report` | the audit trail |

Stages 12 and 16 (IntaRNA) sit **between** R19/R14 and R21: R19 and R14 write
the FASTAs, the two bash stages predict the duplexes, R21 reads them back.
`suz_pipeline.bash` handles that ordering.

### Figure map

`metadata/figure_script_map.tsv` has all 41 panels with their content and
notes. In summary:

| Figure | Panels | Source |
|---|---|---|
| 1, 2A, S1 | workflow, microscopy, qPCR, western | experimental; not this pipeline |
| 2B-C | biotype inventory, lncRNA inventory | R06, R07 |
| 3A-C | DE counts, temporal programs, heatmaps | R09, R10, R11, R12, R13 |
| 4A-B | curated KEGG GSEA, TF-by-pathway | R05 |
| 5A-C | guilt-by-association, miRNA–mRNA maps | R09, R14 |
| 6A-B | the two integrated ceRNA networks | R21 → R22 → R24 → R25 |
| S2A-C | lncRNA exon count, isoform complexity, length | R07 |
| S3A-C | volcanoes, directional overlap | R16, R17, R20, R18, R22 |
| S4A | curated KEGG among lncRNA–mRNA correlations | R09 |
| S5A-F | phase trajectories, five phase networks | R23, R24 |

---

## Layout

```
suz-nrvm-hypertrophy2026/
├── README.md                     this file
├── suz_pipeline.bash             orchestrator
├── Dockerfile                    one image, all five environments
├── docs/
│   ├── CHANGES_FROM_ORIGINALS.md every change to the original scripts
│   └── OPEN_ISSUES.md            what is unresolved
├── envs/                         5 conda environments (+ .lock.yml after setup)
├── metadata/                     12 tables: sample sheets, manifests,
│                                 the figure map, the stage I/O contract,
│                                 the reference lock, checksums
├── data/frozen/                  2 inputs that cannot be regenerated
├── references/                   populated by stage 00, not committed
├── scripts/
│   ├── bash/                     stages 00-16 + 2 helpers
│   ├── R/                        stages 01-26 + lib/
│   ├── checks/                   4 repository self-checks
│   └── commands/                 2 checksum helpers
└── archive/                      every package received, verbatim,
                                  never executed
```

### `metadata/`

| File | Holds |
|---|---|
| `reference_release.tsv` | the active Ensembl release — the one switch |
| `reference_manifest.tsv` | reference file names, URLs, checksums |
| `longRNA_samples.tsv`, `smallRNA_samples.tsv` | the two sample sheets |
| `fastq_manifest.tsv` | 48 FASTQ rows with paths, sizes, MD5s |
| `external_inputs.tsv` | 29 inputs no stage can produce, with who has them |
| `figure_script_map.tsv` | 41 panels → producing stage |
| `supplementary_tables_map.tsv` | 19 workbook sheets → producing stage |
| `stage_io.tsv` | the I/O contract `check_io_graph.R` enforces |
| `software_tool_inventory.tsv` | every tool and version, with its evidence |
| `curated_kegg_allowlist.tsv` | the curated pathway allowlist, read by R09 and R24 |
| `frozen_inputs_checksums.tsv` | 104 files under `data/` and `archive/` |

### `archive/`

Everything received, byte-for-byte, never executed and never on any runner's
path. `archive/ARCHIVE_NOTES.md` says what each directory is and where each
script went.

| Directory | Files | What |
|---|---|---|
| `raw_scripts_original/` | 26 | the raw script collection, including `bash_scripts/bash_scripts.txt`, the command log |
| `final_package_v3_original/` | 44 | the 18 Aug 2026 final package: final figure scripts, circRNA provenance, IntaRNA commands, the audit |
| `circrna_addendum_v1_original/` | 6 | the 31 Jul 2026 circRNA addendum |
| `starter_v1_original/` | 11 | the 7 Jul 2026 starter package |
| `transcriptome_metadata_original/` | 3 | the three shell scripts of the metadata chain |

---

## Conventions every stage follows

**Inputs resolve by basename, not by path.** `find_input()` in
`scripts/R/lib/io.R` searches `--in-dir` recursively, then the working
directory. Pass a full path and it goes through unchanged. Browser-download
duplicates (`name(1).csv`) are rejected outright.

**Each stage owns its output directory.** `results/R<NN>_<stage>/`. No stage
writes into another's. Every output directory gets a
`00_INPUT_manifest.csv` recording which file was resolved for each role, with
its size and mtime.

**Figures are a PDF/PNG pair.** `save_figure()` writes `<stem>.pdf` (cairo
where available) and `<stem>_<dpi>dpi.png`. Default 600 dpi, `--dpi` to change.

**One Ensembl connection.** `connect_ensembl()` in `scripts/R/lib/ensembl.R` is
the only entry point, pinned by release rather than by `version=`.
`check_biomart_pinning.sh` fails the build on any other call.

**Shared palettes.** `scripts/R/lib/plotting.R` holds the direction, RNA-class
and temporal-phase colours, so panels assembled side by side agree.

---

## Things worth knowing before you run anything

**Two stages are reconstructions.** `09b` rebuilds the step producing the
`exons/` layout; `R04` + `R04b` rebuild the circRNA count-matrix assembly, tool
merge and export. Neither original was ever supplied. Both say so in their
headers and in their console output.
`scripts/checks/check_feelnc_split_counts.sh` verifies 09b against the
historical transcript counts and currently passes exactly.

**Sample size is n = 3, not n = 4.** Replicate 3 is excluded from every final
analysis — the mRNA/lncRNA DESeq2, the miRNA matrices and the circRNA volcanoes
all use three per group. The sample sheets list 16 long-RNA libraries; 12 were
analysed. Which physical library "rep3" was is not recorded anywhere.

**gawk is required** by stage 09e. Under mawk the lncRNA filter silently passes
every line through and produces an unfiltered GTF that looks correct. The stage
checks and aborts.

**circRNA differential expression is exploratory.** Nominal p, not FDR.

**Figures are not pixel-reproducible.** Multi-panel assembly happened in
Affinity and the Cytoscape v3.10.3 layouts were adjusted by hand. The pipeline
reproduces the numerical content, the selected features, each source panel as
its own PDF/PNG pair, and every node and edge attribute — not the page
composition or the node coordinates.

**No full pipeline run has been executed.** The repository is verified
structurally by the eight checks above, and numerically only for Figure 2B,
which reproduces exactly from the frozen tables.
