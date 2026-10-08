# suz-nrvm-hypertrophy2026

Reproducibility pipeline for the NE-NRVM manuscript: **norepinephrine versus
control in neonatal rat ventricular myocytes, at 6 h and 24 h**, integrating
mRNA, lncRNA, miRNA and circRNA into ceRNA networks.

- **6 main figures + 5 supplementary figures** — see `metadata/figure_script_map.tsv`
- **Linux-first.** Windows via WSL2 or Docker, not native paths.
- **Reference lock:** `Rattus_norvegicus.GRCr8.dna.toplevel.fa` +
  `Rattus_norvegicus.GRCr8.115.gtf`, Ensembl release 115, biomaRt pinned to
  `https://sep2025.archive.ensembl.org`.

> **Before treating any output as final, read [`docs/MANUSCRIPT_ALIGNMENT.md`](docs/MANUSCRIPT_ALIGNMENT.md)
> and [`docs/OPEN_ITEMS.md`](docs/OPEN_ITEMS.md).**
> Two pipeline steps are reconstructions of scripts that were never supplied,
> and several thresholds in the Methods text do not match the code.

---

## Quick start

```bash
# 1. environments
conda env create -f envs/environment_qc_align.yml
conda env create -f envs/environment_mirge.yml
conda env create -f envs/environment_intarna.yml
conda env create -f envs/environment_circrna.yml
conda env create -f envs/environment_r.yml

# 2. check the repository is internally consistent (seconds, no data needed)
bash scripts/checks/check_biomart_pinning.sh
bash scripts/checks/check_hardcoded_paths.sh
Rscript scripts/checks/check_io_graph.R

# 3. confirm the pinned Ensembl archive is reachable and is release 115
Rscript scripts/R/01_transcriptome_metadata_annotation.R --dry-run-biomart

# 4. run
./suz_pipeline.bash \
  --longrna-samples metadata/longRNA_samples.tsv \
  --smallrna-samples metadata/smallRNA_samples.tsv \
  --fastq-manifest metadata/fastq_manifest.tsv \
  --threads 16 --out results --run-r-analysis
```

Add `--skip-circrna` to run everything except the circRNA module, which needs
external resources this repository cannot ship (see
`metadata/circRNA_input_manifest.tsv`).

Every stage has `--help`.

---

## Pipeline

### Bash stages

| # | Script | Env | Purpose |
|---|---|---|---|
| 00 | `00_download_references.sh` | qc_align | Fetch GRCr8 toplevel FASTA + GRCr8.115 GTF |
| 01 | `01_qc_raw.sh` | qc_align | FastQC on raw reads |
| 02 | `02_trim_fastp.sh` | qc_align | fastp trimming |
| 03 | `03_qc_trimmed.sh` | qc_align | FastQC on trimmed reads |
| 04 | `04_align_hisat2.sh` | qc_align | HISAT2 alignment |
| 05 | `05_sort_index_bam.sh` | qc_align | samtools sort/index |
| 06 | `06_stringtie_assemble_merge.sh` | qc_align | StringTie assembly and merge |
| 07 | `07_gffcompare.sh` | qc_align | Annotate against the reference; produces the `.tmap` |
| 08 | `08_gffread_extract_fasta.sh` | qc_align | Transcript FASTA |
| 09 | `09_feelnc_classify.sh` | qc_align | FEELnc filter → codpot → classifier |
| **09b** | `09b_feelnc_split_exon_gtfs.sh` | qc_align | **RECONSTRUCTED** — the `exons/` layout |
| 09c | `09c_extract_transcript_ids.sh` | qc_align | Transcript ID lists |
| 09d | `09d_subset_merged_gtfs.sh` | qc_align | Full records for those IDs |
| 09e | `09e_filter_out_lncRNAs.sh` | qc_align | lncRNA-depleted merged GTF (**needs gawk**) |
| 10 | `10_salmon_index_quant.sh` | qc_align | Salmon quantification |
| 11 | `11_mirge3_run.sh` | mirge | miRge3.0 (`-db mirbase`) |
| 12 | `12_intarna_run.sh` | intarna | lncRNA–miRNA duplexes (4 commands) |
| 13 | `13_circrna_circexplorer.sh` | circrna | CIRCexplorer2 + CLEAR |
| 14 | `14_circrna_ciriquant.sh` | circrna | CIRIquant |
| 15 | `15_circrna_ciri_as.sh` | circrna | BWA-MEM → CIRI2 → CIRI-AS |
| 16 | `16_circrna_intarna.sh` | intarna | circRNA–miRNA duplexes (4 commands) |

Stages 09b–09e carry the transcript-metadata chain that previous versions of
this repository omitted entirely; `scripts/R/01` cannot run without them.

### R stages

Numbered in dependency order. `scripts/checks/check_io_graph.R` enforces that
no stage reads a file a later stage produces.

| # | Stage | Panels |
|---|---|---|
| 01 | `transcriptome_metadata_annotation` | — (produces the `meta_*` tables; **the only biomaRt call**) |
| 02 | `deseq2_mRNA_lncRNA` | — (DE tables, temporal phases, correlations, hubs, KEGG) |
| 03 | `mirna_deseq2` | — (miRNA DE and temporal tables) |
| 04 | `circRNA_merge_export` | — **RECONSTRUCTED** |
| 05 | `tf_kegg_integration` | **Fig 4A-B, S4A** |
| 06 | `transcriptome_overview` | **Fig 2B** |
| 07 | `lncRNA_characterization_gtf` | **Fig 2C, S2A-C** |
| 08 | `feelnc_panels` | S2 support |
| 09 | `de_counts_heatmaps_gba` | **Fig 3A, 3C miRNA, 5A, S4B** |
| 10 | `figure3B_temporal_programs` | **Fig 3B** |
| 11 | `figure3C_mRNA_heatmap` | **Fig 3C mRNA** |
| 12 | `figure3C_lncRNA_heatmap` | **Fig 3C lncRNA** |
| 13 | `figure3C_circRNA_heatmap` | **Fig 3C circRNA** — NEW |
| 14 | `mirna_mRNA_integration` | **Fig 5B-C**; emits the miRNA FASTAs |
| 15 | `lncRNA_hubs` | **Fig S4C** |
| 16 | `volcanoes_lncRNA` | **Fig S3A-B lncRNA** |
| 17 | `volcanoes_mRNA_miRNA` | **Fig S3A-B mRNA/miRNA** — NEW |
| 18 | `directional_overlap` | **Fig S3C** |
| 19 | `extract_lncRNA_for_intarna` | — (IntaRNA targets) |
| 20 | `circRNA_volcanoes` | **Fig S3A-B circRNA** |
| 21 | `ceRNA_coherent_pairs_triplets` | Fig 6B-F input |
| 22 | `circRNA_cross_reference` | **Fig S3C circRNA**, cross-ncRNA modules |
| 23 | `ceRNA_phase_aware_plots` | **Fig 6A** |
| 24 | `ceRNA_cytoscape_export` | **Fig 6B-F, S5A-B** |
| 25 | `cytoscape_key_validation` | QC for stage 24 |
| 26 | `results_report` | Audit trail |

---

## Layout

```
├── suz_pipeline.bash              top-level orchestrator
├── envs/                          5 conda environments
├── metadata/                      sample sheets, manifests, the figure map,
│                                  the stage I/O contract, checksums
├── data/frozen/                   inputs that cannot be regenerated (see its README)
├── references/                    populated by stage 00; not committed
├── scripts/
│   ├── bash/                      stages 00-16 + helpers
│   ├── R/                         stages 01-26 + lib/
│   ├── checks/                    repository self-checks
│   ├── commands/                  checksum helpers
│   └── dev/                       the archive→stage porter
├── provenance/                    upstream scripts, verbatim, never executed
├── archive/                       everything as received, verbatim
└── docs/                          OPEN_ITEMS, REFERENCE_AUDIT, METHODS_NOTES,
                                   METHODS_NOTES, MANUAL_ARTIFACTS
```

---

## Things worth knowing before you run anything

**Two stages are reconstructions.** `09b` rebuilds the step that produces the
`exons/` layout; `R04` rebuilds the circRNA count-matrix assembly, tool merge
and export. Neither original was ever supplied. Both say so loudly in their
headers and in their output. `scripts/checks/check_feelnc_split_counts.sh`
verifies 09b against the historical transcript counts and currently passes
exactly.

**Sample size is n = 3, not n = 4.** Replicate 3 is excluded from every final
analysis — the mRNA/lncRNA DESeq2, the miRNA matrices and the circRNA volcanoes
all use three per group. The sample sheets list 16 long-RNA libraries; 12 were
analysed. Which physical library "rep3" was is not recorded anywhere
(`docs/OPEN_ITEMS.md`, A3 and A4).

**circRNA DE is exploratory.** Nominal p, not FDR.

**gawk is required** by stage 09e. Under mawk the lncRNA filter silently passes
every line through and produces an unfiltered "filtered" GTF. The stage checks
and aborts.

**Figures are not pixel-reproducible.** Affinity assembly and manual Cytoscape
layout sit between these outputs and the manuscript. See
`docs/MANUAL_ARTIFACTS.md`.

---

## Self-checks

```bash
bash scripts/checks/check_biomart_pinning.sh     # one pinned Ensembl connection, no version= 
bash scripts/checks/check_hardcoded_paths.sh     # no foreign paths, no name(1).csv duplicates
Rscript scripts/checks/check_io_graph.R          # every input has a declared source
bash scripts/checks/check_feelnc_split_counts.sh # the 09b reconstruction still matches
bash scripts/commands/make_checksums.sh --verify # frozen inputs and provenance intact

for f in scripts/R/[0-9]*.R; do Rscript "$f" --help >/dev/null || echo "FAIL $f"; done
for f in scripts/bash/*.sh suz_pipeline.bash; do bash -n "$f" || echo "FAIL $f"; done
```

---

## What is not here

- **FASTQ files.** On the server and an external drive. Paths, sizes and MD5s
  for the 32 long-RNA files are in `metadata/fastq_manifest.tsv`; the 16
  small-RNA rows are still `TO_FILL`.
- **Large references.** Downloaded by stage 00 or supplied via `--genome-fasta`
  and `--gtf`.
- **External resources** the circRNA branch and the miRNA integration need:
  `rn8.yml`, `rn8.refflat.nochr.txt`, a STAR index, `CIRI2.pl`, `CIRI_AS.pl`,
  `TFLink_All.tsv`, the 12 miRWalk target files, miRBase `mature.fa`. All listed
  in `metadata/circRNA_input_manifest.tsv` and `metadata/external_inputs.tsv`.
- **The Quarto notebooks** from the metadata analysis. Their logic is ported to
  `scripts/R/01`; the originals stay in `SUZ2026/Babi analysis/`. See
  `archive/ARCHIVE_NOTES.md`.
