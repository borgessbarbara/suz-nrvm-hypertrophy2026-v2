# Pending and unsolved

What is not finished, grouped by whether it blocks a run, blocks the
manuscript, or is a loose end. Each item says who can close it.

Reference release is **not** on this list: release 115 with the GRCr8 primary
assembly is a settled forward decision (B. Borges, 8 Oct 2026) — everything
will be re-run and the results compared.

---

## 1. Blocks running the circRNA branch

Nothing in the circRNA module can run end-to-end until these arrive. Everything
downstream of them is also blocked: Figure 2B's circRNA bar (8,736 loci), all
circRNA panels in Figures 3 and S3, Supplementary Tables S7A-B, and the 9
circRNA triplets in Figure 6B.

| # | Missing | Needed by | Who |
|---|---|---|---|
| 1.1 | `rn8.yml` — CIRIquant config | bash 14 (fails loudly) | J. Pablo |
| 1.2 | `rn8.refflat.nochr.txt` | bash 13 | J. Pablo |
| 1.3 | STAR index `STAR_rn8/` | bash 13 | J. Pablo |
| 1.4 | `CIRI2.pl` v2.0.6, `CIRI_AS_v1.2.pl` (not conda packages) | bash 15 | J. Pablo |
| 1.5 | `rat_primary_assembly_orderr8.fa` — verify it is the GRCr8 primary assembly, by checksum | bash 13, 15 | J. Pablo |

## 2. Blocks running the rest of the pipeline

| # | Missing | Needed by | Who |
|---|---|---|---|
| 2.1 | The 12 `miRWalk_miRNA_Targets_*.csv` files | R14, R24 | S. Urquiza |
| 2.2 | miRBase `mature.fa` + the exact release used | R14 | S. Urquiza |
| 2.3 | `TFLink_All.tsv` + version/date | R05 | S. Urquiza |
| 2.4 | The three Cytoscape GRN edge exports (`GRN_NE06_*`, `GRN_NE24_*`, `Merged_Network_*`) — manual exports, not reproducible from the pipeline | R05 | S. Urquiza |
| 2.5 | `DESmetadata_24h.csv` (miRge `-mdt`) | bash 11 | S. Urquiza |
| 2.6 | Small-RNA FASTQ paths, sizes and MD5s — 16 rows still `TO_FILL` in `metadata/fastq_manifest.tsv` | bash 11 | S. Urquiza |
| 2.7 | `miRDB target prediction data_4miRNAs_correlation.xls.xlsx` (optional) | historical only | — |

## 3. Commands that were never recorded

These were reconstructed or inferred. They are the places where a re-run could
legitimately diverge from the original.

| # | Gap | What the repo does | Confidence |
|---|---|---|---|
| 3.1 | **The `exons/` split** — the step producing `<id>.annotated{,.lncRNA,.otherRNAs}.gtf` is in neither the Babi folder nor `bash_scripts.txt` | Reconstructed as bash **09b** from the FEELnc_codpot outputs | **High.** `check_feelnc_split_counts.sh` reproduces 46,030 / 42,127 / 3,822 exactly, with the 87 classifier orphans and 81 noORF residuals, and the sets are disjoint |
| 3.2 | **`FEELnc_codpot` for the final rn8 run** — `bash_scripts.txt` has the heading `M2` with **nothing under it** (lines 295-297). Only Module 1 and Module 3 were recorded. | Stage 09 uses the parameters from the four earlier recorded codpot runs (`--mode=intergenic -b transcript_biotype=protein_coding -b transcript_status=KNOWN`) | **Medium.** Parameters were identical across all four earlier runs, so carrying them forward is reasonable — but unverified for the final one |
| 3.3 | **circRNA count matrix / tool merge / export** — nothing produces `raw_count_matrix_all_*`, `0{6,7}_circRNA_sequences_de*_DE_complete_table.txt` or the four `circRNA_*_all_records.fa` | Reconstructed as **R04** + **R04b** from `03_tool_comparation.Rmd` | **Medium.** Asserts the recorded 46,302 / 3,419 / 23,739 locus counts and warns on mismatch |
| 3.4 | **rn8 HISAT2 index build** — only the rn7 builds are recorded (`:55`, `:58`) | Stage 04 builds its own index | Low risk |
| 3.5 | **rn8 StringTie merge + gffcompare** — only the rn7-era commands are recorded (`:138`, `:140`, `:158`, `:190`) | Stages 06-07 use the documented guided parameters | Low risk |
| 3.6 | **Which candidate GTF was final** — four successive FEELnc runs on `candidate_lncRNA_NE6_NE24.gtf` → `candidate_lncRNAs_NE6_NE24_oct.gtf` → `cand_lncRNAs_NE6_NE24_annotated.gtf` → `cand_lncRNAs_ALLPE_annotated.gtf` | The rn8 run (`:293`) uses `merged_CNE624_8_rn8_H.annotated.gtf`, which supersedes all four | Resolved for rn8 |
| 3.7 | **Primary-assembly concatenation order** — `rat_primary_assembly_orderr8.fa` implies a specific one | `1..20 X Y MT`, `nonchromosomal` excluded by decision | Affects BAM header order, not alignment |
| 3.8 | **No 6 h miRge command** — only the 24 h invocation survives, though `miR.Counts_C6_NE6_3v3_new.csv` is a required input | Stage 11 parameterises the recorded command | Medium |

## 4. Decisions only Seba and Juan Pablo can make

| # | Question | Consequence |
|---|---|---|
| 4.1 | **BSJ threshold: ≥2 or ≥3?** Methods say two, code says `readNumber > 2` | Changes the circRNA inventory, Figure 2B's 8,736, and every circRNA count in Figures 3 and S3. `--bsj-min` |
| 4.2 | **Was ComBat_seq used in the final circRNA DE?** Notebooks yes, final volcano script no, Methods silent | `--use-combat`. Note batch maps 1:1 to replicate, so dropping replicate 3 removes Batch3 |
| 4.3 | **Should CIRI-AS cover all samples or only SQ13-SQ16?** The original loop is SQ13-SQ16 only | Gates the four circRNA IntaRNA target FASTAs. `--all --yes-widen` |
| 4.4 | **Which physical library is "rep3"?** No file maps `Ctrl6_rep1`-style Salmon aliases to `H1_ln`/`SQ13_ln`/… | The Methods should ideally say which library was excluded and why |
| 4.5 | **lncRNA volcano: padj alone, or padj + \|log2FC\| ≥ 0.5?** The script says the former, the Methods the latter | Changes Figure S3A-B lncRNA |
| 4.6 | **Is "CIRCexplorer3" the right name** for CIRCexplorer2 + CLEAR? | Methods wording |

## 5. Manuscript fixes needed

| # | Fix |
|---|---|
| 5.1 | `Rattus_norvegicus.mRatBN7.2.115.gtf` → `Rattus_norvegicus.GRCr8.115.gtf` (×2) |
| 5.2 | circRNA reference "Ensembl rat rn7 genome" → GRCr8 |
| 5.3 | FEELnc window "200 bp to 100 kb" → "1 kb to 100 kb" |
| 5.4 | miRge flags "(`-lib` and `-rat`)" → `-lib`, `-on rat`, `-db mirbase` |
| 5.5 | **Regenerate Figure 3B's lncRNA panel** — published 145/275/328, correct 148/258/336 |
| 5.6 | Figure 4 and 5 legends cite the wrong supplementary tables |
| 5.7 | "Supplementary Table 10" → S8 |
| 5.8 | Remove the three Spanish editorial notes; fix "Differential.expression" |
| 5.9 | Fill `[FINAL GITHUB REPOSITORY URL]` ×2 and `[SEQUENCING ACCESSION(S)]` |
| 5.10 | Soften "no additional interaction tiers were used", or note the tiers are intermediate audit output only |

## 6. Repository loose ends

| # | Item |
|---|---|
| 6.1 | All `checksum` cells in `metadata/reference_manifest.tsv` are `TO_FILL` |
| 6.2 | Affinity sources (`.afdesign`/SVG) and the Cytoscape `.cys` or node coordinates, if pixel-level reproduction is wanted |
| 6.3 | Conda environments are templates; replace with `conda env export` once built |
| 6.4 | `sessionInfo()` from the machine that produces the deposited figures |
| 6.5 | No full pipeline run has been executed — the repository is verified structurally (8 checks) and functionally only for Figure 2B |

## 7. Explicitly out of scope

| Item | Why |
|---|---|
| The STAR branch of the metadata analysis | HISAT2 only, confirmed by S. Urquiza |
| The Babi Quarto notebooks | Logic ported to R01; originals stay in `SUZ2026/Babi analysis/` |
| Pixel-identical figures | Affinity assembly and manual Cytoscape layout sit between pipeline output and the manuscript |
| Figure content validation against the drafts | Needs the final figures; only Figure 2B has been numerically verified |
