# Manuscript Methods vs the original scripts

Every point where the v14 Methods text and the scripts Seba actually ran
disagree, with the evidence for each and what the repository does about it.

Sources compared:
- `NE_NRVM_CDDIS_adapted_v14-SL.docx` — Methods and figure legends
- `provenance/longRNA_command_notebook/bash_scripts.txt` — the raw command log
- `provenance/circRNA_upstream/*` — the circRNA shell scripts and notebooks
- `provenance/mirge/miRge_3.txt`, `provenance/intarna/*`
- `archive/raw_scripts_original/`, `archive/final_package_v3_original/`

Line references are to those files.

---

## A. Contradictions

### A1. The GTF named in the Methods does not exist

| | |
|---|---|
| **Methods** | StringTie "guided by `Rattus_norvegicus.mRatBN7.2.115.gtf`"; the same filename again for FEELnc Module 1 |
| **Scripts** | `Rattus_norvegicus.GRCr8.114.gtf` — `bash_scripts.txt:293` (FEELnc filter), `:300` (classifier), `CIRI_2_run_all_seba.sh:6`, `vis_run_ciri2_seba.sh:2` |
| **Why** | Ensembl release 115 ships **GRCr8**, not mRatBN7.2. The assembly switched at release 114. The filename combines a release number and an assembly name from two different sources. |
| **Fix** | The GTF name is generated from the active release: `Rattus_norvegicus.GRCr8.${RELEASE}.gtf`. Release set in `metadata/reference_release.tsv`. |
| **Status** | Methods text still needs correcting. |

### A2. Release 114 vs 115

| | |
|---|---|
| **Methods** | Implies release 115 throughout |
| **Scripts** | Every executed command uses **114**: FEELnc (`:293`, `:300`), CIRI2 (`:6`), biomaRt `version = 114` in all three Babi notebooks. The Salmon work is dated `12/08/25`; release 115 was published 2025-09-03. |
| **Fix** | One switch (`metadata/reference_release.tsv`) drives the GTF, the FASTA and the biomaRt archive together. |
| **Status** | **Resolved by decision**: the analysis will be re-run at release 115, so 115 is the forward-looking lock and the Methods become correct. The frozen 114 tables in `data/frozen/` then need regenerating, and Figure 2B/2C and the Figure S2 lncRNA inventory should be re-checked against the published numbers. |

### A3. Genome FASTA: primary assembly, not toplevel — *the Methods were right*

| | |
|---|---|
| **Methods** | `Rattus_norvegicus.GRCr8.dna.primary_assembly` |
| **Scripts** | Agree. `bash_scripts.txt:27-28` "Ensembl primary assembly with all concate"; `:46` marks the toplevel index **`(NO USAR)`**; `:49-50` points at the primary-assembly index. circRNA branch uses `rat_primary_assembly_orderr8.fa`. |
| **Who was wrong** | This repository. It had been downloading `dna.toplevel.fa`, on my incorrect claim that Ensembl publishes no primary assembly for GRCr8 and that toplevel is equivalent. The Ensembl README says the opposite: primary-assembly files exist only when haplotype/patch regions do, and toplevel includes "N padded haplotype/patch regions". |
| **Fix** | Stage 00 builds the primary assembly — single `dna.primary_assembly.fa.gz` if published, otherwise the per-chromosome parts concatenated. It refuses to fall back to toplevel silently. |
| **Decision** | `nonchromosomal` is **excluded** (B. Borges, 8 Oct 2026): parts are `1..20 X Y MT`. Override with `PRIMARY_PARTS`. |

### A4. circRNA reference called "rn7"

| | |
|---|---|
| **Methods** | CIRCexplorer3 and CIRIquant "both referencing the Ensembl rat **rn7** genome" |
| **Scripts** | Every circRNA path is rn8/GRCr8: `rat_primary_assembly_orderr8.fa`, `rn8.refflat.nochr.txt`, `Rattus_norvegicus.GRCr8.114.gtf`, `STAR_rn8`, HISAT2 index `Rn8` |
| **Fix** | Documentation only; the stages take the reference as an argument. |
| **Status** | Methods text needs correcting. |

### A5. FEELnc classifier window

| | |
|---|---|
| **Methods** | "a 200 bp to 100 kb window (`--window=200 --maxwindow=100000`)" |
| **Scripts** | The final 6+24 h rn8 run, `bash_scripts.txt:300`: `--window=1000 --maxwindow=100000`, output `lncRNA_classes_CNE624_rn8_8_annotated_1kb_100kb.txt` — the exact file the pipeline consumes |
| **Why** | `--window=200` appears only in six earlier exploratory runs (`:199, :207, :209, :219, :227`), all with maxwindow 10 000 or 50 000, all producing differently named files that nothing reads. The Methods pair the window from the discarded runs with the maxwindow from the final one. |
| **Fix** | Stage 09 uses `1000 / 100000`, with the evidence quoted in the script header. |
| **Status** | Methods text needs correcting to "1 kb to 100 kb". |

### A6. Salmon index — *the Methods were right*

| | |
|---|---|
| **Methods** | "a genome-aware decoy/gentrome index" |
| **Scripts** | Agree. `bash_scripts.txt:232-235`: `salmon index -t gentrome_CNE624_8_H.fa -i transcripts_CNE624_H_index --decoys decoys.txt -p 6 -k 31` |
| **Who was wrong** | This repository — it built a transcriptome-only index, which gives different quantifications. |
| **Fix** | Stage 10 builds a gentrome (transcripts + genome) with a decoy list, same `--decoys` and `-k 31`. `BUILD_DECOY=0` restores the old behaviour and warns. |

### A7. BSJ read-support threshold

| | |
|---|---|
| **Methods** | "only circRNAs with **at least two** BSJ counts were carried forward" |
| **Scripts** | `readNumber > 2` — i.e. **at least three** — `01_circexp3.Rmd:91`, `02_ciriquant.Rmd:85` |
| **Impact** | Changes the circRNA inventory, so Figure 2B's 8,736 loci and every circRNA count in Figure 3 and S3. |
| **Fix** | `scripts/R/04b_circRNA_de_temporal.R --bsj-min`, logged on every run, defaulting to 3 (the code). |
| **Status** | **Unresolved — blocking.** One of the two must change. |

### A8. CIRCexplorer3 vs CIRCexplorer2 + CLEAR

| | |
|---|---|
| **Methods** | "CIRCexplorer3 v1.0.0" |
| **Scripts** | `CIRCexplorer2 parse` and `CIRCexplorer2 annotate`, with CLEAR's `run3.py` and `circ_quant.py`. The string `circexplorer3` appears only as a directory name. |
| **Fix** | The actual commands are preserved verbatim in stage 13; both packages are pinned in `envs/environment_circrna.yml`. |
| **Status** | Defensible as written if CLEAR is what "CIRCexplorer3" means, but the Methods should name what ran. |

### A9. CIRI-AS sample coverage

| | |
|---|---|
| **Methods** | Implies circRNA sequences were reconstructed for the analysis as a whole |
| **Scripts** | `CIRI_2_run_all_seba.sh:2` — `for m in {13..16}`, i.e. **SQ13–SQ16 only**. No CIRI2/CIRI-AS loop exists for H1–H4 or SQ1–SQ8. |
| **Fix** | Stage 15 defaults to those four samples and refuses `--all` without `--yes-widen`. |
| **Status** | **Unresolved.** This gates the IntaRNA circRNA target FASTAs. |

### A10. ComBat_seq

| | |
|---|---|
| **Methods** | Not mentioned |
| **Scripts** | The historical notebooks apply `ComBat_seq` before DE (`01_circexp3.Rmd:132`, `:219`); the final volcano script `25_` does not. |
| **Note** | Batch maps 1:1 to replicate index, so excluding replicate 3 removes Batch3 entirely. |
| **Fix** | `--use-combat` in R04b, **off** by default. |
| **Status** | **Unresolved** — confirm the final tables came from the no-ComBat route. |

### A11. lncRNA volcano fold-change cutoff

| | |
|---|---|
| **Methods** | lncRNAs: adjusted p < 0.05 **and \|log2FC\| ≥ 0.5** |
| **Scripts** | `07_Figure4_clean_volcanoes.R` sets `absfc_cutoff <- 0` — significance on padj alone |
| **Fix** | R16 left as ported (it produces a published panel; changing it would change the figure). The new R17 (mRNA/miRNA) applies the Methods thresholds: `--lfc-cut-mrna 1`, `--lfc-cut-mirna 0.5`. |
| **Status** | **Unresolved.** Figure S3's legend says the dashed lines mark "the thresholds used for classification", so panel and Methods should agree. |

### A12. "No additional interaction tiers were used"

| | |
|---|---|
| **Methods** | "Only interactions with MFE ≤ −83.7 kJ/mol (−20 kcal/mol) were retained; no additional interaction tiers were used." |
| **Scripts** | `12_..._ALL_RECORDS_FINAL.R` defines three tiers — `energy_high -20`, `energy_moderate -15`, `energy_exploratory -10` — and sets `minimum_pair_energy <- energy_exploratory`, so intermediate outputs retain pairs down to −10. |
| **Reality** | The **final** networks do use −20 (`HIGH_CONFIDENCE_ceRNA_triplets_E_le_minus20.csv`, and the MFE20 naming through R24/R25), so the figures match the Methods. The sentence is just not literally true of the script. |
| **Fix** | Documented; the tiers are retained because they produce audit outputs. |
| **Status** | Soften the sentence, or state that tiers exist only as intermediate audit output. |

### A13. miRge3.0 invocation

| | |
|---|---|
| **Methods** | "using the rat miRBase library (`-lib` and `-rat`)" |
| **Scripts** | `-lib <path> -on rat -db mirbase` — there is no `-rat` flag; the organism flag is `-on rat`, and the database flag `-db mirbase` is what makes it miRBase |
| **Fix** | Stage 11 uses the actual flags. |
| **Status** | Minor; correct the flag names in the Methods. |

### A14. Sample size vs the sample sheets

| | |
|---|---|
| **Methods / legends** | "N = 3 independent biological replicates per condition" |
| **Scripts** | Agree — 3 v 3 everywhere: `Integrated_Bulk_v3_6_24h.R:30-33` uses `rep1, rep2, rep4`; the miRNA matrices are `*_3v3_new.csv`; `25_*.R:284-291` drops `_3` |
| **But** | The sample sheets list 16 long-RNA and 16 small-RNA libraries, 4 per group. 12 of each were analysed. |
| **Fix** | `metadata/longRNA_samples.tsv` marks the four replicate-3 rows `include_final_analysis = no_excluded_from_final_DE`. |
| **Status** | Consistent. Open only as A-list item A4 in `OPEN_ITEMS.md`: *which physical library* "rep3" was is unrecorded. |

### A15. Parameters the Methods omit

Not contradictions, but the Methods are not sufficient to reproduce these:

| Step | In the scripts, absent from the Methods |
|---|---|
| StringTie | `-f 0.05 -m 200 --conservative -t -s 3.0` (`bash_scripts.txt:107`, `:152-158`) |
| FEELnc codpot | `--mode=intergenic`, `-b transcript_status=KNOWN` (`:195`, `:205`, `:217`, `:225`) |
| fastp | `-3` and `-w` thread count; the CIRI2 branch additionally uses `-b 140 -B 140 -l 140` |
| Salmon index | `-k 31` |
| BWA-MEM | `-t 20` alongside the documented `-T 19` |
| STAR | `--chimSegmentMin 10 --runThreadN 20` |

All are carried in the corresponding stage with the source line cited.

---

## B. Verified as matching — no action

Checked explicitly and in agreement between the Methods and the scripts:

| Item | Evidence |
|---|---|
| fastp `-f 10 -F 10 -3 --cut_tail_mean_quality 28 -q 30` | `bash_scripts.txt:66`; stage 02 |
| HISAT2 `--dta-cufflinks`, `-S` | `:78-90`; stage 04 |
| Salmon quant `-l A --validateMappings` | `:241-288`; stage 10 |
| Salmon v1.10.3 | `:280` |
| StringTie guided `-G`, then merge, then gffcompare | `:132-140`, `:170`; stages 06-07 |
| FEELnc three modules, `--monoex=1` | `:193-227`, `:293-300`; stage 09 |
| FEELnc maxwindow 100 000 | `:300` |
| Correlation: **Pearson**, \|rho\| ≥ 0.8, p < 0.05, VST | `cor(..., method="pearson")` in `Integrated_Bulk_v3_6_24h.R:1491`; `Hmisc::rcorr(..., type="pearson")` in the archived scripts; `rho_cutoff <- 0.80` |
| mRNA DE padj < 0.05, \|log2FC\| ≥ 1 | R02 |
| circRNA DE nominal p < 0.05, \|log2FC\| ≥ 0.5 | `25_*.R:26-27` |
| IntaRNA `--mode=M --outMode=C`, MFE ≤ −20 kcal/mol | `IntaRNA_commands_as_used.sh`; all 8 output names match the `Interact_*.csv` glob in R21 |
| miRWalk score ≥ 0.95, 3'UTR/5'UTR/CDS, miRDB support | 12 files in `metadata/external_inputs.tsv` |
| miRge adapters `AGATCGGAAGAGCACACGTCT` / `GTTCAGAGTTCTACAGTCCGACGATC` | `miRge_3.txt` |
| BWA-MEM `-T 19` | `CIRI_2_run_all_seba.sh:98` |
| GRN: TFLink, expressed > 0 in ≥ 1 replicate, mean normalised count > 10 | R05 |
| Tool versions: HISAT2 2.2.1, STAR 2.7.11b, CIRIquant 1.1, CIRI2 2.0.6, CIRI-AS 1.2, fastp 0.22.0, StringTie 3.0.0, gffcompare 0.12.9, samtools 1.21, FEELnc 0.2, IntaRNA 3.4.1, Cytoscape 3.10.3 | `metadata/software_tool_inventory.tsv`, each with its evidence |

---

## C. Paper-internal inconsistencies (not script-related)

| # | Issue |
|---|---|
| C1 | **Figure 3B's lncRNA panel is stale.** Published 145 / 275 / 328 (Early/Sustained/Late); the Results text and Figure S3C both give 148 / 258 / 336. Derived from S3C by summing the annotated and novel matrices, every number matches the text. mRNA, miRNA and circRNA reconcile everywhere. Regenerate with R10. |
| C2 | Figure legends cite Supplementary Tables 4-5 and 7-8 where the Results correctly cite S2/S3 and S5/S6. |
| C3 | "Supplementary Table 10" is cited; the workbook contains S1-S8 (should be S8). |
| C4 | Three Spanish editorial notes remain in the draft, plus the typo "Differential.expression". |
| C5 | `[FINAL GITHUB REPOSITORY URL]` ×2 and `[SEQUENCING ACCESSION(S)]` ×1 unfilled. |
| C6 | Merge-plan item D cites "Table 4", which does not exist in v14. Assumed to mean S4. |
