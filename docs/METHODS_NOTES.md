# Methods notes

Points where the manuscript text and the code disagree, or where the code does
something a reader would not assume. Each needs either a correction to the
Methods or a decision.

## 1. miRge3.0 used miRBase, not miRGeneDB

`archive/raw_scripts_original/miRge_scripts/miRge_3.txt` shows `-db mirbase` on all three recorded
invocations. Confirmed by S. Urquiza, 2 Jul 2026: *"I confirmed that the
miRge3.0 command used `-db mirbase`, so I will later correct the manuscript
Methods text accordingly."* The exact miRBase release is still unknown
(open item D3).

Recorded parameters:

```
-lib <miRge3_Lib> -on rat -db mirbase
-a AGATCGGAAGAGCACACGTCT          3' adapter
-g GTTCAGAGTTCTACAGTCCGACGATC     5' adapter
-dex -mdt DESmetadata_24h.csv
```

## 2. CIRCexplorer2 + CLEAR, described as CIRCexplorer3

The methods say "CIRCexplorer3 v1.0.0". The recovered scripts invoke
`CIRCexplorer2 parse` and `CIRCexplorer2 annotate`, together with CLEAR's
`run3.py` and `circ_quant.py`. The string `circexplorer3` appears only as a
directory name.

This is not necessarily an error: CLEAR is the package that provides the
circexplorer3-era quantification workflow, and describing the combination as
"CIRCexplorer3" is defensible. But the Methods should name what was actually
run. `envs/environment_circrna.yml` pins both.

## 3. Sample size is n = 3, not n = 4

Every final analysis uses three replicates per group — replicate 3 is excluded
throughout. See open item A3 for the evidence. The sample sheets list 16 long-RNA
libraries; 12 were analysed.

The Methods must say n = 3 per group, and ideally say which library was excluded
and why. That second part is currently unanswerable (open item A4).

## 4. circRNA differential expression is exploratory

The final circRNA volcanoes use a **nominal** p < 0.05 with |log2FC| >= 0.5, not
an FDR threshold (`25_recreate_circRNA_volcanoes_with_all_NS.R:26-27`). The
phase-aware ceRNA script states explicitly that "circRNA trajectories remain
exploratory". Figure 3C circRNA and Figure S3A-B circRNA should be described
that way.

By contrast the mRNA, lncRNA and miRNA volcanoes use `padj`.

## 5. The lncRNA volcanoes have no fold-change cutoff

`07_Figure4_clean_volcanoes.R` sets `absfc_cutoff <- 0`. Significance for the
lncRNA panels is adjusted p alone. The figure legend should say so, because a
reader will assume a fold-change threshold was applied.

## 6. BSJ support: at least 3 reads, not 2

See open item A1. Whichever rule is adopted, the Methods and the code must agree.

## 7. "rn8" means GRCr8

Used as a nickname throughout the historical scripts and directory names. The
circRNA addendum README calls it "Ensembl rat rn7 genome", which is wrong for
every path the scripts actually use. See `docs/REFERENCE_AUDIT.md`.

## 8. MSTRG/StringTie transcripts are treated as novel

`37c_..._MSTRG_FIXED.R` (ported as R12) reclassifies MSTRG-prefixed entries as
novel. `15_lncRNA_characterization_from_GTF.R` (R07) infers the same from the
`^MSTRG` / `^ENSRNOG` prefixes, and its own summary adds a caveat worth
reusing in the manuscript:

> Because this is a transcriptome reconstruction framework, novel MSTRG entries
> are best described conservatively as novel reconstructed lncRNA
> loci/transcripts rather than definitively novel lncRNA genes unless
> orthogonal evidence supports locus-level novelty.

## 9. The FEELnc `isBest == 1` filter does nothing here

Every row of the shipped classifier output already has `isBest == 1`. Do not
describe it as a filtering step. See open item B5.

## 10. IntaRNA pairing is deliberately anti-correlated

All eight runs pair UP ncRNAs against DOWN miRNAs and vice versa, because a
sponging relationship predicts opposite directions. Downstream the manuscript
retains duplexes with MFE <= -20 kcal/mol. All eight commands used
`--mode=M --outMode=C`.

## 11. The circRNA sequence files are not FASTA

`03_tool_comparation.Rmd` writes `circRNA_sequences_de{6,24}_2.fa` with
`write.table(sep="\t", col.names=TRUE)`. They are tab-separated tables with a
`.fa` extension and are read back with `read.table(header=TRUE)`. The rebuilt
stage names them `.tsv`.
