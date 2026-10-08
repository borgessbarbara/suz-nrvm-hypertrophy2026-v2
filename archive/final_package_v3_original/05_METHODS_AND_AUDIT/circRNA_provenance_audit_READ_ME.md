# CircRNA provenance audit — read before parameterization

The newly recovered files substantially close the missing upstream circRNA branch, but they should **not be run unchanged** in the final repository.

## What is now recovered

- Three CIRCexplorer-family shell scripts cover H1-H4, SQ1-SQ8, and SQ13-SQ16.
- Three CIRIquant shell scripts cover the same 16 long-RNA samples.
- A BWA-MEM/CIRI2/CIRI-AS script provides the reconstruction command pattern.
- R notebooks document annotation, count-matrix construction, differential-expression experiments, and merging of CIRCexplorer/CIRIquant calls.
- Exact IntaRNA commands are now preserved for all four circRNA-miRNA direction/state combinations and all four lncRNA-miRNA direction/state combinations.

## Important discrepancies that must remain visible

1. **Reference files:** historical circRNA scripts include `GRCr8.114.gtf`, `rn8.refflat.nochr.txt`, `rn8.yml`, and rn8 index names, whereas the final project reference lock has been stated as GRCr8 primary assembly plus `Rattus_norvegicus.mRatBN7.2.115.gtf`. Do not silently replace the historical paths in the provenance files. The clean parameterized workflow should make the reference explicit and validation should determine which version reproduces the deposited circRNA results.

2. **BSJ filter:** the prose record says two or more BSJ reads, but `01_circexp3.Rmd` and `02_ciriquant.Rmd` use `readNumber > 2`. These are not equivalent.

3. **Historical DE thresholds:** the old notebooks classify circRNAs using nominal p < 0.05 and |log2FC| > 1 after ComBat_seq. The final manuscript/figure script uses nominal p < 0.05 and |log2FC| >= 0.5. Treat the notebook DE sections as historical unless confirmed as the exact source of final tables.

4. **Replicates:** the old notebooks build 4-control/4-NE matrices for each state. The current final circRNA volcano script uses three controls and three NE samples. This must be reconciled before the repository is frozen.

5. **CIRI-AS coverage:** the supplied CIRI2/CIRI-AS shell loop is written only for SQ13-SQ16. The generic command can be parameterized for all samples, but only after confirming that this reflects what was actually run.

6. **CIRCexplorer nomenclature:** the methods record calls the workflow CIRCexplorer3, while the recovered shell scripts explicitly invoke `CIRCexplorer2 parse` and `CIRCexplorer2 annotate` together with CLEAR `run3.py` and `circ_quant.py`. Preserve the actual commands and confirm the environment/package version.

7. **Final table lineage:** `03_tool_comparation.Rmd` clearly merges CIRCexplorer/CIRIquant calls and adds sequence/annotation/count information, but the exact step producing the later filenames `raw_count_matrix_all_6.txt`, `raw_count_matrix_all_24.txt`, `06_circRNA_sequences_de6_DE_complete_table.txt`, and `07_circRNA_sequences_de24_DE_complete_table.txt` is not fully explicit in the recovered files.

## Manual figure provenance

Some final multi-panel figures were assembled or cosmetically edited in Affinity. Therefore the analytical scripts should be expected to reproduce the underlying source plots/data, not necessarily the exact final page composition.

The Cytoscape network layouts were manually adjusted. The node/edge tables reproduce network membership and attributes; exact visual reproduction additionally requires the saved Cytoscape session or final node coordinates if available.
