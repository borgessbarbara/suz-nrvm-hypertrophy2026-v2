# circRNA methodology supplied for repository reconstruction

circRNA identification was conducted using two tools: CIRCexplorer3 v1.0.0 and CIRIquant v1.1, both referencing the Ensembl rat rn7/rn8-era genome resource as recorded in the historical workflow, following the approach used in Silva et al. 2025.

For the CIRCexplorer3 branch, sequencing reads were aligned to the rat genome using HISAT2 v2.2.1, and unmapped reads were aligned with STAR v2.7.11b to detect back-splicing junctions and extract readNumber values.

For the CIRIquant branch, trimmed clean reads were processed through CIRIquant; the historical auxiliary reconstruction script also shows BWA-MEM with `-T 19`, followed by CIRI2 and CIRI-AS.

Only circRNAs passing the BSJ read-support rule were carried forward for differential-expression analysis. The exact threshold requires final audit because the prose record states two or more BSJ reads, whereas the historical R notebooks implement `readNumber > 2`.

circRNA sequences were reconstructed using CIRI-AS. The final manuscript treats circRNA differential expression as exploratory and uses nominal p < 0.05 and |log2FC| >= 0.5.
