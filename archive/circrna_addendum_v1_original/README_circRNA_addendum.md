# NE-NRVM circRNA repository addendum

This addendum describes the circRNA module to add to the existing private repository for the NE-induced NRVM hypertrophy manuscript.

## Purpose

The circRNA branch is being reintroduced into the manuscript/repository. It should be implemented as an independent module that does not disrupt the existing mRNA/lncRNA/miRNA/TF workflow.

## Proposed added modules

- `12_circRNA_identification_quantification/`
- `13_circRNA_miRNA_IntaRNA/`
- `14_circRNA_ceRNA_cross_reference/`
- `15_cross_ncRNA_networks/`

## Analytical scope

1. circRNA detection and quantification from long RNA-seq data using CIRCexplorer3 and CIRIquant.
2. Filtering of circRNAs for differential expression using BSJ read support.
3. circRNA sequence reconstruction with CIRI-AS.
4. circRNA–miRNA interaction prediction using IntaRNA.
5. Integration with miRNA DE, miRNA–mRNA target evidence, mRNA DE, and existing lncRNA–miRNA–mRNA triplets.
6. Export of strict circRNA–miRNA–mRNA triplets and cross-ncRNA convergence tables for downstream figures/supplementary material.

## Important reference consistency note

The main RNA-seq workflow uses Ensembl `Rattus_norvegicus.GRCr8.dna.primary_assembly` and `Rattus_norvegicus.mRatBN7.2.115.gtf`. The circRNA methods text currently says “Ensembl rat rn7 genome”. Before freezing the repository, confirm whether this refers to the same final rat reference build used in the rest of the manuscript or whether circRNA detection used a distinct rn7-specific reference.

## Script included here

- `scripts/24_FINAL_circRNA_lncRNA_miRNA_mRNA_cross_reference.R`

This script contains the final cross-reference logic but still has local Windows paths. It should be parameterized before being used in the automated repository workflow.

## Expected primary outputs

- `00_INPUT_manifest.csv`
- `03_HC_circRNA_miRNA_duplexes_locus_collapsed.csv`
- `08_CROSS_ncRNA_convergent_modules.csv`
- `11_ALL_putative_circRNA_miRNA_mRNA_relationships.csv`
- `12_STRICT_DE_coherent_circRNA_triplets.csv`
- `13_SHARED_with_lncRNA_triplet_axes.csv`
- `14_NEW_strict_circRNA_triplets_not_in_lncRNA_analysis.csv`
- `15_EXPRESSED_nonDE_circRNA_miRNA_targets.csv`
- `16_AUDIT_circRNA_target_filtering.csv`
- `17_CYTOSCAPE_cross_ncRNA_nodes.csv`
- `18_CYTOSCAPE_cross_ncRNA_edges.csv`
- `19_circRNA_temporal_direction_by_locus.csv`
- `20_circRNA_directional_overlap_counts.csv`
- `21_circRNA_temporal_class_counts.csv`
- `22_RUN_summary_counts.csv`
