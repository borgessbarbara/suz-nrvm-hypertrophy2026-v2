# references/

Genome FASTA and GTF annotation. **Nothing in this directory is committed** —
the files are too large, and they are reproducible from Ensembl.

## What belongs here

| File | Source |
|---|---|
| `Rattus_norvegicus.GRCr8.dna.primary_assembly.fa` | Ensembl, concatenated from the per-chromosome files (`1..20 X Y MT`) |
| `Rattus_norvegicus.GRCr8.116.gtf` | Ensembl release 116 (114 and 115 selectable) |

Both are fetched by `scripts/bash/00_download_references.sh`, which reads the
active release from `metadata/reference_release.tsv`. Names, URLs and checksums
are tracked in `metadata/reference_manifest.tsv`.

**Primary assembly, not toplevel.** Toplevel carries N-padded haplotype and
patch regions; the primary assembly excludes them, which is what Ensembl's own
FASTA README recommends for alignment. `nonchromosomal` is excluded by
decision. Override the parts with `PRIMARY_PARTS`.

## Using your own files

```bash
./suz_pipeline.bash \
  --genome-fasta /path/to/genome.fa \
  --gtf /path/to/annotation.gtf \
  ...
```

The defaults above are used only when these flags are omitted.

## Also not committed

`TFLink_All.tsv` is needed by `scripts/R/05_tf_kegg_integration.R` and is
supplied separately. It and the other 28 inputs no stage can produce are listed
in `metadata/external_inputs.tsv`.
