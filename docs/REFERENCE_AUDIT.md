# Reference audit

What the project uses, what the historical scripts used, and which of the
apparent conflicts are real.

## The lock

| Resource | File | Release |
|---|---|---|
| Genome FASTA | `Rattus_norvegicus.GRCr8.dna.toplevel.fa` | Ensembl 115 |
| GTF | `Rattus_norvegicus.GRCr8.115.gtf` | Ensembl 115 |
| biomaRt | `https://sep2025.archive.ensembl.org` | Ensembl 115 |

Downloaded by `scripts/bash/00_download_references.sh` from:

```
https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz
https://ftp.ensembl.org/pub/release-115/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.115.gtf.gz
```

## Which rat assembly each Ensembl release ships

Checked directly against the Ensembl FTP during this rebuild:

| Release | Rat GTF | Assembly |
|---|---|---|
| 112 | `Rattus_norvegicus.mRatBN7.2.112.gtf.gz` | mRatBN7.2 (rn7) |
| 113 | `Rattus_norvegicus.mRatBN7.2.113.gtf.gz` | mRatBN7.2 (rn7) |
| **114** | `Rattus_norvegicus.GRCr8.114.gtf.gz` | **GRCr8 (rn8)** |
| **115** | `Rattus_norvegicus.GRCr8.115.gtf.gz` | **GRCr8 (rn8)** |

The switch happened at release 114. This single fact resolves most of the
apparent reference conflict in the historical material.

## `sep2025.archive.ensembl.org` is release 115

Two independent confirmations:

1. The release-115 FTP directory is dated **2025-09-03**, and Ensembl names its
   archives after the release month.
2. The archive's own BioMart registry says so:

```
$ curl -s https://sep2025.archive.ensembl.org/biomart/martservice?type=registry
<MartRegistry>
  <MartURLLocation database="ensembl_mart_115" default="1"
                   displayName="Ensembl Genes 115" ... />
```

A live connection through `connect_ensembl()` reports the dataset version as
`GRCr8`, matching the GTF.

### Mind the trailing slash

`https://sep2025.archive.ensembl.org/` — with a trailing slash — **fails**,
because biomaRt appends `/biomart/martservice` and produces a double slash.
The error is the unhelpful `Error: Can't retrieve empty body.`
`suz_normalize_host()` in `scripts/R/lib/ensembl.R` strips it, so both forms work.

## Conflicts in the historical material, and their status

| Appears in | Reference | Status |
|---|---|---|
| Starter `reference_manifest.tsv`, circRNA addendum | `Rattus_norvegicus.mRatBN7.2.115.gtf` | **Does not exist.** Release 115 has no mRatBN7.2 build. A transcription error. Corrected to `GRCr8.115.gtf`. |
| `bash_scripts.txt` (FEELnc, ×4) | `Rattus_norvegicus.mRatBN7.2.112.gtf` | rn7-era historical test. Not used in the final pipeline (S. Urquiza, 2 Jul 2026). |
| `bash_scripts.txt` (StringTie/gffcompare, ×5) | `Rattus_norvegicus.mRatBN7.2.ensembl.gtf` | Unversioned rn7-era file. Historical. |
| `CIRI_2_run_all_seba.sh:6`, `vis_run_ciri2_seba.sh:2` | `Rattus_norvegicus.GRCr8.114.gtf` | **Same assembly**, one release earlier. Annotation-release gap only. Open item C3. |
| `Babi analysis` metadata notebook | biomaRt `version = 114` | **Same assembly**, one release earlier. This is what the frozen tables were built with. Open item C1. |
| CIRCexplorer/CIRIquant scripts | `rat_primary_assembly_orderr8.fa`, `rn8.refflat.nochr.txt`, index `Rn8`, `STAR_rn8` | "rn8" is GRCr8. Not supplied; verify by checksum. Open item A8. |
| circRNA addendum README | "Ensembl rat rn7 genome" | Contradicted by the scripts, which are all rn8/GRCr8. Methods text needs correcting. |

**Summary:** there is no assembly conflict in anything that reaches the final
pipeline. Every path that matters is GRCr8. What remains is a one-release
annotation gap (114 vs 115) in two places, both documented above.

## Why the toplevel FASTA rather than primary_assembly

Ensembl publishes a separate `dna.primary_assembly` file only when the toplevel
set contains haplotype or patch regions. For GRCr8 it does not, so the two are
equivalent and Ensembl ships only the toplevel. This was raised with
S. Urquiza on 14 Jul 2026 and accepted.

Note that `bash_scripts.txt` marks a *different* toplevel file as
`#Para Hisat2 rata y ensembl toplevel (NO USAR)` — that was the **rn7**
`Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa`, a different file from a
different assembly. It does not apply to the GRCr8 toplevel.

## Why release 115 and not 116

Release 116 was available by 14 Jul 2026. S. Urquiza chose to stay on 115:
*"Creo que sale mejor mantener la 115."* The pipeline is pinned accordingly; a
future change requires editing `SUZ_ENSEMBL_HOST_DEFAULT` in
`scripts/R/lib/ensembl.R` **and** the URLs in `00_download_references.sh`, and
recording it here.

## Still to fill

`metadata/reference_manifest.tsv` has `TO_FILL` in every `checksum` cell and
most `location_to_provide` cells. The large resources (genome FASTA, GTF,
TFLink, Salmon index) are manifest-only by design; their checksums are not.
