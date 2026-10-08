# Reference audit

## The headline

**Every executed command in this project used Ensembl release 114. No command
anywhere uses release 115.**

The "release 115" lock that appears in the starter package, the circRNA
addendum, the v3 audit documents and the manuscript Methods is traceable to a
single source — the 7 Jul starter package's `reference_manifest.tsv` — and it
names a file that does not exist.

This section is the evidence. The release the pipeline actually runs with is
set in `metadata/reference_release.tsv`, and currently reads **115**, which is
what the project was instructed to lock to. Changing it is a one-line edit that
moves the GTF, the FASTA and the biomaRt archive together. See "What to do"
at the end.

---

## Which rat assembly each Ensembl release ships

Checked directly against the Ensembl FTP:

| Release | Rat GTF | Assembly | biomaRt archive |
|---|---|---|---|
| 112 | `Rattus_norvegicus.mRatBN7.2.112.gtf.gz` | mRatBN7.2 (rn7) | — |
| 113 | `Rattus_norvegicus.mRatBN7.2.113.gtf.gz` | mRatBN7.2 (rn7) | — |
| **114** | `Rattus_norvegicus.GRCr8.114.gtf.gz` | **GRCr8 (rn8)** | `may2025.archive.ensembl.org` → `ensembl_mart_114` |
| **115** | `Rattus_norvegicus.GRCr8.115.gtf.gz` | **GRCr8 (rn8)** | `sep2025.archive.ensembl.org` → `ensembl_mart_115` |

Both archive hosts verified live via `/biomart/martservice?type=registry`.
The rat assembly switched from mRatBN7.2 to GRCr8 at release 114.

**`Rattus_norvegicus.mRatBN7.2.115.gtf` has never existed.** Release 115 ships
GRCr8. That filename is the tell: wherever it appears, the release number and
the assembly name have been combined from two different sources.

---

## The evidence for release 114

Five independent lines, all pointing the same way:

### 1. The final FEELnc run

`provenance/longRNA_command_notebook/bash_scripts.txt`, under the heading
*"Para 6 y 24 hrs"* — the 6+24 h rn8 run whose output filename this pipeline
consumes:

```
:293  FEELnc_filter.pl     -i merged_CNE624_8_rn8_H.annotated.gtf
                           -a .../Rattus_norvegicus.GRCr8.114.gtf ...
:300  FEELnc_classifier.pl -i cand_lncRNAs_CNE624_rn8_8_annotated.gtf.lncRNA.gtf
                           -a .../Rattus_norvegicus.GRCr8.114.gtf
                           --window=1000 --maxwindow=100000
                           > lncRNA_classes_CNE624_rn8_8_annotated_1kb_100kb.txt
```

### 2. The circRNA branch

```
CIRI_2_run_all_seba.sh:6    reference_annotation=.../Rattus_norvegicus.GRCr8.114.gtf
vis_run_ciri2_seba.sh:2     reference_annotation=.../Rattus_norvegicus.GRCr8.114.gtf
```

### 3. biomaRt in all three metadata notebooks

```
Babi analysis/2_nrvm_transcriptome_metadata_handling.qmd:535   version = 114
Babi analysis/nrvm_transcriptome_metadata_handling.qmd:262     version=114
Babi analysis/.Rhistory:415                                    version = 114
```

### 4. The frozen metadata tables

`data/frozen/hisat2_6_24_meta_ann.csv.gz` was produced by that notebook, at 114.

### 5. Figure 2B reproduces exactly from them

`scripts/R/06_transcriptome_overview.R` run against the frozen 114 tables gives
21,981 protein-coding / 16,149 annotated lncRNA / 840 novel lncRNA / 1,778
snoRNA / 1,517 snRNA / 204 miRNA — every number in the published panel. A
different annotation release would not do that.

### Chronology corroborates it

`bash_scripts.txt` dates the Salmon index work **12/08/25** (12 August 2025).
Release 115 was published **2025-09-03**. The long-RNA quantification predates
release 115 by about three weeks.

---

## Where "115" actually comes from

Every occurrence is prose or a manifest. None is an executed command.

| File | Form |
|---|---|
| `NE_NRVM_coauthor_reproducibility_starter_v1/metadata/reference_manifest.tsv:3` | `mRatBN7.2.115`, "release 115" — **the origin** |
| …`/README.md:26`, `/README.html:395`, `/docs/message_to_coauthor_draft.txt:11`, `/commands/example_future_pipeline_call.sh:10` | `mRatBN7.2.115` |
| `NE_NRVM_circRNA_repo_addendum_v1/.../circRNA_reference_manifest_addendum.tsv:3,4`, `README_circRNA_addendum.md:27` | `mRatBN7.2.115` |
| `NE_NRVM_FINAL_repository_update_v3/04_METADATA/reference_resource_audit.tsv:3-4`, `remaining_open_items.tsv:2`, `05_METHODS_AND_AUDIT/circRNA_provenance_audit_READ_ME.md:15`, `06_NOTES_FOR_COLLABORATOR/handoff_summary.txt:10` | quoting the starter package |
| Manuscript v14 Methods | `Rattus_norvegicus.mRatBN7.2.115.gtf`, twice |

The v3 audit documents did not invent it — they correctly flagged the conflict
("historical circRNA scripts include `GRCr8.114.gtf` … whereas the final
project reference lock has been stated as … `mRatBN7.2.115.gtf`") and asked for
validation. This is that validation, and it resolves in favour of 114.

S. Urquiza's recollection (WhatsApp, 14 Jul 2026: *"el GTF fue la versión 115"*,
and *"creo que sale mejor mantener la 115"*) is the only other support, and was
given a year after the commands were run.

---

## Other reference facts

| Appears in | Reference | Status |
|---|---|---|
| `bash_scripts.txt` (19×) | `Rattus_norvegicus.mRatBN7.2.112.gtf` | rn7-era historical tests. Superseded by the rn8 runs. |
| `bash_scripts.txt` (5×) | `mRatBN7.2.112.chr.gtf` | featureCounts, historical. |
| `bash_scripts.txt:107,132-140` | `mRatBN7.2.ensembl.gtf` (unversioned) | Early StringTie/gffcompare. Historical. |
| CIRCexplorer scripts | `rat_primary_assembly_orderr8.fa`, `rn8.refflat.nochr.txt`, `Rn8`, `STAR_rn8` | "rn8" = GRCr8. Not supplied; verify by checksum. |
| circRNA addendum README | "Ensembl rat rn7 genome" | **Wrong.** Every circRNA path is rn8/GRCr8. |
| Manuscript Methods | `GRCr8.dna.primary_assembly` | **Correct — use this, not toplevel.** See below. |

### Primary assembly vs toplevel — the Methods are right

An earlier version of this document claimed Ensembl publishes no separate
`dna.primary_assembly` for GRCr8 and that toplevel is equivalent. **Both claims
were wrong.** Ensembl publishes `Rattus_norvegicus.GRCr8.dna.primary_assembly.*`,
split per chromosome.

The Ensembl FASTA README defines the two:

> **TOPLEVEL** — "...includes chromosomes, regions not assembled into
> chromosomes and **N padded haplotype/patch regions**."
>
> **PRIMARY ASSEMBLY** — "...all toplevel sequence regions **excluding
> haplotypes and patches**. This file is best used for performing sequence
> similarity searches where patch and haplotype sequences would confuse
> analysis. If the primary assembly file is not present, that indicates that
> there are no haplotype/patch regions, and the 'toplevel' file is equivalent."

Ensembl publishes `dna.primary_assembly.*` for GRCr8, so primary assembly and
toplevel are **not** equivalent here — toplevel carries N-padded patch and
haplotype sequence, which is exactly what the README warns against for
alignment.

The raw command notebook shows this was understood at the time:

```
bash_scripts.txt:27-28   #Ensembl primary assembly with all concate
                         .../genome_primary_RN7.fa
bash_scripts.txt:46      #Para Hisat2 rata y ensembl toplevel (NO USAR)
bash_scripts.txt:49-50   #Para Hisat2 rata y ensembl primary assembly
                         .../index/hisat2/Rat/primary_assembly/Rn72
```

The toplevel index is explicitly marked **"NO USAR"** — do not use — and the
genome in use is the **concatenated** primary assembly. The circRNA branch is
consistent: `rat_primary_assembly_orderr8.fa`.

An earlier version of this document dismissed the `(NO USAR)` line as applying
only to the rn7 toplevel and therefore not carrying over to GRCr8. That was a
misreading: the principle is about toplevel versus primary assembly, and it
applies to any assembly where the two differ.

`scripts/bash/00_download_references.sh` now builds the primary assembly,
preferring a single `dna.primary_assembly.fa.gz` and otherwise downloading the
per-chromosome parts and concatenating them in the order
`1..20 X Y MT nonchromosomal` (override with `PRIMARY_PARTS`). It refuses to
fall back to toplevel silently.

**Open:** the exact concatenation order Seba used is not recorded — the circRNA
filename `rat_primary_assembly_orderr8.fa` implies a specific one. Order affects
BAM header order and some tools' output, though not alignment itself. Confirm
with S. Urquiza.

---

## How the release is set

One file: **`metadata/reference_release.tsv`**. Its `active` column drives

- `scripts/bash/00_download_references.sh` — the FTP URLs for the GTF and FASTA
- `scripts/bash/06`, `07`, `09` and `suz_pipeline.bash` — the default GTF path
- `scripts/R/lib/ensembl.R` — the biomaRt archive host

so the GTF and the biomaRt annotation can never drift apart. Override for a
single run with `ENSEMBL_RELEASE=114`.

```bash
# verify which release the pipeline will use
Rscript scripts/R/01_transcriptome_metadata_annotation.R --dry-run-biomart

# try the other one without editing anything
ENSEMBL_RELEASE=114 Rscript scripts/R/01_transcriptome_metadata_annotation.R --dry-run-biomart
```

---

## What to do

**Recommended:** set `active = yes` on release **114** in
`metadata/reference_release.tsv`, and correct the manuscript Methods to
`Rattus_norvegicus.GRCr8.114.gtf` + `Rattus_norvegicus.GRCr8.dna.toplevel.fa`.

That makes the repository, the frozen metadata tables, the FEELnc run, the
circRNA branch and the published Figure 2B all consistent with each other and
with what was actually executed — and it removes the impossible filename from
the Methods.

**Currently active: 115**, because that is the lock the project was instructed
to adopt. Keeping it is defensible only as a deliberate decision to re-run the
annotation at a newer release, and it would require regenerating
`data/frozen/` and re-checking Figure 2B's counts, Figure 2C, and the lncRNA
inventory in Figure S2 — because those numbers currently come from 114.

**This is S. Urquiza's call, not an edit to make silently.** It is the top item
in `docs/OPEN_ITEMS.md`.
