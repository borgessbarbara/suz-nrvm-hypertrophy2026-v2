# NE NRVM reproducibility starter package v1

This package is intended to be sent with the raw script archive (`NE_NRVM_all_scripts_raw.zip`) so a co-author can start preparing the conda/docker runtime environment and repository structure.

## What this package contains

- `metadata/longRNA_samples.tsv`: preliminary long RNA-seq sample sheet.
- `metadata/smallRNA_samples.tsv`: preliminary small RNA-seq/miRNA sample sheet.
- `metadata/fastq_manifest_template.tsv`: template to fill with exact FASTQ filenames, paths, file sizes and md5 checksums.
- `metadata/reference_manifest.tsv`: manifest of final reference resources and processed resources required by the project.
- `metadata/software_tool_inventory.tsv`: preliminary external tool list for conda/docker.
- `metadata/R_package_inventory_from_uploaded_scripts.tsv`: R packages detected from the uploaded scripts.
- `scripts_inventory/script_inventory_preliminary.tsv`: preliminary classification of uploaded scripts into final-use, extract-logic, or archive.
- `commands/make_fastq_manifest_checksums.sh`: helper command to compute checksums and file sizes from the server FASTQ folder.
- `commands/example_future_pipeline_call.sh`: proposed future interface for the final pipeline.

## Important decisions already made

1. The final manuscript references are:
   - Genome FASTA: `Rattus_norvegicus.GRCr8.dna.primary_assembly.fa`
   - GTF annotation: `Rattus_norvegicus.mRatBN7.2.115.gtf`

2. Older references appearing in the scripts are historical tests and should not be used for the final reproducible pipeline.

3. circRNA and old 8-figure branches should be archived but excluded from the final public workflow.

4. The public workflow should be Linux-first. Windows use should be through WSL2 or Docker, not native Windows paths.

5. miRge3.0 was actually run with `-db mirbase` according to the uploaded command notes. The manuscript text should be corrected if this is the final analysis.

## What the co-author should do first

1. Use `metadata/software_tool_inventory.tsv` and `metadata/R_package_inventory_from_uploaded_scripts.tsv` to start preparing the conda/docker environment.
2. Use `metadata/longRNA_samples.tsv` and `metadata/smallRNA_samples.tsv` to understand the expected sample structure.
3. Ask Sebastián to fill `metadata/fastq_manifest_template.tsv` with exact server paths, backup paths, file sizes and checksums.
4. Use `scripts_inventory/script_inventory_preliminary.tsv` to see which scripts are likely final and which are historical/archive.
5. Do not try to run all raw scripts as-is. The next step is to extract the final logic into clean modular scripts.

## What still needs to be provided

- Exact FASTQ filenames and paths from server/hard drive.
- md5 checksums for all FASTQ files.
- Exact local path or download source for TFLink All.
- Exact miRBase version/library used by miRge3.0.
- Final figure-input tables for Figures 3–7.
- R `sessionInfo()` and `conda env export` outputs from the original environment if available.
