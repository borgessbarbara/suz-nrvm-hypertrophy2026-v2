# syntax=docker/dockerfile:1
# =============================================================================
# One image carrying all five pipeline environments.
#
# The base is miniforge3, not miniconda or micromamba, for two reasons:
#   * it provides the `conda` binary, which scripts/bash/_activate_env.sh calls
#     directly (`conda info --base`, then `conda activate`). micromamba does
#     not, and every bash stage would fail to activate.
#   * it defaults to conda-forge and ships no `defaults` channel, matching the
#     channel set in envs/*.yml.
#
# Build (always target linux/amd64 - bioconda publishes linux-64, so an ARM
# host such as an Apple-silicon Mac must emulate):
#
#     docker build --platform linux/amd64 -t suz-nrvm:latest .
#
# This installs five conda environments including R and Bioconductor. Expect a
# long first build and an image of several GB. To build a subset:
#
#     docker build --platform linux/amd64 --build-arg ENVS="qc_align r" \
#                  -t suz-nrvm:qc .
#
# Run, mounting your data and writing results back to the host as yourself:
#
#     docker run --rm -it --platform linux/amd64 \
#       --user "$(id -u):$(id -g)" \
#       -v /path/to/fastq:/data/fastq:ro \
#       -v /path/to/references:/data/references:ro \
#       -v "$PWD/results:/opt/suz/results" \
#       suz-nrvm:latest bash
#
# Then inside the container, exactly as on a host:
#
#     bash scripts/commands/create_envs.sh --verify-only
#     FASTQ_DIR=/data/fastq bash scripts/bash/01_qc_raw.sh
#     conda activate env_r && bash scripts/bash/_run_r_analysis.sh
#
# NOT IN THIS IMAGE. CIRI2 (CIRI2.pl) and CIRI-AS (CIRI_AS_v1.2.pl) are
# standalone Perl scripts that are not redistributable through conda; mount them
# and pass --ciri2-path / --ciri-as-path to stage 15. The references, the FASTQs
# and the other external inputs in metadata/external_inputs.tsv are mounted too
# - see docs/OPEN_ISSUES.md section 7.4.
# =============================================================================

ARG BASE_IMAGE=condaforge/miniforge3:latest
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="suz-nrvm-hypertrophy2026" \
      org.opencontainers.image.description="Reproducibility pipeline for the NE-NRVM manuscript: mRNA, lncRNA, miRNA and circRNA integration in neonatal rat ventricular myocytes" \
      org.opencontainers.image.source="https://github.com/integrativebioinformatics/suz-nrvm-hypertrophy2026" \
      org.opencontainers.image.licenses="NOASSERTION"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# git is needed to install CLEAR from source; procps and the rest are what the
# stages and their tools expect to find on a normal Linux box.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
        procps \
        less \
        unzip \
        bzip2 \
        xz-utils \
    && rm -rf /var/lib/apt/lists/*

# Make `conda activate` work in non-interactive shells. BASH_ENV covers
# `bash -c` and anything the stages invoke; bash.bashrc covers `docker run -it`.
ENV BASH_ENV=/opt/conda/etc/profile.d/conda.sh
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/bash.bashrc

# bioconda requires strict channel priority, which a YAML file cannot express.
RUN conda config --system --set channel_priority strict \
    && conda config --system --set always_yes true

WORKDIR /opt/suz

# --- environments ------------------------------------------------------------
# Copied and built before the rest of the repository, so editing a stage script
# does not invalidate the (very expensive) environment layers.
COPY envs/ /opt/suz/envs/

ARG ENVS="qc_align mirge intarna circrna r"
RUN for e in ${ENVS}; do \
        echo "=== creating env_${e} ===" && \
        mamba env create -f "envs/environment_${e}.yml" || exit 1; \
    done \
    && mamba clean -afy \
    && find /opt/conda -follow -type f -name '*.a'    -delete \
    && find /opt/conda -follow -type f -name '*.pyc'  -delete

# --- CLEAR -------------------------------------------------------------------
# No bioconda package exists, so it is built from source into env_circrna. The
# commit is recorded because nothing else pins it (docs/OPEN_ISSUES.md 15.6).
#
# pip --no-deps: pysam, pybedtools, CIRCexplorer2, HISAT2 and StringTie all come
# from conda, and letting pip resolve them would build them a second time.
# pip is tried before `setup.py install`, which current setuptools has removed.
ARG CLEAR_REPO=https://github.com/YangLab/CLEAR
ARG CLEAR_REF=master
RUN if echo " ${ENVS} " | grep -q ' circrna '; then \
        . /opt/conda/etc/profile.d/conda.sh && conda activate env_circrna && \
        git clone "${CLEAR_REPO}" /opt/CLEAR && \
        cd /opt/CLEAR && git checkout --quiet "${CLEAR_REF}" && \
        git rev-parse HEAD > /opt/CLEAR/INSTALLED_COMMIT && \
        { pip install --no-deps --no-build-isolation . \
          || python ./setup.py install ; } && \
        echo "CLEAR installed at commit $(cat /opt/CLEAR/INSTALLED_COMMIT)"; \
    else \
        echo "env_circrna not built; skipping CLEAR"; \
    fi
ENV CLEAR_SRC=/opt/CLEAR

# --- the pipeline itself -----------------------------------------------------
COPY . /opt/suz

# Fail the build rather than ship a broken image: the repository's own checks,
# plus a binary-presence check inside every environment that was built.
RUN . /opt/conda/etc/profile.d/conda.sh \
    && for f in scripts/bash/*.sh scripts/checks/*.sh scripts/commands/*.sh suz_pipeline.bash; do \
           bash -n "$f" || exit 1; \
       done \
    && bash scripts/checks/check_biomart_pinning.sh \
    && bash scripts/checks/check_hardcoded_paths.sh \
    && bash scripts/commands/create_envs.sh --verify-only \
         --only "$(echo ${ENVS} | tr ' ' ',')" \
    && if echo " ${ENVS} " | grep -q ' r '; then \
           conda activate env_r && Rscript scripts/checks/check_io_graph.R; \
       fi

# Record what actually resolved, so the image carries its own version manifest.
RUN . /opt/conda/etc/profile.d/conda.sh \
    && for e in ${ENVS}; do \
           conda env export --no-builds -n "env_${e}" \
             > "/opt/suz/envs/environment_${e}.lock.yml"; \
       done \
    && echo "lock files written to /opt/suz/envs/"

# results/ is the only path the pipeline writes to by default. Pre-create it
# world-writable so `docker run --user $(id -u)` works without a bind mount.
RUN mkdir -p /opt/suz/results && chmod 0777 /opt/suz/results

CMD ["bash"]
