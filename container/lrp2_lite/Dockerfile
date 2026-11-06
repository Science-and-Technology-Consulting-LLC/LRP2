# File: container/lrp2_lite/Dockerfile
# Our R and Python scripts and the libraries they load, built the same way as nsforest, scsilhouette and cellxgene-harvester:
# micromamba base, packages from an environment file, a build-time check, no ENTRYPOINT.

FROM mambaorg/micromamba:1.5.6

LABEL maintainer="Science-and-Technology-Consulting-LLC"

USER root:root

# procps: Nextflow needs ps to collect task metrics
RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --chown=mambauser:mambauser context/lrp2_lite.yml /app/env.yaml

USER mambauser:mambauser

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:$PATH \
    DEBIAN_FRONTEND=noninteractive \
    R_LIBS_USER=

RUN micromamba install -y -n base -f /app/env.yaml && \
    micromamba clean --all --yes

# the build fails unless every package the scripts load is there and edgeR and DRIMSeq are the versions that ran before
RUN Rscript -e 'pk <- c("GenomicRanges","tidyverse","rtracklayer","Biostrings","GenomicFeatures","data.table","magrittr","optparse","igraph","Matrix","matrixStats","edgeR","DRIMSeq","BiocParallel","viridis","patchwork"); for (p in pk) { suppressMessages(library(p, character.only = TRUE)); cat(p, as.character(packageVersion(p)), "\n") }; stopifnot(as.character(packageVersion("edgeR")) == "4.0.16", as.character(packageVersion("DRIMSeq")) == "1.30.0")' && \
    python -c "import numpy, pandas, scipy, Bio; print('numpy', numpy.__version__, 'pandas', pandas.__version__, 'scipy', scipy.__version__, 'biopython', Bio.__version__)"

# the scripts; modules call them at /opt/lrp2_lite/<script>
COPY --chown=mambauser:mambauser context/src/lrp2_lite/ /opt/lrp2_lite/

WORKDIR /data

# No ENTRYPOINT: Lifebit (AWS Batch) starts a task as "bash -o pipefail -c ...", and an entrypoint would receive that as arguments.
CMD ["Rscript", "--version"]
