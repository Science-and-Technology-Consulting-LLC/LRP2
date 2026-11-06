# File: container/lrp2_lite/Dockerfile
# Our R and Python scripts and the libraries they load. Versions match the jtllab/lrp2-lite:1.4 image that ran before:
# Ubuntu 24.04, R 4.3.3 with Bioconductor 3.18, Python 3.12.

FROM ubuntu:24.04

LABEL maintainer="Science-and-Technology-Consulting-LLC"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    R_LIBS_USER=

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential ca-certificates curl git procps \
        libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev libbz2-dev liblzma-dev \
        python3.12 python3.12-dev python3.12-venv python3-pip \
        r-base r-base-dev \
        r-cran-tidyverse r-cran-data.table r-cran-magrittr r-cran-igraph r-cran-matrix r-cran-matrixstats \
        r-bioc-genomicranges r-bioc-rtracklayer r-bioc-biostrings && \
    rm -rf /var/lib/apt/lists/* && \
    ln -sf /usr/bin/python3.12 /usr/local/bin/python && \
    ln -sf /usr/bin/python3.12 /usr/local/bin/python3

COPY context/r_packages.R /tmp/r_packages.R
RUN Rscript /tmp/r_packages.R && rm /tmp/r_packages.R

COPY context/requirements.txt /tmp/requirements.txt
RUN python3 -m pip install --no-cache-dir --break-system-packages -r /tmp/requirements.txt && rm /tmp/requirements.txt

# the build fails unless every package the scripts load is there, and the log records the versions
RUN Rscript -e 'pk <- c("GenomicRanges","tidyverse","rtracklayer","Biostrings","data.table","magrittr","optparse","igraph","Matrix","matrixStats","edgeR","DRIMSeq","BiocParallel","viridis","patchwork"); for (p in pk) { suppressMessages(library(p, character.only = TRUE)); cat(p, as.character(packageVersion(p)), "\n") }; stopifnot(as.character(packageVersion("edgeR")) == "4.0.16", as.character(packageVersion("DRIMSeq")) == "1.30.0")' && \
    python -c "import numpy, pandas, scipy, Bio; print('numpy', numpy.__version__, 'pandas', pandas.__version__, 'scipy', scipy.__version__, 'biopython', Bio.__version__)"

WORKDIR /opt/lrp2_lite
COPY context/src/lrp2_lite/ /opt/lrp2_lite/
RUN chmod -R a+rX /opt/lrp2_lite

WORKDIR /data
CMD ["Rscript", "--version"]
