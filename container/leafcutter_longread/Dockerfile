# File: container/leafcutter_longread/Dockerfile
# LeafCutter differential splicing for long reads: the leafcutter Python package (pyro, torch) and our scripts.
# The leafcutter package is copied from context/leafcutter and installed onto PYTHONPATH.

FROM ubuntu:24.04

LABEL maintainer="Science-and-Technology-Consulting-LLC"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    R_LIBS_USER=

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential ca-certificates procps zlib1g-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
        python3.12 python3.12-dev python3.12-venv python3-pip \
        r-base r-base-dev \
        r-cran-tidyverse r-cran-data.table r-cran-magrittr r-cran-igraph r-cran-matrix r-cran-optparse \
        r-bioc-rtracklayer && \
    rm -rf /var/lib/apt/lists/* && \
    ln -sf /usr/bin/python3.12 /usr/local/bin/python && \
    ln -sf /usr/bin/python3.12 /usr/local/bin/python3

COPY context/requirements.txt /tmp/requirements.txt
RUN python3 -m pip install --no-cache-dir --break-system-packages -r /tmp/requirements.txt \
        --extra-index-url https://download.pytorch.org/whl/cpu && \
    rm /tmp/requirements.txt

WORKDIR /opt/leafcutter_longread
COPY context/leafcutter/ /opt/leafcutter_longread/leafcutter/
COPY context/src/leafcutter_longread/ /opt/leafcutter_longread/
ENV PYTHONPATH=/opt/leafcutter_longread
RUN chmod -R a+rX /opt/leafcutter_longread

# the build fails unless the packages the scripts import are there
RUN python -c "import numpy, pandas, scipy, sklearn, torch, pyro, tqdm, rpy2, leafcutter; from leafcutter.differential_splicing.differential_splicing import differential_splicing; print('torch', torch.__version__, 'pyro', pyro.__version__)" && \
    Rscript -e 'for (p in c("data.table","tidyverse","magrittr","rtracklayer","Matrix","igraph","optparse")) { suppressMessages(library(p, character.only = TRUE)); cat(p, as.character(packageVersion(p)), "\n") }'

WORKDIR /data
CMD ["python", "/opt/leafcutter_longread/leafcutter_ds.py", "--help"]
