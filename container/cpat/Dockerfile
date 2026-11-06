# File: container/cpat/Dockerfile
# CPAT 3.0.4 with the human and mouse hexamer tables and logit models, so the run needs no reference files.

FROM ubuntu:24.04

LABEL maintainer="Science-and-Technology-Consulting-LLC"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential ca-certificates procps zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
        python3.12 python3.12-dev python3.12-venv python3-pip \
        r-base-core && \
    rm -rf /var/lib/apt/lists/* && \
    ln -sf /usr/bin/python3.12 /usr/local/bin/python && \
    ln -sf /usr/bin/python3.12 /usr/local/bin/python3

# same pins as the image that ran before
RUN python3 -m pip install --no-cache-dir --break-system-packages \
        numpy==1.26.4 pandas==2.2.3 scipy==1.11.4 biopython==1.81 pysam==0.22.1 cpat==3.0.4

# the build fails unless CPAT is the version we expect
RUN python -c "import importlib.metadata as m; v = m.version('CPAT'); print('CPAT', v); assert v == '3.0.4', v" && \
    cpat.py --version

# reference files; modules read them at /opt/cpat/
COPY context/Human_Hexamer.tsv context/Human_logitModel.RData context/Mouse_Hexamer.tsv context/Mouse_logitModel.RData /opt/cpat/
RUN chmod -R a+rX /opt/cpat

WORKDIR /data
CMD ["cpat.py", "--help"]
