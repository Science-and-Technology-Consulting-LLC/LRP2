# File: container/isocall/Dockerfile
# The Isocall image that ran before, plus the Isocall configuration, so the run needs no config file.

FROM docker.io/jtllab/isocall:1.3.0-nextflow-fix

LABEL maintainer="Science-and-Technology-Consulting-LLC"

COPY context/isocall_config.toml /opt/isocall/isocall_config.toml
USER root
RUN chmod -R a+rX /opt/isocall
