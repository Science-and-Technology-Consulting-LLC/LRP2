# File: container/sqanti3_protein/Dockerfile
# SQANTI3 v5.5.4 with our protein script, already patched for the gtfToGenePred binary that v5.5.4 ships.
# The module used to copy the script into each task and run sed on it; the patch is now done once, here.

FROM docker.io/anaconesalab/sqanti3:v5.5.4

LABEL maintainer="Science-and-Technology-Consulting-LLC"

# the base image runs as "user"; the patch needs root to write into /opt
USER root

COPY context/04_sqanti3_protein.py /opt/sqanti3_protein/04_sqanti3_protein.py

# v5.5.4 has gtfToGenePred-linux-x86_64 instead of gtfToGenePred; the build fails if the patch does not apply
RUN sed -i 's|GTF2GENEPRED_PROG = os.path.join(sqanti_path, "src", "utilities", "gtfToGenePred")|GTF2GENEPRED_PROG = os.path.join(sqanti_path, "src", "utilities", "gtfToGenePred-linux-x86_64")|g' /opt/sqanti3_protein/04_sqanti3_protein.py && \
    grep -q 'gtfToGenePred-linux-x86_64' /opt/sqanti3_protein/04_sqanti3_protein.py && \
    chmod -R a+rX /opt/sqanti3_protein

USER user
