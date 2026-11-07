# File: container/fragpipe/Dockerfile
# FragPipe 24.0 with the default workflow files and the licence texts the setup step shows.
# The module used to read these from the repository's assets/ folder; they are inside the image now.

FROM docker.io/fcyucn/fragpipe:24.0

LABEL maintainer="Science-and-Technology-Consulting-LLC"

COPY context/DIA_SpecLib_Quant.workflow context/LFQ-MBR.workflow context/licensing_information.txt context/msfragger_licensing_information.txt /opt/fragpipe/

# the build fails if a file did not arrive
RUN test -s /opt/fragpipe/DIA_SpecLib_Quant.workflow && test -s /opt/fragpipe/LFQ-MBR.workflow && \
    test -s /opt/fragpipe/licensing_information.txt && test -s /opt/fragpipe/msfragger_licensing_information.txt && \
    chmod -R a+rX /opt/fragpipe
