# File: container/metamorpheus/Dockerfile
# MetaMorpheus 1.1.14 with its settings, Mods and CustomAminoAcids folders writable.
# MetaMorpheus writes into /MetaMorpheus at run time; the module used to copy those folders out of the
# image and bind-mount them back. Making them writable here removes that step.

FROM docker.io/smithchemwisc/metamorpheus:1.1.14

LABEL maintainer="Science-and-Technology-Consulting-LLC"

USER root

RUN touch /MetaMorpheus/settings.toml && \
    chmod -R a+rwX /MetaMorpheus/Mods /MetaMorpheus/CustomAminoAcids /MetaMorpheus/settings.toml && \
    test -f /MetaMorpheus/Mods/ProteaseMods.txt && test -f /MetaMorpheus/CMD.dll
