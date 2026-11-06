# Packages that Ubuntu's r-base does not ship. Bioconductor release follows R 4.3.3 (3.18).
options(timeout = 600)
install.packages(c("BiocManager", "optparse", "patchwork", "viridis"), repos = "https://cloud.r-project.org/")
BiocManager::install(c("edgeR", "DRIMSeq", "BiocParallel", "GenomicFeatures"), update = FALSE, ask = FALSE)
