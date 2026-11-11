#' review.R
#'
#' Helpers for review.ipynb, the one-page summary of a comparison. They only show what the other
#' notebooks saved (figures/ and tables/) and count rows of the pipeline's edgeR tables. They decide nothing.

suppressPackageStartupMessages(library(data.table))

# Show figures/<comparison>_<name>.png in the notebook, or say which notebook makes it when it is not there.
#   made_by  the notebook that saves the figure (named in the message when the file is missing)
show_fig <- function(name, made_by, dir = "figures") {
  f <- paste0(figure_stem(name, dir), ".png")
  if (file.exists(f)) {
    if (requireNamespace("IRdisplay", quietly = TRUE)) IRdisplay::display_png(file = f) else message("saved at ", f)
  } else {
    message("Not found: ", f, ". Run ", made_by, " first.")
  }
  invisible(f)
}

# Read tables/<comparison>_<name>.tsv (as saved by the other notebooks); NULL with a message when it is missing.
read_saved_table <- function(name, made_by, dir = "tables") {
  f <- paste0(figure_stem(name, dir), ".tsv")
  if (!file.exists(f)) { message("Not found: ", f, ". Run ", made_by, " first."); return(NULL) }
  fread(f, sep = "\t", data.table = FALSE)
}

# For genes, transcript isoforms and ORFs: how many rows edgeR tested and how many it found different.
#   alpha  FDR cut-off   lfc  fold-change cut-off (log2) for the "and at least this large" columns
#   up_in  name of the group edgeR called the second group (logFC above 0 means higher there)
#   down_in  name of the other group
de_summary <- function(loc, comparison, alpha = 0.05, lfc = 1, up_in = "case", down_in = "control") {
  files <- c(genes = file.path(loc$de_gene, paste0(comparison, "_DGE_edgeR_results.txt")),
             isoforms = file.path(loc$de_transcript, paste0(comparison, "_DTE_edgeR_results.txt")),
             ORFs = file.path(loc$de_orf, paste0(comparison, "_DE_ORF_edgeR_results.txt")))
  rows <- lapply(names(files), function(l) {
    if (!file.exists(files[[l]])) { message("Not found: ", files[[l]]); return(NULL) }
    r <- fread(files[[l]], sep = "\t", select = c("logFC", "FDR"), data.table = FALSE)
    r <- r[!is.na(r$FDR), ]
    sig <- r$FDR < alpha
    out <- data.frame(level = l, tested = nrow(r),
                      significant = sum(sig),
                      up = sum(sig & r$logFC > 0), down = sum(sig & r$logFC < 0),
                      significant_and_large = sum(sig & abs(r$logFC) >= lfc))
    names(out)[4:5] <- c(paste0("higher_in_", up_in), paste0("higher_in_", down_in))
    names(out)[6] <- paste0("significant_and_abs_log2FC_", lfc, "_or_more")
    out
  })
  do.call(rbind, rows)
}
