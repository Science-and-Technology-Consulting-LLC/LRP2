#' volcano.R
#'
#' Volcano plot of an edgeR results table (logFC against -log10 FDR).
#' Used for the gene (DGE) and transcript isoform (DTE) tables written by S4.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# Read one edgeR results table written by the pipeline (tab separated).
read_edger_results <- function(path) {
  res <- fread(path, sep = "\t", data.table = FALSE)
  stopifnot(all(c("logFC", "FDR") %in% names(res)))
  res
}

# Colours of the novelty classes in the volcano plot (changed points only; unchanged points stay grey).
volcano_novelty_colours <- c(known = "#0072B2", novel = "#CC79A7", "no ORF" = "#E69F00")

# Volcano plot.
#   res        data frame with logFC and FDR
#   label_col  column that names each point (gene_name, transcript_name)
#   fdr, lfc   cut-offs that colour a point as changed
#   n_label    how many of the most significant changed points get a name
#   novelty    (optional) one entry per row of `res`: "known", "novel" or "no ORF" (see novelty_of in the notebooks).
#              Changed points are then coloured by it instead of by direction, and the legend counts them.
volcano_plot <- function(res, label_col, fdr = 0.05, lfc = 1, n_label = 10, title = NULL,
                         group_names = c("down", "up"), novelty = NULL) {
  stopifnot(label_col %in% names(res))
  res$neglog10fdr <- -log10(res$FDR)
  res$change <- "not changed"
  res$change[res$FDR < fdr & res$logFC >=  lfc] <- group_names[2]
  res$change[res$FDR < fdr & res$logFC <= -lfc] <- group_names[1]
  res$label <- ""
  changed <- which(res$change != "not changed")
  top <- changed[order(res$FDR[changed])]
  top <- head(top, n_label)
  lab <- res[[label_col]]
  lab <- ifelse(is.na(lab) | lab == "", res[[1]], lab)   # novel isoforms have no name: use the id
  res$label[top] <- lab[top]

  if (is.null(novelty)) {
    res$colour_by <- res$change
    colours <- setNames(c("#0072B2", "#D55E00", "grey70"), c(group_names, "not changed"))
  } else {
    stopifnot(length(novelty) == nrow(res))
    novelty <- as.character(novelty)
    cnt <- table(factor(novelty[changed], names(volcano_novelty_colours)))
    lev <- paste0(names(cnt), " (", cnt, ")")                 # legend: how many changed points of each class
    res$colour_by <- ifelse(res$change == "not changed", "not changed", lev[match(novelty, names(cnt))])
    colours <- c(setNames(volcano_novelty_colours, lev), "not changed" = "grey80")
  }
  res <- res[order(res$colour_by != "not changed"), ]      # draw the changed points on top
  ggplot(res, aes(logFC, neglog10fdr, colour = colour_by, label = label)) +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", colour = "grey50") +
    geom_point(size = 2) +
    geom_text_repel(max.overlaps = Inf, size = 3.2, show.legend = FALSE, colour = "black") +
    scale_colour_manual(values = colours, name = NULL) +
    labs(x = "log2 fold change", y = "-log10 FDR", title = title) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank())
}
