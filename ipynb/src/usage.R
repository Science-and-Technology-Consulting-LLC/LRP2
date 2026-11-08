#' usage.R
#'
#' Differential usage (DRIMSeq) of transcript isoforms and ORFs: for a gene, the share each
#' isoform takes of that gene's expression in every sample, drawn as stacked bars by group.
#' Reads the *_DTU_transcript_DRIMSeq_summary.txt and *_DU_ORF_DRIMSeq_summary.txt tables from S4.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

# Read a DRIMSeq summary table (tab separated, NA as text).
read_usage <- function(path) {
  fread(path, sep = "\t", na.strings = "NA", data.table = FALSE)
}

# Add the columns the plots need: gene (name or id) and the name of each isoform or ORF.
#   kind "transcript": gene_name and transcript_name; kind "orf": gene taken from the orf_base_id.
prepare_usage <- function(usage, kind = c("transcript", "orf")) {
  kind <- match.arg(kind)
  if (kind == "transcript") {
    usage$gene  <- ifelse(is.na(usage$gene_name), usage$gene_id, usage$gene_name)
    usage$label <- ifelse(is.na(usage$transcript_name), usage$isoform_id, usage$transcript_name)
    usage$adj_p <- usage$adj_pvalue_gene
  } else {
    usage$gene  <- sub("::.*$", "", usage$orf_base_id)
    usage$label <- usage$orf_base_id
    usage$adj_p <- usage$adj_pvalue_gene
  }
  usage
}

# One row per gene with its adjusted p value, most significant first (genes DRIMSeq did not test come last).
gene_usage_table <- function(usage) {
  tab <- unique(usage[, c("gene", "adj_p")])
  tab[order(tab$adj_p, na.last = TRUE), ]
}

# Stacked bars for one gene: each column a sample (ordered by group), each colour an isoform.
usage_plot <- function(usage, groups, gene) {
  d <- usage[usage$gene == gene, ]
  props <- paste0(groups$sample, "_prop")
  long <- do.call(rbind, lapply(seq_along(props), function(i) {
    data.frame(label = d$label, sample = groups$sample[i], group = groups$group[i],
               proportion = d[[props[i]]])
  }))
  long$sample <- factor(long$sample, levels = groups$sample)
  p <- unique(d$adj_p)[1]
  ptxt <- if (is.na(p)) "not tested" else paste0("adj. p = ", signif(p, 2))
  ggplot(long, aes(sample, proportion, fill = label)) +
    geom_col(width = 0.85) +
    facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    scale_fill_viridis_d(name = NULL, option = "turbo", end = 0.95) +
    labs(x = NULL, y = "share of the gene's expression", title = paste0(gene, " (", ptxt, ")")) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid = element_blank(),
          legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 8))
}

# The n most significant genes side by side.
usage_panel <- function(usage, groups, n = 4) {
  genes <- head(gene_usage_table(usage)$gene, n)
  wrap_plots(lapply(genes, function(g) usage_plot(usage, groups, g)), ncol = 2)
}
