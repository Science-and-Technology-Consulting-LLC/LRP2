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
    usage$p     <- usage$pvalue_gene
  } else {
    usage$gene  <- sub("::.*$", "", usage$orf_base_id)
    usage$label <- usage$orf_base_id
    usage$adj_p <- usage$adj_pvalue_gene
    usage$p     <- usage$pvalue_gene
  }
  usage
}

# One row per gene with its raw and adjusted p value, most significant first (genes DRIMSeq did not test come last).
gene_usage_table <- function(usage) {
  tab <- unique(usage[, c("gene", "p", "adj_p")])
  tab[order(tab$adj_p, tab$p, na.last = TRUE), ]
}

# How many genes were tested, how many pass, and the smallest p values.
usage_summary <- function(usage) {
  tab <- gene_usage_table(usage)
  tested <- tab[!is.na(tab$adj_p), ]
  data.frame(genes_in_table = nrow(tab), genes_tested = nrow(tested),
             adj_p_below_0.05 = sum(tested$adj_p < 0.05), adj_p_below_0.1 = sum(tested$adj_p < 0.1),
             smallest_p = signif(min(tested$p), 3), smallest_adj_p = signif(min(tested$adj_p), 3))
}

# Histogram of the raw p value of every tested gene. Roughly flat means the groups do not differ in
# isoform usage at this sample size; a pile of genes at the left (p near 0) means real differences.
pvalue_histogram <- function(usage, title = NULL) {
  tab <- gene_usage_table(usage)
  tab <- tab[!is.na(tab$p), ]
  ggplot(tab, aes(p)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "#0072B2", colour = "white") +
    labs(x = "raw p value (one per tested gene)", y = "genes", title = title) +
    theme_bw(base_size = 12) + theme(panel.grid = element_blank())
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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), panel.grid = element_blank(),
          legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 8))
}

# Only genes DRIMSeq found different between the groups (adjusted p below `alpha`), one above the other (ncol = 1
# keeps each plot wide enough to read the sample names and the group labels). If more than n pass, the n with the
# smallest adjusted p are drawn and the rest are reported. Returns NULL, with a message, when no gene passes.
usage_panel <- function(usage, groups, alpha = 0.05, n = 4, ncol = 1) {
  tab <- gene_usage_table(usage)
  sig <- tab$gene[!is.na(tab$adj_p) & tab$adj_p < alpha]
  if (length(sig) == 0) {
    message("no gene has an adjusted p below ", alpha, ": nothing to draw")
    return(invisible(NULL))
  }
  message(length(sig), " genes have an adjusted p below ", alpha, if (length(sig) > n) paste0("; drawing the ", n, " with the smallest"))
  wrap_plots(lapply(head(sig, n), function(g) usage_plot(usage, groups, g)), ncol = ncol)
}
