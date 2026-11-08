#' splicing.R
#'
#' Matrix of alternative splicing events from the long-read leafcutter output of S4.
#' A leafcutter cluster is a set of introns that share splice sites; each intron (a splicing
#' variant) is one row. For every sample we show PSI: the share of the cluster's reads that use
#' that intron. Samples that use exons differently show up as columns that cluster apart.

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
})

# Intron counts per sample for every splicing variant (the *.lr_leafcutter.subisoform_clusters.txt table).
read_subisoforms <- function(path) {
  fread(path, sep = "\t", data.table = FALSE)
}

# Cluster significance (*.ds_cluster_significance.txt) and effect sizes (*.ds_effect_sizes.txt).
read_cluster_significance <- function(path) fread(path, sep = "\t", na.strings = "NA", data.table = FALSE)
read_effect_sizes <- function(path) fread(path, sep = "\t", na.strings = "NA", data.table = FALSE)

# PSI matrix: one row per splicing variant, one column per sample (in the order of the groups table).
#   Returns a list: psi (matrix), info (one row per variant: cluster, gene, deltaPSI, adjusted p).
subisoform_psi <- function(subs, groups, sig, eff) {
  counts <- as.matrix(subs[, paste0(groups$sample, "_counts"), drop = FALSE])
  colnames(counts) <- groups$sample
  cluster <- paste0(subs$chr, ":clu_", subs$cluster_idx)
  total <- rowsum(counts, cluster)[cluster, , drop = FALSE]
  psi <- counts / total
  psi[total == 0] <- NA                       # the cluster has no reads in that sample
  gene <- sub("::.*$", "", sub(",.*$", "", subs$isoforms))
  starts <- sub(",.*$", "", subs$intron_starts); ends <- sub("^.*,", "", subs$intron_ends)
  n_introns <- lengths(strsplit(subs$intron_starts, ","))
  label <- paste0(gene, " ", cluster, " ", format(as.numeric(starts), scientific = FALSE, trim = TRUE),
                  "-", format(as.numeric(ends), scientific = FALSE, trim = TRUE),
                  ifelse(n_introns > 1, paste0(" (", n_introns, " introns)"), ""))
  rownames(psi) <- make.unique(label)
  delta <- eff[match(subs$subisoform_id, eff$intron), grep("^deltapsi_", names(eff))[1]]
  # every gene named by the isoforms that use the variant (a cluster can span more than one gene)
  genes <- vapply(strsplit(subs$isoforms, ","), function(v) paste(unique(sub("::.*$", "", v)), collapse = ";"), "")
  info <- data.frame(cluster = cluster, gene = gene, genes = genes, delta_psi = delta,
                     adj_p = sig$p.adjust[match(cluster, sig$cluster)], row.names = rownames(psi))
  list(psi = psi, info = info)
}

# The PSI matrix for the samples left after removing those named in `exclude` (none removed: all samples).
# Returns the matrix and the groups table that goes with it. The leafcutter test results (cluster p values and
# deltaPSI) were computed by the pipeline with all samples; only the PSI shown here changes.
psi_set <- function(subs, groups, sig, eff, exclude = character(0)) {
  g <- if (length(exclude) > 0) drop_samples(groups, exclude) else groups   # drop_samples() is in heatmap.R
  list(x = subisoform_psi(subs, g, sig, eff), groups = g)
}

# Keep the variants of clusters that leafcutter found different between the groups (adjusted p below `alpha`),
# and only clusters that have at least two variants. Clusters leafcutter did not test (little coverage, too
# few samples with reads) have no p value and are left out, even with alpha = 1.
select_variants <- function(x, alpha = 0.05) {
  multi <- ave(seq_len(nrow(x$info)), x$info$cluster, FUN = length) > 1
  keep <- multi & !is.na(x$info$adj_p) & x$info$adj_p < alpha
  list(psi = x$psi[keep, , drop = FALSE], info = x$info[keep, , drop = FALSE])
}

# Every splicing variant of the named genes, whether or not the cluster passed the test (or was tested at all),
# in clusters that have at least two variants. Genes with nothing to show are reported.
select_genes <- function(x, genes) {
  in_gene <- function(g) vapply(strsplit(x$info$genes, ";"), function(v) g %in% v, TRUE)
  hits <- vapply(genes, function(g) any(in_gene(g)), TRUE)
  multi <- ave(seq_len(nrow(x$info)), x$info$cluster, FUN = length) > 1
  keep <- multi & Reduce(`|`, lapply(genes, in_gene))
  if (any(!hits)) message("not in any isoform of this run: ", paste(genes[!hits], collapse = ", "))
  once <- genes[hits & !vapply(genes, function(g) any(in_gene(g) & multi), TRUE)]
  if (length(once) > 0) message("only in clusters with a single variant (nothing to compare): ", paste(once, collapse = ", "))
  out <- list(psi = x$psi[keep, , drop = FALSE], info = x$info[keep, , drop = FALSE])
  o <- order(out$info$gene, out$info$cluster, rownames(out$psi))
  list(psi = out$psi[o, , drop = FALSE], info = out$info[o, , drop = FALSE])
}

# The heatmap: variants in rows, samples in columns split by group, colour = PSI.
#   Row side shows the change in PSI between the groups (deltaPSI, second group minus first).
splicing_heatmap <- function(x, groups, title = NULL, row_font = NULL, split_by_gene = FALSE) {
  lay <- row_layout(nrow(x$psi), row_font)   # row_layout() is in heatmap.R, which the notebook loads first
  levels_g <- unique(groups$group)
  group_colours <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(levels_g)], levels_g)
  top <- HeatmapAnnotation(group = groups$group, col = list(group = group_colours), annotation_name_side = "left")
  lim <- max(0.1, max(abs(x$info$delta_psi), na.rm = TRUE))
  side <- rowAnnotation(deltaPSI = x$info$delta_psi,
                        col = list(deltaPSI = colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B"))))
  Heatmap(x$psi, name = "PSI", top_annotation = top, right_annotation = side,
          col = colorRamp2(c(0, 1), c("white", "#35264C")), na_col = "grey85",
          cluster_rows = nrow(x$psi) > 1 && !split_by_gene, cluster_columns = TRUE,
          row_split = if (split_by_gene) factor(x$info$gene, levels = unique(x$info$gene)) else NULL,
          column_split = factor(groups$group, levels = levels_g),
          show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = lay$font),
          height = unit(lay$body_in, "in"), column_title = title)
}
