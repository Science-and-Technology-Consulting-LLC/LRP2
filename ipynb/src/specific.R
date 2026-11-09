#' specific.R
#'
#' Splice variants that belong to one group and not the other. For every transcript isoform, ORF and
#' leafcutter splicing variant we count in how many samples of each group it is detected (reads at or
#' above a minimum). A feature is "specific" to a group when it is detected in at least k samples of
#' that group and in none of the other. This is a presence/absence view; the tests of different usage
#' (DRIMSeq, leafcutter) are in usage.ipynb and splicing.ipynb.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

# Colour of each group, the same as on the other heatmaps.
group_palette <- function(levels_g) setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(levels_g)], levels_g)

# Counts and a description of every feature, from the table the pipeline wrote.
#   kind "isoforms": *_DTE_edgeR_raw_CPM_matrix.txt   kind "orfs": *_DE_ORF_edgeR_raw_CPM_matrix.txt
#   kind "variants": *.lr_leafcutter.subisoform_clusters.txt (intron counts)
# Returns list(counts = matrix, info = data.frame). The columns of counts follow the groups table.
read_features <- function(kind = c("isoforms", "orfs", "variants"), path, groups) {
  kind <- match.arg(kind)
  tab <- fread(path, sep = "\t", data.table = FALSE)
  cols <- if (kind == "variants") paste0(groups$sample, "_counts") else groups$sample
  stopifnot(all(cols %in% names(tab)))
  counts <- as.matrix(tab[, cols, drop = FALSE])
  colnames(counts) <- groups$sample
  if (kind == "isoforms") {
    id <- tab$isoform_id
    gene <- ifelse(is.na(tab$gene_name) | tab$gene_name == "", sub("::.*$", "", id), tab$gene_name)
    label <- ifelse(is.na(tab$transcript_name) | tab$transcript_name == "", id, tab$transcript_name)
    info <- data.frame(id = id, gene = gene, label = label)
  } else if (kind == "orfs") {
    id <- tab$orf_base_id
    info <- data.frame(id = id, gene = sub("::.*$", "", id), label = id)
  } else {
    gene <- sub("::.*$", "", sub(",.*$", "", tab$isoforms))
    starts <- format(as.numeric(sub(",.*$", "", tab$intron_starts)), scientific = FALSE, trim = TRUE)
    ends <- format(as.numeric(sub("^.*,", "", tab$intron_ends)), scientific = FALSE, trim = TRUE)
    n_introns <- lengths(strsplit(tab$intron_starts, ","))
    label <- paste0(gene, " ", tab$chr, ":clu_", tab$cluster_idx, " ", starts, "-", ends, ifelse(n_introns > 1, paste0(" (", n_introns, " introns)"), ""))
    info <- data.frame(id = tab$subisoform_id, gene = gene, label = label, chr = tab$chr, strand = tab$strand,
                       intron_starts = tab$intron_starts, intron_ends = tab$intron_ends)
  }
  rownames(counts) <- make.unique(as.character(info$label))
  list(counts = counts, info = info)
}

# For each feature: in how many samples of each group it is detected (n_<group>), its mean count in each
# group (mean_<group>), and `specific_to`: the group it is specific to, or "not specific".
#   min_count  reads needed in a sample to call the feature detected
#   k          samples of a group needed (never more than the group has)
# Needs exactly two groups in `groups`.
specificity <- function(feat, groups, min_count = 1, k = 2) {
  lev <- unique(groups$group)
  stopifnot(length(lev) == 2)
  cols <- setNames(lapply(lev, function(g) groups$sample[groups$group == g]), lev)
  kk <- min(k, lengths(cols))
  det <- feat$counts >= min_count
  n <- do.call(cbind, lapply(cols, function(s) rowSums(det[, s, drop = FALSE])))
  mean_count <- do.call(cbind, lapply(cols, function(s) rowMeans(feat$counts[, s, drop = FALSE])))
  spec <- ifelse(n[, 1] >= kk & n[, 2] == 0, lev[1], ifelse(n[, 2] >= kk & n[, 1] == 0, lev[2], "not specific"))
  pal <- col2rgb(group_palette(lev))
  rgb_txt <- setNames(apply(pal, 2, paste, collapse = ","), lev)
  out <- cbind(feat$info,
               setNames(as.data.frame(n), paste0("n_", lev)),
               setNames(as.data.frame(mean_count), paste0("mean_", lev)))
  out$specific_to <- factor(spec, levels = c(lev, "not specific"))
  out$colour <- ifelse(spec %in% lev, rgb_txt[spec], "")
  rownames(out) <- rownames(feat$counts)
  attr(out, "k") <- kk
  out
}

# How many features are specific to each group, shared or not specific, at each level.
specific_summary <- function(tabs) {
  do.call(rbind, lapply(names(tabs), function(l) {
    t <- as.data.frame(table(specific_to = tabs[[l]]$specific_to), responseName = "features")
    cbind(level = l, t)
  }))
}

# Bar chart of that summary: one panel per level, one bar per group, grey for not specific.
specific_barplot <- function(summary, title = NULL) {
  lev <- setdiff(levels(summary$specific_to), "not specific")
  cols <- c(group_palette(lev), "not specific" = "grey70")
  ggplot(summary, aes(specific_to, features, fill = specific_to)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = features), vjust = -0.3, size = 3.5) +
    facet_wrap(~ level, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "features", title = title) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
}

# The specific features of the chosen genes, all levels together.
specific_genes <- function(tabs, genes) {
  rows <- lapply(names(tabs), function(l) {
    t <- tabs[[l]]
    t <- t[t$gene %in% genes & t$specific_to != "not specific", , drop = FALSE]
    if (nrow(t) == 0) return(NULL)
    cbind(level = l, t[, c("gene", "label", "specific_to", grep("^n_", names(t), value = TRUE), grep("^mean_", names(t), value = TRUE))])
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) message("none of ", paste(genes, collapse = ", "), " has a specific feature")
  else rownames(out) <- NULL
  out
}

# Genes with the most specific features, with how many are specific to each group.
specific_by_gene <- function(tab) {
  t <- tab[tab$specific_to != "not specific", , drop = FALSE]
  if (nrow(t) == 0) return(data.frame())
  by_gene <- as.data.frame.matrix(table(t$gene, droplevels(t$specific_to)))
  by_gene <- cbind(gene = rownames(by_gene), by_gene, total = rowSums(by_gene))
  by_gene[order(-by_gene$total, by_gene$gene), ]
}

# Heatmap of the specific features only: rows grouped by the group they belong to, colour = log2(count + 1).
#   novelty (optional): one entry per row of `tab` ("known", "novel", "no ORF"), drawn as a bar at the right
specific_heatmap <- function(feat, tab, groups, title = NULL, row_font = NULL, novelty = NULL) {
  keep <- which(tab$specific_to != "not specific")
  if (length(keep) == 0) {
    message("nothing to draw: no specific features (", if (is.null(title)) "no title" else title, ")")
    return(invisible(NULL))
  }
  ord <- keep[order(tab$specific_to[keep], tab$gene[keep], tab$label[keep])]
  m <- log2(feat$counts[ord, groups$sample, drop = FALSE] + 1)
  lay <- row_layout(nrow(m), row_font)   # row_layout() is in heatmap.R
  levels_g <- unique(groups$group)
  top <- HeatmapAnnotation(group = groups$group, col = list(group = group_palette(levels_g)), annotation_name_side = "left")
  side <- if (is.null(novelty)) NULL else
    rowAnnotation(novelty = factor(novelty[ord], c("known", "novel", "no ORF")), col = list(novelty = novelty_colours))
  Heatmap(m, name = "log2count", top_annotation = top, right_annotation = side,
          heatmap_legend_param = list(title = "log2(count + 1)"),
          col = colorRamp2(c(0, max(m, 1)), c("white", "#35264C")),
          cluster_rows = FALSE, cluster_columns = TRUE,
          row_split = droplevels(tab$specific_to[ord]),
          column_split = factor(groups$group, levels = levels_g),
          show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = lay$font),
          height = unit(lay$body_in, "in"), column_title = title)
}

# Counts of the same features restricted to some of the samples (used to leave one out).
use_samples <- function(feat, groups) {
  feat$counts <- feat$counts[, groups$sample, drop = FALSE]
  feat
}
