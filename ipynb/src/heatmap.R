#' heatmap.R
#'
#' Clustered heatmap of an expression matrix written by S4 (normalized CPM),
#' with the samples annotated by group. Used for transcript isoforms and ORFs.

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
})

# Samples and their group from the leafcutter groups file (two columns, no header).
read_groups <- function(path) {
  g <- fread(path, sep = "\t", header = FALSE, col.names = c("sample", "group"), data.table = FALSE)
  g
}

# Expression matrix: one row per isoform or ORF, one column per sample in the groups file.
#   id_col  column that names the rows; label_col (optional) a nicer row name
read_expression <- function(path, groups, id_col, label_col = NULL) {
  tab <- fread(path, sep = "\t", data.table = FALSE)
  stopifnot(id_col %in% names(tab), all(groups$sample %in% names(tab)))
  mat <- as.matrix(tab[, groups$sample, drop = FALSE])
  rn <- tab[[id_col]]
  if (!is.null(label_col)) {   # novel isoforms have no name: keep their id
    nm <- tab[[label_col]]
    rn <- ifelse(is.na(nm) | nm == "", rn, nm)
  }
  rownames(mat) <- make.unique(as.character(rn))
  mat
}

# log2(CPM + 1), then each row centred and scaled, so the colour shows relative expression.
scale_rows <- function(mat) {
  z <- t(scale(t(log2(mat + 1))))
  z[is.nan(z)] <- 0
  z
}

# The heatmap. Rows with no expression in any sample are dropped.
expression_heatmap <- function(mat, groups, title = NULL, row_font = 8) {
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  z <- scale_rows(mat)
  levels_g <- unique(groups$group)
  group_colours <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(levels_g)], levels_g)
  top <- HeatmapAnnotation(group = groups$group, col = list(group = group_colours),
                           annotation_name_side = "left")
  Heatmap(z, name = "zscore", top_annotation = top,
          heatmap_legend_param = list(title = "row z-score of log2(CPM+1)"),
          col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
          cluster_rows = nrow(z) > 1, cluster_columns = TRUE,
          column_split = if (length(levels_g) > 1) factor(groups$group, levels = levels_g) else NULL,
          row_names_gp = grid::gpar(fontsize = row_font),
          column_title = title)
}

# Draw the heatmap to figures/<name>.png and figures/<name>.pdf, then show the PNG in the notebook.
# (Drawing straight onto the notebook's own plot device fails for ComplexHeatmap, so it goes to files first.)
show_heatmap <- function(ht, name, width = 9, height = 7, dir = "figures") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  stem <- file.path(dir, name)
  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 150)
  draw(ht); invisible(dev.off())
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  draw(ht); invisible(dev.off())
  if (requireNamespace("IRdisplay", quietly = TRUE)) IRdisplay::display_png(file = paste0(stem, ".png"))
  invisible(stem)
}
