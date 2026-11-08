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

# Remove named samples (for example an outlier) from the groups table. Names that are not found are reported.
drop_samples <- function(groups, samples) {
  missing <- setdiff(samples, groups$sample)
  if (length(missing) > 0) warning("not in the groups table: ", paste(missing, collapse = ", "))
  groups[!groups$sample %in% samples, , drop = FALSE]
}

# Number of isoforms or ORFs with any expression in each sample, lowest first: a very low sample is a coverage outlier.
detected_per_sample <- function(mat) sort(colSums(mat > 0))

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

# Size of the rows. Each row gets 1.2 times its font size in points, so the labels do not overlap and
# every row can be read. A PDF page is at most 200 inches, so with too many rows for 190 inches the
# font is made smaller to fit. The default font is 8 points up to 60 rows and 4 points above that.
row_layout <- function(n, row_font = NULL) {
  if (is.null(row_font)) row_font <- if (n <= 60) 8 else 4
  row_pt <- row_font * 1.2
  max_pt <- 190 * 72
  if (n * row_pt > max_pt) {
    row_pt <- max_pt / n
    row_font <- row_pt / 1.2
    message(n, " rows: the font is shrunk to ", round(row_font, 1), " points to fit one page; zoom into the PDF to read it")
  }
  list(font = row_font, row_pt = row_pt, body_in = max(n * row_pt, 250) / 72)
}

# The heatmap. Rows with no expression in any sample are dropped. Every row has its label, and the
# heatmap is as tall as the rows need (see row_layout); `row_font` in points overrides the default.
expression_heatmap <- function(mat, groups, title = NULL, row_font = NULL) {
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  z <- scale_rows(mat)
  lay <- row_layout(nrow(z), row_font)
  levels_g <- unique(groups$group)
  group_colours <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(levels_g)], levels_g)
  top <- HeatmapAnnotation(group = groups$group, col = list(group = group_colours),
                           annotation_name_side = "left")
  Heatmap(z, name = "zscore", top_annotation = top,
          heatmap_legend_param = list(title = "row z-score of log2(CPM+1)"),
          col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
          cluster_rows = nrow(z) > 1, cluster_columns = TRUE,
          column_split = if (length(levels_g) > 1) factor(groups$group, levels = levels_g) else NULL,
          show_row_names = TRUE,
          row_names_gp = grid::gpar(fontsize = lay$font),
          height = unit(lay$body_in, "in"),
          column_title = title)
}

# Draw the heatmap to figures/<name>.png and figures/<name>.pdf, then show the PNG in the notebook.
# (Drawing straight onto the notebook's own plot device fails for ComplexHeatmap, so it goes to files first.)
show_heatmap <- function(ht, name, width = 10, height = NULL, dir = "figures") {
  # a tall heatmap sets its own height: the rows (see row_layout) plus 3.5 inches for title, group bar and sample names
  if (is.null(height)) {
    lay <- row_layout(nrow(ht@matrix), ht@row_names_param$gp$fontsize)
    height <- lay$body_in + 3.5
  }
  res <- max(50, min(150, 30000 / height))   # 150 dpi, or less for a very tall figure to stay under 30000 pixels
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  stem <- file.path(dir, name)
  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = res)
  draw(ht); invisible(dev.off())
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  draw(ht); invisible(dev.off())
  if (requireNamespace("IRdisplay", quietly = TRUE)) IRdisplay::display_png(file = paste0(stem, ".png"))
  invisible(stem)
}
