#' classification.R
#'
#' How many transcripts and proteins fall in each SQANTI structural category.
#' Reads the *.best_ORF_SQANTI_classification.tsv table written by S3 (M3_SQANTI_PROTEIN).
#'   FSM  full-splice_match         same splice junctions as a reference transcript or protein
#'   ISM  incomplete-splice_match   a truncated version of a reference
#'   NIC  novel_in_catalog          known splice sites in a new combination
#'   NNC  novel_not_in_catalog      at least one new splice site

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

class_levels <- c("FSM", "ISM", "NIC", "NNC")
class_colours <- c(FSM = "#0072B2", ISM = "#009E73", NIC = "#E69F00", NNC = "#D55E00")

read_classification <- function(path) fread(path, sep = "\t", data.table = FALSE)

# SQANTI writes the long names; use the standard short ones.
short_class <- function(x) {
  factor(c("full-splice_match" = "FSM", "incomplete-splice_match" = "ISM",
           "novel_in_catalog" = "NIC", "novel_not_in_catalog" = "NNC")[as.character(x)],
         levels = class_levels)
}

# Counts per category for the transcript (tx_cat) and the protein (pr_splice_cat) of every isoform.
class_counts <- function(cls) {
  rbind(
    data.frame(level = "Transcript", table(category = short_class(cls$tx_cat))),
    data.frame(level = "Protein",    table(category = short_class(cls$pr_splice_cat)))
  )
}

# Bar chart: one panel for transcripts, one for proteins; the bar label is the count.
class_barplot <- function(counts, title = NULL) {
  counts$level <- factor(counts$level, levels = c("Transcript", "Protein"))
  ggplot(counts, aes(category, Freq, fill = category)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = Freq), vjust = -0.4, size = 4.2) +
    facet_wrap(~ level) +
    scale_fill_manual(values = class_colours, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "SQANTI structural category", y = "number of isoforms", title = title) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA))
}

# Transcript category against protein category: how many novel transcripts still make a known protein, and so on.
class_crosstab <- function(cls) {
  table(transcript = short_class(cls$tx_cat), protein = short_class(cls$pr_splice_cat))
}
