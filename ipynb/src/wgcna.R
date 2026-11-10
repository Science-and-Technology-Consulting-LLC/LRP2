#' wgcna.R
#'
#' Weighted gene co-expression network analysis (WGCNA) of transcript isoforms: which isoforms rise and
#' fall together across the samples, and which of those sets of isoforms follow the group.
#' The steps and settings follow the endotypes-transcriptomics notebooks 12 to 16
#' (https://github.com/NIH-NLM/endotypes-transcriptomics), which use the WGCNA authors' recommendations.
#' Every function here does one step, and the notebook calls them in order.
#'
#' References
#'   Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis.
#'     BMC Bioinformatics 2008;9:559. doi:10.1186/1471-2105-9-559
#'   WGCNA FAQ (Langfelder and Horvath), https://edo98811.github.io/WGCNA_official_documentation/faq.html
#'   Booeshaghi AS, Hallgrimsdottir IB, Galvez-Merchan A, Pachter L. Depth normalization for single-cell
#'     genomics count data. bioRxiv 2022. doi:10.1101/2022.05.06.490859

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Install WGCNA when it is missing. It needs three Bioconductor packages, which BiocManager installs.
need_wgcna <- function() {
  if (requireNamespace("WGCNA", quietly = TRUE)) return(invisible(TRUE))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(c("GO.db", "impute", "preprocessCore"), ask = FALSE, update = FALSE)
  install.packages("WGCNA", repos = "https://cloud.r-project.org")
  invisible(requireNamespace("WGCNA"))
}

# Remove rows with low counts (WGCNA FAQ rule for RNA-seq counts): drop a row when fewer than `min_count`
# reads are found in more than `max_share` of the samples. Returns the filtered counts.
filter_low_counts <- function(counts, min_count = 10, max_share = 0.90) {
  counts[rowMeans(counts < min_count) <= max_share, , drop = FALSE]
}

# Proportional fitting: scale every sample so that its total equals the mean total over samples.
prop_fit <- function(X) t(t(X) / colSums(X) * mean(colSums(X)))

# PFlog1pPF (Booeshaghi et al. 2022): proportional fitting, log(1 + x), proportional fitting again.
pflog1ppf <- function(counts) prop_fit(log1p(prop_fit(counts)))

# Power for the soft threshold, chosen as in endotypes-transcriptomics notebook 13: the lowest power whose
# scale-free fit (R squared) is at least 0.90 and has flattened; if none does, the WGCNA FAQ value for a
# signed network by number of samples (under 20: 18, 20-30: 16, 30-40: 14, over 40: 12).
# Returns list(power, fit = table of the fit for each power, rule = which rule chose it).
choose_power <- function(datExpr, powers = c(1:10, seq(12, 30, 2))) {
  invisible(capture.output(sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)))
  fit <- data.frame(power = sft$fitIndices$Power,
                    scale_free_R2 = round(-sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq, 3),
                    mean_connectivity = round(sft$fitIndices$mean.k., 1))
  flat <- fit$power[fit$scale_free_R2 >= 0.90 & c(FALSE, abs(diff(fit$scale_free_R2)) < 0.01)]
  n <- nrow(datExpr)
  faq <- if (n < 20) 18 else if (n <= 30) 16 else if (n <= 40) 14 else 12
  list(power = if (length(flat)) min(flat) else faq, fit = fit,
       rule = if (length(flat)) "lowest power with R2 >= 0.90 where the curve has flattened" else paste0("WGCNA FAQ value for a signed network with ", n, " samples"))
}

# Modules of isoforms that behave alike (signed network and signed TOM, as in endotypes-transcriptomics notebook 14).
# `datExpr` has one row per sample and one column per isoform. Returns the WGCNA result (colors = module of each isoform).
find_modules <- function(datExpr, power, min_module_size = 30, merge_cut_height = 0.25, seed = 1) {
  set.seed(seed)
  WGCNA::blockwiseModules(datExpr, power = power, networkType = "signed", TOMType = "signed",
                          minModuleSize = min_module_size, mergeCutHeight = merge_cut_height, maxBlockSize = 5000,
                          numericLabels = FALSE, pamRespectsDendro = FALSE, verbose = 0)
}

# Test every module eigengene (the first principal component of the module, one value per sample) against a 0/1 trait:
# Pearson correlation (the point-biserial correlation, the same as a t test of the eigengene between the groups),
# its p value, and the Benjamini-Hochberg q value over the modules. r above 0 means higher in the trait = 1 group.
module_trait_table <- function(ME, trait) {
  ME <- ME[, colnames(ME) != "MEgrey", drop = FALSE]
  out <- do.call(rbind, lapply(colnames(ME), function(me) {
    ct <- cor.test(ME[, me], trait)
    data.frame(module = sub("^ME", "", me), r = unname(ct$estimate), p = ct$p.value)
  }))
  out$q <- p.adjust(out$p, "BH")
  out[order(out$p), ]
}

# Hub isoforms: for each module the isoforms most correlated with its eigengene (kME), best first.
hub_table <- function(datExpr, ME, colors, labels, genes, n = 5) {
  kME <- cor(datExpr, ME)
  mods <- setdiff(unique(colors), "grey")
  do.call(rbind, lapply(mods, function(mo) {
    g <- names(colors)[colors == mo]
    k <- kME[g, paste0("ME", mo)]
    top <- head(order(k, decreasing = TRUE), n)
    data.frame(module = mo, rank = seq_along(top), isoform = labels[g[top]], gene = genes[g[top]], kME = round(k[top], 3), row.names = NULL)
  }))
}

# Bar chart of the module-group correlations: one bar per module, coloured by the module, starred when q < alpha.
module_trait_plot <- function(tab, alpha = 0.05, title = NULL, case_name = "case", control_name = "control") {
  tab$module <- factor(tab$module, levels = rev(tab$module))
  tab$star <- ifelse(tab$q < alpha, "*", "")
  ggplot(tab, aes(r, module, fill = as.character(module))) +
    geom_col(width = 0.7) +
    geom_text(aes(label = star), hjust = ifelse(tab$r >= 0, -0.3, 1.3), size = 6) +
    scale_fill_identity() +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    geom_vline(xintercept = 0) +
    labs(x = paste0("correlation with group (above 0: higher in ", case_name, ", below 0: higher in ", control_name, ")"),
         y = NULL, title = title, caption = paste0("* q < ", alpha)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank())
}
