#' paths.R
#'
#' Where a run's results are, set once.
#' `paths.json` holds the run folder, the comparison name and the dataset name; the folders
#' inside the run are written there as short relative names. Every notebook calls
#' `run_locations()` and uses names such as `loc$leafcutter` or `loc$de_gene`.
#' To point the notebooks at another run, edit `paths.json` and nothing else.
#' (The environment variable LRP2_PATHS can name a different config file.)

# Returns a list: comparison, dataset_name, and the full path of every folder in `folders`.
run_locations <- function(config = Sys.getenv("LRP2_PATHS", "paths.json")) {
  cfg <- jsonlite::fromJSON(config)
  options(lrp2_figure_prefix = cfg$comparison)   # saved figures start with the comparison name (see figure_stem)
  fill <- function(x) file.path(cfg$run_dir, gsub("{comparison}", cfg$comparison, x, fixed = TRUE))
  c(list(run_dir = cfg$run_dir, comparison = cfg$comparison, dataset_name = cfg$dataset_name), lapply(cfg$folders, fill))
}

# Where a figure is saved: figures/<comparison>_<name> (no extension). The comparison name is set by run_locations().
figure_stem <- function(name, dir = "figures") {
  prefix <- getOption("lrp2_figure_prefix", "")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  file.path(dir, paste0(if (nzchar(prefix)) paste0(prefix, "_"), name))
}

# Save a ggplot as PNG and PDF under its comparison-named file and return the plot, so the notebook still shows it.
# Size defaults to the notebook's plot size (repr.plot.width and repr.plot.height).
save_plot <- function(p, name, width = getOption("repr.plot.width", 7), height = getOption("repr.plot.height", 5),
                      dir = "figures") {
  stem <- figure_stem(name, dir)
  ggplot2::ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 150)
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width, height = height)
  p
}
