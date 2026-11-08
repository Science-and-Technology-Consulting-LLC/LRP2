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
  fill <- function(x) file.path(cfg$run_dir, gsub("{comparison}", cfg$comparison, x, fixed = TRUE))
  c(list(comparison = cfg$comparison, dataset_name = cfg$dataset_name), lapply(cfg$folders, fill))
}
