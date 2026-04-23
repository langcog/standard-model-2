## Registry of named longitudinal datasets.
##
## Each entry specifies:
##   bundle        - filename (relative to PATHS$fits_dir) of the Stan-ready
##                   bundle produced by the corresponding prepare_*.R
##                   script. The bundle list must have $stan_data,
##                   $word_info, $class_levels, $admin_info, $df, $language.
##   prepare_script - path to the R script that builds the bundle (so
##                   `make` can depend on it).
##   description   - short human-readable description.
##
## Add a new longitudinal dataset by:
##   1) writing a prepare_*.R script that outputs a bundle in the
##      standard format to PATHS$fits_dir.
##   2) adding an entry here.
## Fit scripts and analyze scripts accept a dataset key as argument.

DATASETS <- list(
  english = list(
    bundle         = "long_subset_data.rds",
    prepare_script = "model/scripts/prepare_longitudinal_data.R",
    description    = "English (American) CDI:WS longitudinal subset"
  ),
  norwegian = list(
    bundle         = "long_subset_data_nor.rds",
    prepare_script = "model/scripts/prepare_longitudinal_norwegian.R",
    description    = "Norwegian CDI:WS longitudinal subset"
  )
  # Add e.g. 'peekbank_en' here when we build that pipeline.
)

## Resolve a dataset key to (key, bundle path, bundle object).
## If the dataset key is NULL or empty, defaults to "english".
get_dataset <- function(key = NULL) {
  if (is.null(key) || !nzchar(key)) key <- "english"
  if (!key %in% names(DATASETS)) {
    stop(sprintf("Unknown dataset '%s'. Known: %s",
                 key, paste(names(DATASETS), collapse = ", ")))
  }
  entry <- DATASETS[[key]]
  path <- file.path(PATHS$fits_dir, entry$bundle)
  list(key = key, path = path, entry = entry)
}

## Load a dataset bundle by key (or path if explicitly passed).
load_dataset_bundle <- function(key = NULL, path = NULL) {
  if (!is.null(path) && nzchar(path)) {
    readRDS(path)
  } else {
    info <- get_dataset(key)
    if (!file.exists(info$path))
      stop(sprintf("Dataset bundle '%s' not found at %s. Run %s to build it.",
                   info$key, info$path, info$entry$prepare_script))
    readRDS(info$path)
  }
}
