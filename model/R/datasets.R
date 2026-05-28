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
  english_pilot50 = list(
    bundle         = "long_subset_data_50.rds",
    prepare_script = "N/A (built inline)",
    description    = "50-child English pilot bundle for prior/model validation"
  ),
  english_I200 = list(
    bundle         = "long_subset_data_I200.rds",
    prepare_script = "model/scripts/prepare_longitudinal_data.R \"English (American)\" 200 200",
    description    = "English I=200 J=198 bundle (May-10-era scale) for the 4-panel architecture demo"
  ),
  norwegian = list(
    bundle         = "long_subset_data_nor.rds",
    prepare_script = "model/scripts/prepare_longitudinal_norwegian.R",
    description    = "Norwegian CDI:WS longitudinal subset"
  ),
  babyview = list(
    bundle         = "babyview_subset_data.rds",
    prepare_script = "model/scripts/prepare_babyview.R",
    description    = "BabyView English: longitudinal CDI + observed video input rates"
  ),
  seedlings = list(
    bundle         = "seedlings_subset_data.rds",
    prepare_script = "model/scripts/prepare_seedlings.R",
    description    = "SEEDLingS (Bergelson): WG CDI at 12+18mo + LENA AWC per recording"
  ),
  stanford_linked = list(
    bundle         = "stanford_linked_subset_data.rds",
    prepare_script = "model/scripts/prepare_stanford_linked.R",
    description    = "Stanford TotLot 3 (Adams 2018) item-level CDI joined with Peekbank LWL processing"
  ),
  am2018 = list(
    bundle         = "io_am2018_subset_data.rds",
    prepare_script = "model/scripts/prepare_io_dataset.R am2018",
    description    = "AM2018 (TL3): WG+WS CDI (13-27mo) + LENA AWC @16,18mo; delta_j anchored"
  ),
  fmw2013 = list(
    bundle         = "io_fmw2013_subset_data.rds",
    prepare_script = "model/scripts/prepare_io_dataset.R fmw2013",
    description    = "FMW2013 (TLO): WG+WS CDI (18-30mo) + LENA AWC @18mo; delta_j anchored"
  )
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
