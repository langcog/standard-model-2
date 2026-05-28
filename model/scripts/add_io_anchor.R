## Add delta_j anchor fields (delta_j_prior_mean / _sd) to an existing IO
## bundle, so it can be fit with log_irt_io_anchored.stan. Anchors to the
## big EN longitudinal fit's posterior-median delta_j (join on item
## string); unanchored items get a weak prior at the EN population mean.
##
## Usage:  Rscript model/scripts/add_io_anchor.R <bundle_basename>
##   e.g.  Rscript model/scripts/add_io_anchor.R seedlings_subset_data

source("model/R/config.R")
suppressPackageStartupMessages({library(dplyr); library(readr)})

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args) >= 1) args[1] else stop("need bundle basename")
path <- file.path(PATHS$fits_dir, paste0(base, ".rds"))
b <- readRDS(path)

en_bundle <- readRDS(file.path(PATHS$fits_dir, "long_subset_data.rds"))
en_psi    <- read_csv(file.path(PATHS$fits_dir, "summaries",
                                 "long_no_freq_slopes_psi.csv"),
                      show_col_types = FALSE)
en_delta <- en_bundle$word_info |>
  select(jj, item) |>
  inner_join(en_psi |> select(jj, delta_j_median), by = "jj") |>
  select(item, delta_anchor = delta_j_median)
en_mu <- median(en_delta$delta_anchor)

wi <- b$word_info |>
  left_join(en_delta, by = "item") |>
  mutate(anchored = !is.na(delta_anchor),
         delta_j_prior_mean = ifelse(anchored, delta_anchor, en_mu),
         delta_j_prior_sd   = ifelse(anchored, 0.10, 5.0))

stopifnot(nrow(wi) == b$stan_data$J)
b$word_info <- wi
b$stan_data$delta_j_prior_mean <- wi$delta_j_prior_mean
b$stan_data$delta_j_prior_sd   <- wi$delta_j_prior_sd

saveRDS(b, path)
cat(sprintf("%s: anchored %d / %d items (%d weak-prior)\n",
            base, sum(wi$anchored), nrow(wi), sum(!wi$anchored)))
