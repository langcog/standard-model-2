## Prepare the stratified Wordbank subsample + input-rate prior and cache
## the Stan data list for downstream fits.
##
## Usage:   Rscript model/scripts/02_prepare_data.R [n_children] [n_items]
## Outputs: model/fits/subset_data.rds (the stan_data + metadata used by
##          every Wordbank-subset-based fit)
##
## Note: all fits that load this file see identical data, making LOO-CV
## comparisons across diagnostic variants valid.

source("model/R/config.R")
source("model/R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
n_children <- as.integer(if (length(args) >= 1) args[1] else 500)
n_items    <- as.integer(if (length(args) >= 2) args[2] else 200)
SEED <- 20250420

message(sprintf("Preparing subsample: %d children, %d items, seed=%d",
                n_children, n_items, SEED))

df <- load_wordbank_data()
prior_r <- load_input_rate_prior()
message(sprintf("Input prior: log_r ~ N(%.3f, %.3f^2) from n=%d samples",
                prior_r$mu_r, prior_r$sigma_r, prior_r$n))

sub <- subsample_wordbank(df, n_children = n_children,
                         n_items = n_items, seed = SEED)
bundle <- build_stan_data(sub, prior_r)

message(sprintf("Subsample: %d children, %d items, %d obs (mean produces=%.3f)",
                bundle$stan_data$I, bundle$stan_data$J,
                bundle$stan_data$N, mean(bundle$stan_data$y)))

saveRDS(bundle, file.path(PATHS$fits_dir, "subset_data.rds"))
cat("\nSaved to model/fits/subset_data.rds\n")
