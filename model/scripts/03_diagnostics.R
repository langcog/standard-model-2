## 2x2 diagnostic design: vary whether delta and s are free or fixed.
##
## variant         | delta        | s
## ----------------|--------------|-------------
## baseline        | free         | free
## fix_delta       | fixed at 0   | free
## fix_s           | free         | fixed at 2
## both_fixed      | fixed at 0   | fixed at 2
##
## Usage:   Rscript model/scripts/03_diagnostics.R
## Reads:   model/fits/subset_data.rds (from 02_prepare_data.R)
## Outputs: model/fits/wordbank_{variant}.rds for each variant,
##          model/fits/loo_compare.rds

source("model/R/config.R")
source("model/R/helpers.R")

data_file <- file.path(PATHS$fits_dir, "subset_data.rds")
if (!file.exists(data_file))
  stop("Run 02_prepare_data.R first to produce subset_data.rds.")
bundle <- readRDS(data_file)
base_data <- bundle$stan_data

variants <- c("baseline", "fix_delta", "fix_s", "both_fixed")

fits <- list()
for (nm in variants) {
  overrides <- variant_hyperpriors(nm)
  stan_data <- modifyList(base_data, overrides)
  fits[[nm]] <- fit_variant(stan_data, tag = sprintf("wordbank_%s", nm))
  cat(sprintf("\n--- %s: scalar posteriors ---\n", nm))
  print(summarize_fit(fits[[nm]]), digits = 3)
}

## Scalar-posterior summary table across variants (LOO-CV deferred until
## the Stan model's log_lik is re-enabled).
cat("\n===== Scalar posteriors across variants =====\n")
summary_tbl <- do.call(rbind, lapply(names(fits), function(nm) {
  s <- summarize_fit(fits[[nm]])
  s$variant <- nm
  s
}))
summary_tbl <- summary_tbl[, c("variant", "param", "median", "lo95", "hi95",
                               "n_eff", "Rhat")]
print(summary_tbl, digits = 3)

saveRDS(summary_tbl, file.path(PATHS$fits_dir, "summary_tbl.rds"))
cat("\nDone. Saved model/fits/summary_tbl.rds\n")
