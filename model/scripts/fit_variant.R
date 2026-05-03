## Fit a single named variant on the cached Wordbank subsample.
##
## Usage:   Rscript model/scripts/fit_variant.R <variant_name>
## Variants: baseline | fix_delta | fix_s | both_fixed | 2pl | 2pl_fix_delta
## Reads:   fits/subset_data.rds
## Writes:  fits/wordbank_<variant_name>.rds

source("model/R/config.R")
source("model/R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Pass a variant name.")
nm <- args[1]

data_file <- file.path(PATHS$fits_dir, "subset_data.rds")
if (!file.exists(data_file))
  stop("Run `make data` first to produce subset_data.rds.")
bundle <- readRDS(data_file)

overrides <- variant_hyperpriors(nm)
stan_data <- modifyList(bundle$stan_data, overrides)

cat(sprintf("\n===== Fitting variant: %s =====\n", nm))
cat("Hyperprior overrides:\n")
if (length(overrides) == 0) cat("  (none)\n")
for (k in names(overrides)) cat(sprintf("  %-25s = %s\n", k, overrides[[k]]))

fit <- fit_variant(stan_data, tag = sprintf("wordbank_%s", nm))

cat("\n--- Scalar posteriors ---\n")
pars <- c("sigma_alpha", "s", "delta", "pi_alpha", "sigma_xi")
if (grepl("2pl", nm))     pars <- c(pars, "sigma_lambda")
if (grepl("slopes", nm))  pars <- c(pars, "sigma_zeta")
print(summarize_fit(fit, pars = pars), digits = 3)

cat(sprintf("\nSaved: fits/wordbank_%s.rds\n", nm))
