## Fit the longitudinal log-linear IRT accumulator model.
##
## Usage:   Rscript model/scripts/fit_longitudinal.R [variant]
## Variants (via sigma_lambda_prior_sd and sigma_zeta_prior_sd):
##   long_2pl         -- 2PL, no per-child slopes
##   long_2pl_slopes  -- 2PL + per-child slopes (default)
##   long_slopes      -- Rasch + per-child slopes
##   long_baseline    -- Rasch, no slopes
##
## Reads:   model/fits/long_subset_data.rds
## Writes:  model/fits/<variant>.rds

source("model/R/config.R")
source("model/R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
nm <- if (length(args) >= 1) args[1] else "long_2pl_slopes"

# A trailing "_nor" in the variant name selects the Norwegian bundle;
# otherwise use the English longitudinal bundle. The core variant
# identity (baseline / 2pl / slopes / 2pl_slopes) is derived by
# stripping the language suffix.
if (grepl("_nor$", nm)) {
  data_file <- file.path(PATHS$fits_dir, "long_subset_data_nor.rds")
  base_nm   <- sub("_nor$", "", nm)
} else {
  data_file <- file.path(PATHS$fits_dir, "long_subset_data.rds")
  base_nm   <- nm
}

bundle <- readRDS(data_file)
base_data <- bundle$stan_data
cat(sprintf("Data bundle: %s (language=%s, I=%d, A=%d, J=%d, N=%d)\n",
            data_file,
            ifelse(is.null(bundle$language), "?", bundle$language),
            base_data$I, base_data$A, base_data$J, base_data$N))

overrides <- switch(base_nm,
  long_baseline    = list(),
  long_2pl         = list(sigma_lambda_prior_sd = 1),
  long_slopes      = list(sigma_zeta_prior_sd = 1),
  long_2pl_slopes  = list(sigma_lambda_prior_sd = 1,
                          sigma_zeta_prior_sd = 1),
  stop(sprintf("Unknown longitudinal variant: %s", base_nm))
)

stan_data <- modifyList(base_data, overrides)

cat(sprintf("\n===== Longitudinal variant: %s =====\n", nm))
cat("Overrides:\n"); str(overrides)

cfg <- modifyList(DEFAULT_FIT_CONFIG,
                  list(chains = 4, iter = 1500, warmup = 750,
                       adapt_delta = 0.95))

fit <- fit_variant(stan_data, tag = nm,
                   cfg = cfg,
                   model_path = file.path(PROJECT_ROOT,
                                          "model/stan/log_irt_long.stan"))

pars <- c("sigma_alpha", "pi_alpha", "sigma_xi", "s", "delta",
          "sigma_zeta", "rho_xi_zeta")
if (grepl("2pl", nm)) pars <- c(pars, "sigma_lambda")
print(summarize_fit(fit, pars = pars), digits = 3)

cat(sprintf("\nSaved: model/fits/%s.rds\n", nm))
