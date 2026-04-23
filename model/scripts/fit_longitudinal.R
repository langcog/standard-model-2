## Fit the longitudinal log-linear IRT accumulator model.
##
## Usage:
##   Rscript model/scripts/fit_longitudinal.R <variant> [dataset]
##
## Variants (2PL on/off, slopes on/off):
##   long_baseline    -- Rasch, no per-child slopes
##   long_2pl         -- 2PL, no slopes
##   long_slopes      -- Rasch + per-child slopes
##   long_2pl_slopes  -- 2PL + per-child slopes  (default)
##
## Dataset keys are defined in model/R/datasets.R; default = "english".
## Run tag = <variant>[_<dataset>] unless dataset is "english".
##
## Examples:
##   Rscript model/scripts/fit_longitudinal.R long_2pl_slopes
##   Rscript model/scripts/fit_longitudinal.R long_2pl_slopes norwegian
##   Rscript model/scripts/fit_longitudinal.R long_slopes norwegian

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")

args       <- commandArgs(trailingOnly = TRUE)
variant    <- if (length(args) >= 1) args[1] else "long_2pl_slopes"
dataset    <- if (length(args) >= 2) args[2] else "english"

bundle_info <- get_dataset(dataset)
bundle      <- load_dataset_bundle(dataset)
base_data   <- bundle$stan_data

# Output tag: variant alone for English, variant_<key> otherwise, so
# pre-existing English fits (long_2pl_slopes.rds) don't need renaming.
out_tag <- if (dataset == "english") variant
           else sprintf("%s_%s", variant, dataset)

cat(sprintf("Dataset: %s (%s)  I=%d, A=%d, J=%d, N=%d\n",
            dataset, bundle$language,
            base_data$I, base_data$A, base_data$J, base_data$N))

overrides <- switch(variant,
  long_baseline    = list(),
  long_2pl         = list(sigma_lambda_prior_sd = 1),
  long_slopes      = list(sigma_zeta_prior_sd = 1),
  long_2pl_slopes  = list(sigma_lambda_prior_sd = 1,
                          sigma_zeta_prior_sd = 1),
  stop(sprintf("Unknown longitudinal variant: %s", variant))
)
stan_data <- modifyList(base_data, overrides)

cat(sprintf("\n===== Fitting %s on %s =====\n", variant, dataset))
cat("Hyperprior overrides:\n"); str(overrides)

cfg <- modifyList(DEFAULT_FIT_CONFIG,
                  list(chains = 4, iter = 1500, warmup = 750,
                       adapt_delta = 0.95))

fit <- fit_variant(stan_data, tag = out_tag,
                   cfg = cfg,
                   model_path = file.path(PROJECT_ROOT,
                                          "model/stan/log_irt_long.stan"))

pars <- c("sigma_alpha", "pi_alpha", "sigma_xi", "s", "delta",
          "sigma_zeta", "rho_xi_zeta")
if (grepl("2pl", variant)) pars <- c(pars, "sigma_lambda")
print(summarize_fit(fit, pars = pars), digits = 3)

cat(sprintf("\nSaved: model/fits/%s.rds\n", out_tag))
