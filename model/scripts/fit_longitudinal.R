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
out_tag <- if (dataset == "english") {
  variant
} else {
  sprintf("%s_%s", variant, dataset)
}

cat(sprintf("Dataset: %s (%s)  I=%d, A=%d, J=%d, N=%d\n",
            dataset, bundle$language,
            base_data$I, base_data$A, base_data$J, base_data$N))

# Defer to the shared variant_hyperpriors() registry. The "long_"
# prefix is stripped inside that function so longitudinal and
# cross-sectional pipelines share the same variant grammar.
overrides <- variant_hyperpriors(variant)
# Re-apply DEFAULT_PRIORS at fit time so stale priors baked into older
# bundles (e.g., the s-prior boundary bug) are overridden by config.R.
stan_data <- modifyList(modifyList(base_data, DEFAULT_PRIORS), overrides)
# Some variants need to mutate stan_data structure (e.g., no_class
# collapses class hierarchy by overriding cc / C).
stan_data <- variant_data_overrides(stan_data, variant)

cat(sprintf("\n===== Fitting %s on %s =====\n", variant, dataset))
cat("Hyperprior overrides:\n"); str(overrides)

# Defaults for longitudinal fits; overridable via env vars so SLURM
# scripts can dial them without touching code.
# Defaults tuned for quick exploration: 1000 iter / 500 warmup gives
# ~500 sampling iter x 4 chains = 2000 effective samples, plenty for
# posterior summaries. Bump to 1500/750 for publication-grade runs.
cfg <- modifyList(DEFAULT_FIT_CONFIG, list(
  chains      = as.integer(Sys.getenv("STAN_CHAINS",      unset = "4")),
  iter        = as.integer(Sys.getenv("STAN_ITER",        unset = "1000")),
  warmup      = as.integer(Sys.getenv("STAN_WARMUP",      unset = "500")),
  adapt_delta = as.numeric(Sys.getenv("STAN_ADAPT_DELTA", unset = "0.95"))
))
cat(sprintf("Stan config: chains=%d iter=%d warmup=%d adapt_delta=%.2f\n",
            cfg$chains, cfg$iter, cfg$warmup, cfg$adapt_delta))

# LMM variants use a different Stan file (linear-in-age, no s / delta).
is_lmm <- grepl("^long_lmm", variant) || identical(sub("^long_", "", variant), "lmm") ||
          identical(sub("^long_", "", variant), "lmm_slopes")
stan_file <- file.path(PROJECT_ROOT,
                       if (is_lmm) "model/stan/log_irt_long_lmm.stan"
                       else "model/stan/log_irt_long.stan")
cat(sprintf("Stan model: %s\n", stan_file))

# Backend selection. Default is cmdstanr (faster, supports reduce_sum
# threading); rstan available via STAN_BACKEND=rstan for emergencies.
backend <- Sys.getenv("STAN_BACKEND", unset = "cmdstanr")
fit_fun <- switch(backend,
                  cmdstanr = fit_variant_cmdstanr,
                  rstan    = fit_variant,
                  stop(sprintf("Unknown STAN_BACKEND: %s", backend)))
cat(sprintf("Backend: %s\n", backend))

fit <- fit_fun(stan_data, tag = out_tag,
               cfg = cfg,
               model_path = stan_file)

pars <- c("sigma_alpha", "pi_alpha", "sigma_xi",
          "sigma_zeta", "rho_xi_zeta")
if (is_lmm) {
  pars <- c(pars, "beta_age")
} else {
  pars <- c(pars, "s", "delta")
}
if (grepl("2pl", variant)) pars <- c(pars, "sigma_lambda")
print(summarize_fit(fit, pars = pars), digits = 3)

cat(sprintf("\nSaved: model/fits/%s.rds\n", out_tag))
