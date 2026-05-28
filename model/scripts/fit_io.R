## Fit the input-observed (BabyView/Seedlings) model variant.
##
## Usage:
##   Rscript model/scripts/fit_io.R <variant> [dataset]
##
## Variants:
##   io_2pl_slopes  -- 2PL + per-child slopes  (default)
##   io_2pl         -- 2PL, no slopes
##   io_slopes      -- Rasch + slopes
##   io_baseline    -- Rasch, no slopes
##   io_no_freq_slopes_si_signed -- + per-child trajectory phase s_i
##   io_comp_no_freq_slopes_si_signed -- comp channel + si_signed
##   ...etc; variant grammar from variant_hyperpriors() applies.
##
## Dataset keys: must be a registered "input-observed" dataset
## (e.g., 'babyview', 'seedlings'). Defined in model/R/datasets.R.
##
## Env vars:
##   STAN_BACKEND  cmdstanr (default) | rstan
##   STAN_CHAINS, STAN_ITER, STAN_WARMUP, STAN_ADAPT_DELTA, STAN_MAX_TREEDEPTH
##   STAN_THREADS_PER_CHAIN  cmdstanr only; uses reduce_sum threading
##   STAN_INIT_FROM  cmdstanr only; warm-start tag (see fit_longitudinal.R)

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")

args    <- commandArgs(trailingOnly = TRUE)
variant <- if (length(args) >= 1) args[1] else "io_2pl_slopes"
dataset <- if (length(args) >= 2) args[2] else "babyview"

bundle    <- load_dataset_bundle(dataset)
base_data <- bundle$stan_data

out_tag <- if (dataset == "babyview") {
  variant
} else {
  sprintf("%s_%s", variant, dataset)
}

cat(sprintf("Dataset: %s (%s)  I=%d, A=%d, J=%d, V=%d, N=%d\n",
            dataset, bundle$language,
            base_data$I, base_data$A, base_data$J, base_data$V,
            base_data$N))

# Defer to shared variant_hyperpriors() registry (the "io_" prefix is
# stripped inside, so io_slopes -> slopes etc.). The optional "_anchored"
# suffix only switches the Stan file (delta_j anchoring); strip it for
# the hyperprior lookup so io_no_freq_slopes_anchored reuses the
# io_no_freq_slopes priors.
overrides <- variant_hyperpriors(sub("_anchored$", "", variant))
# Re-apply DEFAULT_PRIORS at fit time so stale priors baked into older
# bundles (e.g., the s-prior boundary bug) are overridden by config.R.
stan_data <- modifyList(modifyList(base_data, DEFAULT_PRIORS), overrides)

cat(sprintf("\n===== Fitting %s on %s =====\n", variant, dataset))
cat("Hyperprior overrides:\n"); str(overrides)

cfg <- modifyList(DEFAULT_FIT_CONFIG, list(
  chains        = as.integer(Sys.getenv("STAN_CHAINS",        unset = "4")),
  iter          = as.integer(Sys.getenv("STAN_ITER",          unset = "1000")),
  warmup        = as.integer(Sys.getenv("STAN_WARMUP",        unset = "500")),
  adapt_delta   = as.numeric(Sys.getenv("STAN_ADAPT_DELTA",   unset = "0.95")),
  max_treedepth = as.integer(Sys.getenv("STAN_MAX_TREEDEPTH", unset = "10"))
))
cat(sprintf("Stan config: chains=%d iter=%d warmup=%d adapt_delta=%.2f max_treedepth=%d\n",
            cfg$chains, cfg$iter, cfg$warmup, cfg$adapt_delta, cfg$max_treedepth))

# Stan-file dispatch. Order matters: si_signed extensions override the
# base comp/std/io files when the variant matches.
is_si_signed <- grepl("_si_signed$", variant)
has_comp_or_std <- grepl("comp|std", variant)
is_anchored  <- grepl("anchored", variant)
stan_file <- if (is_anchored) {
  # delta_j anchored to the EN longitudinal fit (input-uptake revisited).
  # Requires a bundle with delta_j_prior_mean/_sd (io_am2018, io_fmw2013).
  "model/stan/log_irt_io_anchored.stan"
} else if (is_si_signed && has_comp_or_std) {
  "model/stan/log_irt_long_io_comp_si_signed.stan"
} else if (is_si_signed) {
  "model/stan/log_irt_long_io_accel_si_signed.stan"
} else if (has_comp_or_std) {
  "model/stan/log_irt_long_io_comp.stan"
} else {
  "model/stan/log_irt_io.stan"
}
cat(sprintf("Stan model: %s\n", stan_file))

# Backend selection. cmdstanr is default (supports reduce_sum threading);
# rstan available via STAN_BACKEND=rstan as fallback. Mirrors
# fit_longitudinal.R.
backend <- Sys.getenv("STAN_BACKEND", unset = "cmdstanr")
fit_fun <- switch(backend,
                  cmdstanr = fit_variant_cmdstanr,
                  rstan    = fit_variant,
                  stop(sprintf("Unknown STAN_BACKEND: %s", backend)))
cat(sprintf("Backend: %s\n", backend))

# Optional warm-start (cmdstanr only)
init_from <- Sys.getenv("STAN_INIT_FROM", unset = "")
init_list <- NULL
if (nzchar(init_from) && backend == "cmdstanr") {
  cat(sprintf("\nWarm-start init from: %s\n", init_from))
  init_list <- tryCatch(
    build_init_from_summary(init_from, n_chains = cfg$chains),
    error = function(e) {
      message(sprintf("init build failed (%s); using default cmdstanr init",
                      conditionMessage(e)))
      NULL
    }
  )
}

if (backend == "cmdstanr") {
  fit <- fit_fun(stan_data, tag = out_tag,
                 cfg = cfg,
                 model_path = file.path(PROJECT_ROOT, stan_file),
                 init = init_list)
} else {
  fit <- fit_fun(stan_data, tag = out_tag,
                 cfg = cfg,
                 model_path = file.path(PROJECT_ROOT, stan_file))
}

pars <- c("mu_r", "sigma_r", "sigma_alpha", "pi_alpha",
          "beta_react", "reactivity_multiplier", "sigma_within",
          "s", "delta")
if (grepl("2pl",    variant)) pars <- c(pars, "sigma_lambda")
if (grepl("slopes", variant)) pars <- c(pars, "sigma_zeta")
if (is_si_signed)             pars <- c(pars, "sigma_s")
print(summarize_fit(fit, pars = pars), digits = 3)

cat(sprintf("\nSaved: fits/%s.rds\n", out_tag))
