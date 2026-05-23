## Fit the joint vocab + LWL processing model on the Stanford-linked
## subjects (Adams 2018 CDIs joined with Peekbank LWL).
##
## Usage:   Rscript model/scripts/fit_stanford_linked.R [variant]
##   variant: long_proc_slopes (default) | long_proc_baseline |
##            long_proc_no_freq_slopes | long_proc_no_freq_slopes_si_signed | ...
##   The variant grammar from helpers.R applies as usual; "long_proc_"
##   is just another prefix that gets stripped inside variant_hyperpriors().
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
variant <- if (length(args) >= 1) args[1] else "long_proc_slopes"

bundle    <- load_dataset_bundle("stanford_linked")
base_data <- bundle$stan_data

cat(sprintf("Dataset: stanford_linked (%s)\n  I=%d, A=%d, J=%d, N=%d, N_lwl=%d\n",
            bundle$language,
            base_data$I, base_data$A, base_data$J, base_data$N,
            base_data$N_lwl))

# Apply variant hyperprior overrides; the "long_proc_" prefix is
# stripped inside variant_hyperpriors().
overrides <- variant_hyperpriors(variant)
stan_data <- modifyList(modifyList(base_data, DEFAULT_PRIORS), overrides)

cat(sprintf("\n===== Fitting %s =====\n", variant))
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

# Stan-file dispatch: si_signed variants get the s_i-extended Stan file.
is_si_signed <- grepl("_si_signed$", variant)
stan_file <- if (is_si_signed) {
  "model/stan/log_irt_long_proc_si_signed.stan"
} else {
  "model/stan/log_irt_long_proc.stan"
}
cat(sprintf("Stan model: %s\n", stan_file))

# Backend selection (cmdstanr default; rstan fallback)
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
  fit <- fit_fun(stan_data, tag = variant,
                 cfg = cfg,
                 model_path = file.path(PROJECT_ROOT, stan_file),
                 init = init_list)
} else {
  fit <- fit_fun(stan_data, tag = variant,
                 cfg = cfg,
                 model_path = file.path(PROJECT_ROOT, stan_file))
}

pars <- c("sigma_alpha", "pi_alpha", "sigma_xi",
          "sigma_zeta", "sigma_rtslope",
          "s", "delta",
          "mu_rt", "mu_rtslope", "gamma_rt", "sigma_lwl",
          "rho_alpha_zeta", "rho_alpha_rtslope", "rho_zeta_rtslope")
if (grepl("2pl", variant))   pars <- c(pars, "sigma_lambda")
if (is_si_signed)            pars <- c(pars, "sigma_s")
print(summarize_fit(fit, pars = pars), digits = 3)

cat(sprintf("\nSaved: fits/%s.rds\n", variant))
