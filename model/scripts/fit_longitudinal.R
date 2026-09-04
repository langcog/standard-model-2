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

# Optional sigma_r override for sensitivity analysis. Set
# STAN_SIGMA_R_OVERRIDE=<value> to refit at a different externally
# pinned input-rate variance; tag suffix _sigmaR_<value> is appended.
sigma_r_override <- Sys.getenv("STAN_SIGMA_R_OVERRIDE", unset = "")
if (nzchar(sigma_r_override)) {
  sr <- as.numeric(sigma_r_override)
  if (!is.finite(sr) || sr <= 0) stop("Invalid STAN_SIGMA_R_OVERRIDE: ", sigma_r_override)
  cat(sprintf("\nsigma_r override: bundle had %.4f; using %.4f\n",
              stan_data$sigma_r, sr))
  stan_data$sigma_r <- sr
  out_tag <- sprintf("%s_sigmaR_%s", out_tag,
                     sub("\\.", "p", sprintf("%.2f", sr)))
}
# Optional free-form tag suffix (e.g. STAN_TAG_SFX=_2k for a longer
# convergence refit), so a rerun never clobbers the CSV dir / summary of
# an earlier run at other settings. Mirrors fit_joint_io_proc_lean.R.
out_tag <- paste0(out_tag, Sys.getenv("STAN_TAG_SFX", unset = ""))

cat(sprintf("\n===== Fitting %s on %s =====\n", variant, dataset))
cat("Hyperprior overrides:\n"); str(overrides)

# Defaults for longitudinal fits; overridable via env vars so SLURM
# scripts can dial them without touching code.
# Defaults tuned for quick exploration: 1000 iter / 500 warmup gives
# ~500 sampling iter x 4 chains = 2000 effective samples, plenty for
# posterior summaries. Bump to 1500/750 for publication-grade runs.
cfg <- modifyList(DEFAULT_FIT_CONFIG, list(
  chains        = as.integer(Sys.getenv("STAN_CHAINS",        unset = "4")),
  iter          = as.integer(Sys.getenv("STAN_ITER",          unset = "1000")),
  warmup        = as.integer(Sys.getenv("STAN_WARMUP",        unset = "500")),
  adapt_delta   = as.numeric(Sys.getenv("STAN_ADAPT_DELTA",   unset = "0.95")),
  max_treedepth = as.integer(Sys.getenv("STAN_MAX_TREEDEPTH", unset = "10"))
))
cat(sprintf("Stan config: chains=%d iter=%d warmup=%d adapt_delta=%.2f max_treedepth=%d\n",
            cfg$chains, cfg$iter, cfg$warmup, cfg$adapt_delta, cfg$max_treedepth))

# LMM variants use a different Stan file (linear-in-age, no s / delta).
is_lmm <- grepl("^long_lmm", variant) || identical(sub("^long_", "", variant), "lmm") ||
          identical(sub("^long_", "", variant), "lmm_slopes")
# si_corr variants use the trivariate-LKJ (xi, zeta, s_i) Stan file.
is_si_corr <- grepl("_si_corr$", variant)
# si_signed variants use signed-normal s_i (sum-to-zero centered) +
# the (sigma_total, p_zeta) reparam. Fixes both (sigma_zeta, sigma_s)
# AND (sigma_s, delta) mixing ridges. Supersedes the half-normal
# si_reparam variant which is no longer maintained.
is_si_signed <- grepl("_si_signed$", variant)
# D' (input-on-slope): rho_xi_zeta pinned 0 + gamma_in * log_r_dev in slope.
is_dprime <- grepl("_dprime$", variant)
stan_file <- file.path(PROJECT_ROOT,
                       if (is_lmm) "model/stan/log_irt_long_lmm.stan"
                       else if (is_si_corr) "model/stan/log_irt_long_si_corr.stan"
                       else if (is_si_signed) "model/stan/log_irt_long_si_signed.stan"
                       else if (is_dprime) "model/stan/log_irt_long_dprime.stan"
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

# Optional warm-start: STAN_INIT_FROM=<source_tag> uses posterior
# medians from that source fit's summary.rds as chain initialization
# (with small Gaussian jitter per chain). Skips cold-start exploration,
# typically halving warmup requirements. Source must have been
# previously fit and extracted to fits/summaries/<source_tag>.summary.rds.
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
                 model_path = stan_file,
                 init = init_list)
} else {
  # rstan backend doesn't support init in fit_variant; warm-start
  # is cmdstanr-only.
  fit <- fit_fun(stan_data, tag = out_tag,
                 cfg = cfg,
                 model_path = stan_file)
}

pars <- c("sigma_alpha", "pi_alpha", "sigma_xi",
          "sigma_zeta", "rho_xi_zeta")
if (is_lmm) {
  pars <- c(pars, "beta_age")
} else {
  pars <- c(pars, "s", "delta")
}
if (grepl("2pl", variant)) pars <- c(pars, "sigma_lambda")
if (grepl("_si", variant)) pars <- c(pars, "sigma_s")
if (Sys.getenv("STAN_SKIP_SAVE_OBJECT", unset = "0") == "1") {
  # summarize_fit() pulls every draw (log_lik included) into memory -- the
  # same OOM save_object() hits on the big bundles. All six _2k imputed
  # refits died HERE on 2026-08-19/20 after sampling cleanly, because the
  # flag skipped save_object but not this. The slurm streams the summary
  # from the CSVs instead (cluster/sherlock/recover_from_csvs.R).
  cat("STAN_SKIP_SAVE_OBJECT=1: skipping in-R summarize_fit; recover from CSVs\n")
} else {
  print(summarize_fit(fit, pars = pars), digits = 3)
}

cat(sprintf("\nSaved: fits/%s.rds\n", out_tag))
