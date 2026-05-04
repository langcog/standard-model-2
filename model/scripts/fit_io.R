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
##
## Dataset keys: must be a registered "input-observed" dataset
## (e.g., 'babyview', later 'seedlings'). Defined in model/R/datasets.R.

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
# stripped inside, so io_slopes -> slopes etc.)
overrides <- variant_hyperpriors(variant)
# Re-apply DEFAULT_PRIORS at fit time so stale priors baked into older
# bundles (e.g., the s-prior boundary bug) are overridden by config.R.
stan_data <- modifyList(modifyList(base_data, DEFAULT_PRIORS), overrides)

cat(sprintf("\n===== Fitting %s on %s =====\n", variant, dataset))
cat("Hyperprior overrides:\n"); str(overrides)

cfg <- modifyList(DEFAULT_FIT_CONFIG, list(
  chains      = as.integer(Sys.getenv("STAN_CHAINS",      unset = "4")),
  iter        = as.integer(Sys.getenv("STAN_ITER",        unset = "1000")),
  warmup      = as.integer(Sys.getenv("STAN_WARMUP",      unset = "500")),
  adapt_delta = as.numeric(Sys.getenv("STAN_ADAPT_DELTA", unset = "0.95"))
))
cat(sprintf("Stan config: chains=%d iter=%d warmup=%d adapt_delta=%.2f\n",
            cfg$chains, cfg$iter, cfg$warmup, cfg$adapt_delta))

# Stan-file dispatch: comp_* variants need the comp-channel Stan file;
# everything else uses the standard io file. The two are related by
# the comp file being a strict superset; comp_* fits in the original
# file would fail to find the gamma parameters, while the comp file
# with a comp-free bundle (N_comp=0) reduces exactly to the standard
# fit. Keeping them as two files makes the toggle obvious.
stan_file <- if (grepl("comp|std", variant)) {
  "model/stan/log_irt_long_io_comp.stan"
} else {
  "model/stan/log_irt_io.stan"
}
cat(sprintf("Stan model: %s\n", stan_file))

fit <- fit_variant(stan_data, tag = out_tag,
                   cfg = cfg,
                   model_path = file.path(PROJECT_ROOT, stan_file))

pars <- c("mu_r", "sigma_r", "sigma_alpha", "pi_alpha",
          "beta_react", "reactivity_multiplier", "sigma_within",
          "s", "delta")
if (grepl("2pl",    variant)) pars <- c(pars, "sigma_lambda")
if (grepl("slopes", variant)) pars <- c(pars, "sigma_zeta")
print(summarize_fit(fit, pars = pars), digits = 3)

cat(sprintf("\nSaved: fits/%s.rds\n", out_tag))
