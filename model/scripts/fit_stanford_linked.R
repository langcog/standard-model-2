## Fit the joint vocab + LWL processing model on the Stanford-linked
## subjects (Adams 2018 CDIs joined with Peekbank LWL).
##
## Usage:   Rscript model/scripts/fit_stanford_linked.R [variant]
##   variant: long_proc_slopes (default) | long_proc_baseline | long_proc_2pl_slopes ...
##   The variant grammar from helpers.R applies as usual; "long_proc_"
##   is just another prefix that gets stripped inside variant_hyperpriors().

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
  chains      = as.integer(Sys.getenv("STAN_CHAINS",      unset = "4")),
  iter        = as.integer(Sys.getenv("STAN_ITER",        unset = "1000")),
  warmup      = as.integer(Sys.getenv("STAN_WARMUP",      unset = "500")),
  adapt_delta = as.numeric(Sys.getenv("STAN_ADAPT_DELTA", unset = "0.95"))
))
cat(sprintf("Stan config: chains=%d iter=%d warmup=%d adapt_delta=%.2f\n",
            cfg$chains, cfg$iter, cfg$warmup, cfg$adapt_delta))

fit <- fit_variant(stan_data, tag = variant,
                   cfg = cfg,
                   model_path = file.path(PROJECT_ROOT,
                                          "model/stan/log_irt_long_proc.stan"))

pars <- c("sigma_alpha", "pi_alpha", "sigma_xi",
          "sigma_zeta", "sigma_rtslope",
          "s", "delta",
          "mu_rt", "mu_rtslope", "gamma_rt", "sigma_lwl",
          "rho_alpha_zeta", "rho_alpha_rtslope", "rho_zeta_rtslope")
if (grepl("2pl", variant)) pars <- c(pars, "sigma_lambda")
print(summarize_fit(fit, pars = pars), digits = 3)

cat(sprintf("\nSaved: model/fits/%s.rds\n", variant))
