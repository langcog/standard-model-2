## Fit a gamma (input-on-slope) variant of the pooled hierarchical IO model.
##
## Usage: Rscript model/scripts/fit_io_pooled_gamma.R <add|mult>
##   env: STAN_ITER, STAN_WARMUP, STAN_CHAINS, STAN_THREADS_PER_CHAIN
##
## "add"  : slope = 1 + delta + study_delta + gamma*log_r_dev + zeta
## "mult" : slope = (1 + delta + study_delta + zeta)*(1 + gamma*log_r_dev)
##
## Output: fits/io_pooled_gamma_<variant>.rds

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
options(cmdstanr_no_ver_check = TRUE)

variant <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(variant) || !variant %in% c("add", "mult"))
  stop("usage: fit_io_pooled_gamma.R <add|mult>")

# gamma scale differs by parameterization: additive gamma is in slope/kappa
# units (~10), multiplicative gamma is a fraction (~ gamma_add / kappa).
spec <- list(
  add  = list(stan = "model/stan/log_irt_io_pooled_gamma_add.stan",
              gamma_prior_sd = 2.0),
  mult = list(stan = "model/stan/log_irt_io_pooled_gamma_mult.stan",
              gamma_prior_sd = 0.5)
)[[variant]]

bundle <- readRDS(file.path(PATHS$fits_dir, "io_pooled_subset_data.rds"))
ov <- variant_hyperpriors("no_freq_slopes")
stan_data <- modifyList(modifyList(bundle$stan_data, DEFAULT_PRIORS), ov)
stan_data$gamma_prior_sd <- spec$gamma_prior_sd

cfg <- list(
  chains = as.integer(Sys.getenv("STAN_CHAINS",  unset = "4")),
  iter   = as.integer(Sys.getenv("STAN_ITER",    unset = "1000")),
  warmup = as.integer(Sys.getenv("STAN_WARMUP",  unset = "500")),
  threads = as.integer(Sys.getenv("STAN_THREADS_PER_CHAIN", unset = "2"))
)
stan_data$grainsize <- max(1L, as.integer(floor(stan_data$N / (2 * cfg$threads))))
cat(sprintf("Pooled IO gamma-%s: N=%d I=%d J=%d S=%d V=%d | chains=%d iter=%d warmup=%d threads=%d gamma_sd=%g\n",
            variant, stan_data$N, stan_data$I, stan_data$J, stan_data$S, stan_data$V,
            cfg$chains, cfg$iter, cfg$warmup, cfg$threads, spec$gamma_prior_sd))

mod <- cmdstan_model(spec$stan, cpp_options = list(stan_threads = TRUE))
fit <- mod$sample(
  data = stan_data, chains = cfg$chains, parallel_chains = cfg$chains,
  threads_per_chain = cfg$threads,
  iter_warmup = cfg$warmup, iter_sampling = cfg$iter - cfg$warmup,
  adapt_delta = 0.95, max_treedepth = 11, refresh = 100, seed = 2026
)

out <- file.path(PATHS$fits_dir, sprintf("io_pooled_gamma_%s.rds", variant))
fit$save_object(out)
cat(sprintf("Saved %s\n", out))

scal <- fit$summary(c("gamma","kappa_pop","pi_alpha","sigma_r","sigma_alpha",
                      "sigma_zeta","sigma_within","kappa_study",
                      "study_input_mean"))
print(as.data.frame(scal), digits = 3)
