## Wide-delta version of fit_io_pooled_gamma.R: same gamma variants
## (additive / multiplicative) but with delta prior widened from the
## default N(0, 0.5) to N(0, 10) so it's essentially uninformative.
## Needed for an apples-to-apples ladder comparison once the baseline
## refit confirmed the prior was deflating kappa_pop.
##
## Usage: Rscript model/scripts/fit_io_pooled_gamma_widedelta.R <add|mult>
##
## Output: fits/io_pooled_gamma_widedelta_<variant>.rds

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
options(cmdstanr_no_ver_check = TRUE)

variant <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(variant) || !variant %in% c("add", "mult"))
  stop("usage: fit_io_pooled_gamma_widedelta.R <add|mult>")

spec <- list(
  add  = list(stan = "model/stan/log_irt_io_pooled_gamma_add.stan",
              gamma_prior_sd = 2.0),
  mult = list(stan = "model/stan/log_irt_io_pooled_gamma_mult.stan",
              gamma_prior_sd = 0.5)
)[[variant]]

bundle <- readRDS(file.path(PATHS$fits_dir, "io_pooled_subset_data.rds"))
ov <- variant_hyperpriors("no_freq_slopes")
stan_data <- modifyList(modifyList(bundle$stan_data, DEFAULT_PRIORS), ov)
stan_data$gamma_prior_sd   <- spec$gamma_prior_sd
stan_data$delta_prior_mean <- 0
stan_data$delta_prior_sd   <- 10                       # widened
cat(sprintf("delta prior widened: N(0, 10) (was N(0, 0.5))\n"))

cfg <- list(
  chains  = as.integer(Sys.getenv("STAN_CHAINS",            unset = "4")),
  iter    = as.integer(Sys.getenv("STAN_ITER",              unset = "2500")),
  warmup  = as.integer(Sys.getenv("STAN_WARMUP",            unset = "1000")),
  threads = as.integer(Sys.getenv("STAN_THREADS_PER_CHAIN", unset = "4"))
)
stan_data$grainsize <- max(1L, as.integer(floor(stan_data$N / (2 * cfg$threads))))
cat(sprintf("Pooled IO gamma-%s wide-delta: N=%d I=%d J=%d S=%d V=%d | chains=%d iter=%d warmup=%d threads=%d gamma_sd=%g\n",
            variant, stan_data$N, stan_data$I, stan_data$J, stan_data$S, stan_data$V,
            cfg$chains, cfg$iter, cfg$warmup, cfg$threads, spec$gamma_prior_sd))

mod <- cmdstan_model(spec$stan, cpp_options = list(stan_threads = TRUE))
fit <- mod$sample(
  data = stan_data, chains = cfg$chains, parallel_chains = cfg$chains,
  threads_per_chain = cfg$threads,
  iter_warmup = cfg$warmup, iter_sampling = cfg$iter - cfg$warmup,
  adapt_delta = 0.95, max_treedepth = 11, refresh = 100, seed = 2026
)

out <- file.path(PATHS$fits_dir, sprintf("io_pooled_gamma_widedelta_%s.rds", variant))
fit$save_object(out)
cat(sprintf("Saved %s\n", out))

scal <- fit$summary(c("gamma","delta","kappa_pop","pi_alpha","sigma_r",
                      "sigma_alpha","sigma_zeta","sigma_within","kappa_study",
                      "study_input_mean"))
print(as.data.frame(scal), digits = 3)
