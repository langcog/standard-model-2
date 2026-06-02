## Diagnostic refit: pooled hierarchical IO model with delta prior
## widened from the default N(0, 0.5) to N(0, 10) so the prior is
## essentially uninformative and the data fully determines delta.
##
## Diagnoses why the per-kid posterior slope means (~10) don't match
## the reported kappa_pop (~6.4): the tight default delta prior was
## shrinking delta toward zero, with zeta_i posteriors absorbing the
## ~+3 per-study mean shift to keep per-kid trajectories honest.
##
## Output: fits/io_pooled_widedelta.rds

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
options(cmdstanr_no_ver_check = TRUE)

bundle <- readRDS(file.path(PATHS$fits_dir, "io_pooled_subset_data.rds"))
ov <- variant_hyperpriors("no_freq_slopes")
stan_data <- modifyList(modifyList(bundle$stan_data, DEFAULT_PRIORS), ov)

# Widen delta prior: default was N(0, 0.5), now N(0, 10) -- essentially
# uninformative over the plausible kappa range [1, 15].
stan_data$delta_prior_mean <- 0
stan_data$delta_prior_sd   <- 10
cat(sprintf("delta prior widened: N(%.1f, %.1f) (was N(0, 0.5))\n",
            stan_data$delta_prior_mean, stan_data$delta_prior_sd))

cfg <- list(
  chains  = as.integer(Sys.getenv("STAN_CHAINS",            unset = "4")),
  iter    = as.integer(Sys.getenv("STAN_ITER",              unset = "2500")),
  warmup  = as.integer(Sys.getenv("STAN_WARMUP",            unset = "1000")),
  threads = as.integer(Sys.getenv("STAN_THREADS_PER_CHAIN", unset = "4"))
)
stan_data$grainsize <- max(1L, as.integer(floor(stan_data$N / (2 * cfg$threads))))
cat(sprintf("Pooled IO wide-delta: N=%d I=%d J=%d S=%d V=%d | chains=%d iter=%d warmup=%d threads=%d\n",
            stan_data$N, stan_data$I, stan_data$J, stan_data$S, stan_data$V,
            cfg$chains, cfg$iter, cfg$warmup, cfg$threads))

mod <- cmdstan_model("model/stan/log_irt_io_pooled.stan",
                     cpp_options = list(stan_threads = TRUE))
fit <- mod$sample(
  data = stan_data, chains = cfg$chains, parallel_chains = cfg$chains,
  threads_per_chain = cfg$threads,
  iter_warmup = cfg$warmup, iter_sampling = cfg$iter - cfg$warmup,
  adapt_delta = 0.95, max_treedepth = 11, refresh = 100, seed = 2026
)

out <- file.path(PATHS$fits_dir, "io_pooled_widedelta.rds")
fit$save_object(out)
cat(sprintf("Saved %s\n", out))

scal <- fit$summary(c("delta","kappa_pop","pi_alpha","sigma_r","sigma_alpha",
                      "sigma_zeta","sigma_within","kappa_study",
                      "study_input_mean"))
print(as.data.frame(scal), digits = 3)
