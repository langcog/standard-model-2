## Diagnostic refit: pooled hierarchical IO model with delta_j anchor
## RELAXED across all items (prior SD = 5.0 everywhere instead of 0.10
## for matched items). Tests whether the tight longitudinal-anchored
## prior is what's deflating per-study kappa.
##
## If kappa_pop rises from ~6.4 toward ~10 (matching glmer per-study),
## the anchor was the deflator. If it stays low, look elsewhere.
##
## Output: fits/io_pooled_unanchored.rds

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
options(cmdstanr_no_ver_check = TRUE)

bundle <- readRDS(file.path(PATHS$fits_dir, "io_pooled_subset_data.rds"))
ov <- variant_hyperpriors("no_freq_slopes")
stan_data <- modifyList(modifyList(bundle$stan_data, DEFAULT_PRIORS), ov)

# RELAX the delta_j anchor: 5.0 for every item (was 0.10 for matched)
stan_data$delta_j_prior_sd <- rep(5.0, stan_data$J)
cat(sprintf("Relaxed delta_j anchor: all %d items now have prior SD = 5.0\n",
            stan_data$J))

cfg <- list(
  chains  = as.integer(Sys.getenv("STAN_CHAINS",            unset = "4")),
  iter    = as.integer(Sys.getenv("STAN_ITER",              unset = "2500")),
  warmup  = as.integer(Sys.getenv("STAN_WARMUP",            unset = "1000")),
  threads = as.integer(Sys.getenv("STAN_THREADS_PER_CHAIN", unset = "4"))
)
stan_data$grainsize <- max(1L, as.integer(floor(stan_data$N / (2 * cfg$threads))))
cat(sprintf("Pooled IO unanchored: N=%d I=%d J=%d S=%d V=%d | chains=%d iter=%d warmup=%d threads=%d\n",
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

out <- file.path(PATHS$fits_dir, "io_pooled_unanchored.rds")
fit$save_object(out)
cat(sprintf("Saved %s\n", out))

scal <- fit$summary(c("kappa_pop","pi_alpha","sigma_r","sigma_alpha",
                      "sigma_zeta","sigma_within","kappa_study",
                      "study_input_mean"))
print(as.data.frame(scal), digits = 3)
