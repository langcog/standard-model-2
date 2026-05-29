## Fit the pooled hierarchical IO model over all four datasets.
##
## Usage: Rscript model/scripts/fit_io_pooled.R
##   env: STAN_ITER, STAN_WARMUP, STAN_CHAINS, STAN_THREADS_PER_CHAIN
##
## Output: fits/io_pooled.rds  (+ scratch on Sherlock via STANDARD_MODEL_FITS_DIR)

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
options(cmdstanr_no_ver_check = TRUE)

bundle <- readRDS(file.path(PATHS$fits_dir, "io_pooled_subset_data.rds"))
# no_freq_slopes regime: pin s, beta_c, lambda; free sigma_zeta.
ov <- variant_hyperpriors("no_freq_slopes")
stan_data <- modifyList(modifyList(bundle$stan_data, DEFAULT_PRIORS), ov)

cfg <- list(
  chains = as.integer(Sys.getenv("STAN_CHAINS",  unset = "4")),
  iter   = as.integer(Sys.getenv("STAN_ITER",    unset = "1000")),
  warmup = as.integer(Sys.getenv("STAN_WARMUP",  unset = "500")),
  threads = as.integer(Sys.getenv("STAN_THREADS_PER_CHAIN", unset = "2"))
)
# reduce_sum partition: aim for ~2 chunks per thread for TBB load-balancing
stan_data$grainsize <- max(1L, as.integer(floor(stan_data$N / (2 * cfg$threads))))
cat(sprintf("Pooled IO: N=%d, I=%d, J=%d, S=%d, V=%d | chains=%d iter=%d warmup=%d threads=%d grainsize=%d\n",
            stan_data$N, stan_data$I, stan_data$J, stan_data$S, stan_data$V,
            cfg$chains, cfg$iter, cfg$warmup, cfg$threads, stan_data$grainsize))

mod <- cmdstan_model("model/stan/log_irt_io_pooled.stan",
                     cpp_options = list(stan_threads = TRUE))
fit <- mod$sample(
  data = stan_data, chains = cfg$chains, parallel_chains = cfg$chains,
  threads_per_chain = cfg$threads,
  iter_warmup = cfg$warmup, iter_sampling = cfg$iter - cfg$warmup,
  adapt_delta = 0.95, max_treedepth = 11, refresh = 100, seed = 2026
)

out <- file.path(PATHS$fits_dir, "io_pooled.rds")
fit$save_object(out)
cat(sprintf("Saved %s\n", out))

scal <- fit$summary(c("kappa_pop","pi_alpha","sigma_r","sigma_alpha",
                       "sigma_zeta","sigma_within","kappa_study",
                       "study_input_mean"))
print(as.data.frame(scal), digits = 3)
