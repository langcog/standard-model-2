## Diagnostic fits: 2x2 design varying whether s and delta are fixed.
##
## v  ariant         | delta       | s
## ----------------|-------------|-------------
## baseline        | free        | free
## fix_delta       | fixed at 0  | free
## fix_s           | free        | fixed at 2
## both_fixed      | fixed at 0  | fixed at 2
##
## Uses the SAME subsample of Wordbank fit earlier (so LOO-CV comparisons
## are valid). Each fit: 4 chains x 2000 iter (1000 warmup). Runs variants
## sequentially.

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(posterior)
})
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

here <- function(...) file.path("/Users/mcfrank/Projects/standard_model_2", ...)

## Use the same subsample saved during the first Wordbank run
dat_pack <- readRDS(here("model/wordbank_data_subset.rds"))
base_data <- dat_pack$data

## Hyperprior configurations for the 2x2 design
configs <- list(
  baseline    = list(s_prior_mean = 4.5, s_prior_sd = 2,
                     delta_prior_mean = 0, delta_prior_sd = 0.5),
  fix_delta   = list(s_prior_mean = 4.5, s_prior_sd = 2,
                     delta_prior_mean = 0, delta_prior_sd = 0.001),
  fix_s       = list(s_prior_mean = 2,   s_prior_sd = 0.001,
                     delta_prior_mean = 0, delta_prior_sd = 0.5),
  both_fixed  = list(s_prior_mean = 2,   s_prior_sd = 0.001,
                     delta_prior_mean = 0, delta_prior_sd = 0.001)
)

fits <- list()
for (nm in names(configs)) {
  cat(sprintf("\n===== Fitting variant: %s =====\n", nm))
  cfg <- configs[[nm]]
  stan_data <- modifyList(base_data, cfg)

  t0 <- Sys.time()
  fits[[nm]] <- stan(
    file    = here("model/log_irt_v2.stan"),
    data    = stan_data,
    chains  = 4, iter = 2000, warmup = 1000,
    seed    = 20250420,
    control = list(adapt_delta = 0.95, max_treedepth = 11)
  )
  cat(sprintf("%s: sampling time %.1f min\n",
              nm, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  scalar_pars <- c("sigma_alpha", "s", "delta", "pi_alpha", "sigma_xi")
  print(summary(fits[[nm]], pars = scalar_pars)$summary[
          , c("mean", "2.5%", "50%", "97.5%", "n_eff", "Rhat")])
  saveRDS(fits[[nm]], here(sprintf("model/wordbank_fit_%s.rds", nm)))
}

## LOO comparison
suppressPackageStartupMessages(library(loo))
log_liks <- lapply(fits, function(f) extract_log_lik(f, merge_chains = FALSE))
r_effs   <- lapply(log_liks, function(ll) relative_eff(exp(ll)))
loos     <- mapply(function(ll, re) loo(ll, r_eff = re), log_liks, r_effs,
                   SIMPLIFY = FALSE)

cat("\n===== LOO-CV comparison =====\n")
print(loo_compare(loos))

saveRDS(loos, here("model/wordbank_loos.rds"))
cat("\nDone.\n")
