## Fit one rung of the D'0-D'3 processing regression ladder
## (log_irt_long_proc_dp.stan).
##
## Usage: Rscript model/scripts/fit_proc_dp.R <rung> [dataset_tag]
##   rung:        0,1,2,3  ->  D'0 (input only) .. D'3 (+rt1->kappa)
##   dataset_tag: "all" (default, 3 datasets) | "adams" | ...
##
## Ladder = prior SD toggles (tiny -> coefficient pinned at 0):
##   D'0 gamma_in on              D'2 + beta_k0 (rt0->kappa)
##   D'1 + beta_xi (rt0->xi)      D'3 + beta_k1 (rt1->kappa)
##
## Env: STAN_CHAINS STAN_WARMUP STAN_ITER STAN_ADAPT_DELTA STAN_MAX_TREEDEPTH
## Output: fits/summaries/proc_dp<rung>_<tag>.{summary,draws,loo}.rds

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior); library(dplyr) })

args  <- commandArgs(trailingOnly = TRUE)
rung  <- as.integer(if (length(args) >= 1) args[1] else 0)
dstag <- if (length(args) >= 2) args[2] else "all"
stopifnot(rung %in% 0:3)

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR", unset = PATHS$fits_dir)
bundle <- readRDS(file.path(FITS_DIR, sprintf("proc_dp_%s_subset_data.rds", dstag)))
sd <- bundle$stan_data

# ---- ladder: set the four regression prior SDs by rung ----
WIDE <- 1; PIN <- 0.001
sd$gamma_in_prior_sd <- WIDE                              # D'0+
sd$beta_xi_prior_sd  <- if (rung >= 1) WIDE else PIN      # D'1
sd$beta_k0_prior_sd  <- if (rung >= 2) WIDE else PIN      # D'2
sd$beta_k1_prior_sd  <- if (rung >= 3) WIDE else PIN      # D'3

TAG <- sprintf("proc_dp%d_%s", rung, dstag)
cat(sprintf("===== %s =====\n  I=%d A=%d J=%d S=%d N=%d N_lwl=%d V=%d\n",
            TAG, sd$I, sd$A, sd$J, sd$S, sd$N, sd$N_lwl, sd$V))
cat(sprintf("  ladder: gamma_in=%.3f beta_xi=%.3f beta_k0=%.3f beta_k1=%.3f\n",
            sd$gamma_in_prior_sd, sd$beta_xi_prior_sd, sd$beta_k0_prior_sd, sd$beta_k1_prior_sd))

chains  <- as.integer(Sys.getenv("STAN_CHAINS",        unset = "4"))
warmup  <- as.integer(Sys.getenv("STAN_WARMUP",        unset = "750"))
iter    <- as.integer(Sys.getenv("STAN_ITER",          unset = "750"))
adelta  <- as.numeric(Sys.getenv("STAN_ADAPT_DELTA",   unset = "0.95"))
mtd     <- as.integer(Sys.getenv("STAN_MAX_TREEDEPTH", unset = "10"))

# Sensible inits keep warmup off the lkj/bernoulli boundaries (tau_s near the
# observed log-RT level, delta near the main-model slope). Unspecified raws get
# cmdstan's default random init.
init_fun <- function() list(
  tau_s = rep(6.6, sd$S), psi_s = rep(0, sd$S),
  delta = 7, s = 0.01,   # >0 strictly: s has <lower=0>, init on the boundary fails
  sigma_alpha = 1, sigma_zeta = 1, sigma_rt = c(0.2, 0.3), sigma_lwl = 0.2,
  mu_c = rep(8, sd$C), tau_c = rep(1, sd$C), sigma_lambda = 0.001,
  gamma_in = 0, beta_xi = 0, beta_k0 = 0, beta_k1 = 0
)

m <- cmdstan_model("model/stan/log_irt_long_proc_dp.stan")
# Keep the sampler CSVs (with log_lik) so LOO is recoverable offline if the
# in-job loo() ever OOMs again.
csv_dir <- file.path(FITS_DIR, paste0("csvs_", TAG))
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
fit <- m$sample(data = sd, chains = chains, parallel_chains = chains,
                iter_warmup = warmup, iter_sampling = iter, init = init_fun,
                adapt_delta = adelta, max_treedepth = mtd, refresh = 100, seed = 1,
                output_dir = csv_dir)

cat("\n=== diagnostics ===\n")
dg <- fit$diagnostic_summary()
cat(sprintf("divergences=%d  treedepth=%d\n", sum(dg$num_divergent), sum(dg$num_max_treedepth)))

SCALARS <- c("beta_xi","beta_k0","beta_k1","gamma_in","delta","s",
             "sigma_alpha","sigma_zeta","sigma_rt0","sigma_rt1","rho_rt","sigma_lwl",
             "sigma_xi","sigma_kappa",
             "share_input_xi","share_proc_xi","share_resid_xi",
             "var_input_k","var_proc_k","var_resid_k",
             paste0("tau_s[", 1:sd$S, "]"), paste0("psi_s[", 1:sd$S, "]"))
summ <- fit$summary(SCALARS)
SUMM_DIR <- file.path(FITS_DIR, "summaries")
dir.create(SUMM_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(summ, file.path(SUMM_DIR, paste0(TAG, ".summary.rds")))
saveRDS(fit$draws(SCALARS, format = "df"), file.path(SUMM_DIR, paste0(TAG, ".draws.rds")))

# cores=1 avoids fork-copying the (N x draws) log_lik matrix -- the OOM cause at 24G.
loo_res <- tryCatch(fit$loo(cores = 1), error = function(e) { cat("LOO failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(loo_res)) {
  saveRDS(loo_res, file.path(SUMM_DIR, paste0(TAG, ".loo.rds")))
  cat(sprintf("elpd_loo = %.1f (se %.1f)  N=%d\n", loo_res$estimates["elpd_loo","Estimate"],
              loo_res$estimates["elpd_loo","SE"], nrow(loo_res$pointwise)))
} else {
  cat(sprintf("LOO not written; recover offline from CSVs in %s\n", csv_dir))
}
print(as.data.frame(summ |> select(variable, mean, q5, q95, rhat)), digits = 3)
cat(sprintf("\n%s DONE\n", TAG))
