## Fit the JOINT input + processing model (log_irt_long_proc_dp_joint.stan) on
## the all-dataset bundle (AM2018, FM2012, FMW2013 + BabyView, SEEDLingS).
##
## Same D'0-D'3 ladder as fit_proc_dp.R; per-study sigma_lena; LOO skipped
## (variance-partition scalars don't need it; avoids the as_cmdstan_fit OOM).
##
## Usage: Rscript model/scripts/fit_joint_io_proc.R <rung>
## Output: fits/summaries/joint_io_proc_d<rung>.{summary,draws}.rds
source("model/R/config.R"); source("model/R/helpers.R")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
rung <- as.integer(if (length(args) >= 1) args[1] else 1)
stopifnot(rung %in% 0:3)

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR", unset = PATHS$fits_dir)
bundle <- readRDS(file.path(FITS_DIR, "joint_io_proc_subset_data.rds"))
sd <- bundle$stan_data

WIDE <- 1; PIN <- 0.001
sd$gamma_in_prior_sd <- WIDE
sd$beta_xi_prior_sd  <- if (rung >= 1) WIDE else PIN
sd$beta_k0_prior_sd  <- if (rung >= 2) WIDE else PIN
sd$beta_k1_prior_sd  <- if (rung >= 3) WIDE else PIN

TAG <- sprintf("joint_io_proc_d%d", rung)
cat(sprintf("===== %s =====\n  I=%d A=%d J=%d S=%d N=%d N_lwl=%d V=%d sigma_r=%.2f\n",
            TAG, sd$I, sd$A, sd$J, sd$S, sd$N, sd$N_lwl, sd$V, sd$sigma_r))
cat(sprintf("  sigma_lena/study: %s\n", paste(sprintf("%.3f", sd$sigma_lena), collapse=", ")))
cat(sprintf("  ladder: gamma_in=%.3f beta_xi=%.3f beta_k0=%.3f beta_k1=%.3f\n",
            sd$gamma_in_prior_sd, sd$beta_xi_prior_sd, sd$beta_k0_prior_sd, sd$beta_k1_prior_sd))

chains <- as.integer(Sys.getenv("STAN_CHAINS",        unset = "4"))
warmup <- as.integer(Sys.getenv("STAN_WARMUP",        unset = "750"))
iter   <- as.integer(Sys.getenv("STAN_ITER",          unset = "750"))
adelta <- as.numeric(Sys.getenv("STAN_ADAPT_DELTA",   unset = "0.95"))
mtd    <- as.integer(Sys.getenv("STAN_MAX_TREEDEPTH", unset = "10"))

init_fun <- function() list(
  tau_s = rep(6.6, sd$S), psi_s = rep(0, sd$S),
  delta = 7, s = 0.01,
  sigma_alpha = 1, sigma_zeta = 1, sigma_rt = c(0.2, 0.3), sigma_lwl = 0.2,
  mu_c = rep(8, sd$C), tau_c = rep(1, sd$C), sigma_lambda = 0.001,
  gamma_in = 0, beta_xi = 0, beta_k0 = 0, beta_k1 = 0
)

m <- cmdstan_model("model/stan/log_irt_long_proc_dp_joint.stan")
fit <- m$sample(data = sd, chains = chains, parallel_chains = chains,
                iter_warmup = warmup, iter_sampling = iter, init = init_fun,
                adapt_delta = adelta, max_treedepth = mtd, refresh = 100, seed = 1)

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

print(as.data.frame(summ |> select(variable, mean, q5, q95, rhat)), digits = 3)
cat(sprintf("\n%s DONE\n", TAG))
