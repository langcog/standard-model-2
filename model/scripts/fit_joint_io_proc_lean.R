## Fit the JOINT io+proc LEAN model (log_irt_long_proc_dp_joint_lean.stan):
##   - RT measurement model LEVEL-ONLY: per-child rt0 + per-study tau_s + ONE global psi
##     (drops per-child rt1 slope and per-study psi_s) = the glmer detrend, done jointly.
##   - RAW coefficients; common per-SD-scale priors via FIXED reference SDs (no funnel).
## Ladder: D'0 {gamma_in} -> D'1 +beta_xi -> D'2 +beta_k0  (no D'3 -- no rt1).
##
## Usage: Rscript model/scripts/fit_joint_io_proc_lean.R <rung 0-2>
##   STAN_TAU_SLOPE (sweep var, default 1) -> per-SD prior SD for gamma_in/beta_k0
##   STAN_TAU_LEVEL (default 1)            -> per-SD prior SD for beta_xi
## Output: fits/summaries/joint_io_proc_lean_d<rung>[_tau<t>].{summary,draws}.rds
source("model/R/config.R"); source("model/R/helpers.R")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
rung <- as.integer(if (length(args) >= 1) args[1] else 2)
stopifnot(rung %in% 0:2)

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR", unset = PATHS$fits_dir)
BUNDLE   <- Sys.getenv("STAN_MM_BUNDLE", unset = "joint_io_proc_mm_subset_data.rds")
bundle   <- readRDS(file.path(FITS_DIR, BUNDLE))
sd <- bundle$stan_data

PIN <- 0.001
TAU_SLOPE <- as.numeric(Sys.getenv("STAN_TAU_SLOPE", unset = "1"))
TAU_LEVEL <- as.numeric(Sys.getenv("STAN_TAU_LEVEL", unset = "1"))
## Common per-SD-scale priors via FIXED reference SDs (the a-priori channel SDs): a coef's
## prior SD = tau / sigma_ref, so coef*predictor ~ tau per 1 (reference) SD. Fixed ref ->
## raw coef in the likelihood (term vanishes as sigma->0) -> no divide-by-estimated-SD funnel.
SIGMA_R_REF   <- sd$sigma_r_prior_mean      # ~0.44
SIGMA_RT0_REF <- sd$sigma_rt0_prior_mean    # ~0.143
sd$gamma_in_prior_sd <- TAU_SLOPE / SIGMA_R_REF
sd$beta_xi_prior_sd  <- if (rung >= 1) TAU_LEVEL / SIGMA_RT0_REF else PIN
sd$beta_k0_prior_sd  <- if (rung >= 2) TAU_SLOPE / SIGMA_RT0_REF else PIN

## WIDEN the population acceleration prior. The bundle default is N(0,0.5) (a stale
## cross-sectional default); it shrinks delta on weak/subsample data and was the documented
## past backfire (experiments.md ~entry 33: delta 9.6->5.4, zeta compensated). Required for
## valid subsample verification fits. Overridable via STAN_DELTA_SD.
sd$delta_prior_sd <- as.numeric(Sys.getenv("STAN_DELTA_SD", unset = "10"))

tau_sfx <- if (TAU_SLOPE != 1) sprintf("_tau%g", TAU_SLOPE) else ""
b_sfx   <- if (grepl("inputmulti", BUNDLE)) "_imulti" else if (grepl("bilingual", BUNDLE)) "_bi" else ""
STAN_MODEL <- Sys.getenv("STAN_MODEL", unset = "model/stan/log_irt_long_proc_dp_joint_lean.stan")
m_sfx   <- if (grepl("rasch", STAN_MODEL)) "_rasch" else ""
sub_sfx <- if (grepl("sub[0-9]", BUNDLE)) paste0("_", sub(".*(sub[0-9][^.]*)\\.rds$", "\\1", BUNDLE)) else ""
TAG <- sprintf("joint_io_proc_lean%s_d%d%s%s%s", m_sfx, rung, b_sfx, tau_sfx, sub_sfx)
cat(sprintf("===== %s =====\n  I=%d A=%d J=%d S=%d N=%d N_lwl=%d V_obs=%d\n",
            TAG, sd$I, sd$A, sd$J, sd$S, sd$N, sd$N_lwl, sd$V_obs))
cat(sprintf("  refs: sigma_r_ref=%.3f sigma_rt0_ref=%.3f -> raw prior SDs: gamma_in=%.2f beta_xi=%.2f beta_k0=%.2f  [tau_slope=%.2f tau_level=%.2f]\n",
            SIGMA_R_REF, SIGMA_RT0_REF, sd$gamma_in_prior_sd, sd$beta_xi_prior_sd, sd$beta_k0_prior_sd, TAU_SLOPE, TAU_LEVEL))

chains  <- as.integer(Sys.getenv("STAN_CHAINS",        unset = "4"))
warmup  <- as.integer(Sys.getenv("STAN_WARMUP",        unset = "750"))
iter    <- as.integer(Sys.getenv("STAN_ITER",          unset = "750"))
adelta  <- as.numeric(Sys.getenv("STAN_ADAPT_DELTA",   unset = "0.95"))
mtd     <- as.integer(Sys.getenv("STAN_MAX_TREEDEPTH", unset = "10"))
threads <- as.integer(Sys.getenv("STAN_THREADS_PER_CHAIN", unset = "30"))
sd$grainsize <- max(1L, as.integer(floor(sd$N / (2 * threads))))
cat(sprintf("  reduce_sum: %d chains x %d threads = %d cores; grainsize=%d\n",
            chains, threads, chains * threads, sd$grainsize))

m <- cmdstan_model(STAN_MODEL, cpp_options = list(stan_threads = TRUE))
MPARS <- names(m$variables()$parameters)                # model-specific param names
init_full <- list(
  tau_s = rep(6.84, sd$S), psi = -0.35,
  delta = 7, s = 0.01,
  sigma_alpha = 1, sigma_zeta = 1, sigma_rt0 = 0.143, sigma_lwl = 0.24,
  mu_c = rep(8, sd$C), tau_c = rep(1, sd$C), sigma_lambda = 0.001,
  gamma_in = 0, beta_xi = 0, beta_k0 = 0,
  sigma_r = 0.44, mu_r_s = sd$mu_r_s_prior_mean, sigma_meas = rep(0.3, sd$n_instr)
)
init_fun <- function() init_full[intersect(names(init_full), MPARS)]   # drop s/sigma_lambda for rasch
fit <- m$sample(data = sd, chains = chains, parallel_chains = chains,
                threads_per_chain = threads,
                iter_warmup = warmup, iter_sampling = iter, init = init_fun,
                adapt_delta = adelta, max_treedepth = mtd, refresh = 100, seed = 1)

cat("\n=== diagnostics ===\n")
dg <- fit$diagnostic_summary()
cat(sprintf("divergences=%d  treedepth=%d\n", sum(dg$num_divergent), sum(dg$num_max_treedepth)))

SCALARS <- c("gamma_in","beta_xi","beta_k0","eff_input_k","eff_proc_xi","eff_proc_k",
             "delta","s","sigma_r","sigma_alpha","sigma_zeta","sigma_rt0","sigma_lwl","sigma_lambda","psi",
             "sigma_xi","sigma_kappa",
             "share_input_xi","share_proc_xi","share_resid_xi",
             "var_input_k","var_proc_k","var_resid_k",
             paste0("sigma_meas[", 1:sd$n_instr, "]"), paste0("tau_s[", 1:sd$S, "]"))
SCALARS <- SCALARS[sub("\\[.*", "", SCALARS) %in% fit$metadata()$stan_variables]   # drop missing (s, sigma_lambda)
summ <- fit$summary(SCALARS)
SUMM_DIR <- file.path(FITS_DIR, "summaries")
dir.create(SUMM_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(summ, file.path(SUMM_DIR, paste0(TAG, ".summary.rds")))
saveRDS(fit$draws(SCALARS, format = "df"), file.path(SUMM_DIR, paste0(TAG, ".draws.rds")))

print(as.data.frame(summ |> select(variable, mean, q5, q95, rhat)), digits = 3)
cat(sprintf("\n%s DONE\n", TAG))
