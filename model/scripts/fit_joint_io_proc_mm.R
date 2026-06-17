## Fit the JOINT io+proc MEASUREMENT-MODEL variant (log_irt_long_proc_dp_joint_mm.stan)
## on the all-dataset bundle. Differences from fit_joint_io_proc.R:
##   - input is a measurement model: sigma_r ESTIMATED (lit prior), raw per-recording
##     log input, per-instrument sigma_meas, per-study mu_r_s.
##   - frank_etal_2026 informative RT priors (passed in the bundle).
## Same D'0-D'3 ladder; LOO skipped.
##
## Usage: Rscript model/scripts/fit_joint_io_proc_mm.R <rung>
## Output: fits/summaries/joint_io_proc_mm_d<rung>.{summary,draws}.rds
source("model/R/config.R"); source("model/R/helpers.R")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
rung <- as.integer(if (length(args) >= 1) args[1] else 1)
stopifnot(rung %in% 0:3)

FITS_DIR <- Sys.getenv("STANDARD_MODEL_FITS_DIR", unset = PATHS$fits_dir)
## STAN_MM_BUNDLE overridable to test a re-prepped bundle (e.g. run-level LWL aggregation)
BUNDLE <- Sys.getenv("STAN_MM_BUNDLE", unset = "joint_io_proc_mm_subset_data.rds")
bundle <- readRDS(file.path(FITS_DIR, BUNDLE))
sd <- bundle$stan_data

PIN <- 0.001
## Reparam: priors are on PER-SD effects (a common scale across channels).
## tau_slope = common prior SD for the slope-channel effects (a_in, a_k0, a_k1);
## tau_level = prior SD for b_xi (rt0 -> xi). STAN_TAU_SLOPE is the SWEEP variable.
TAU_SLOPE <- as.numeric(Sys.getenv("STAN_TAU_SLOPE", unset = "1"))
TAU_LEVEL <- as.numeric(Sys.getenv("STAN_TAU_LEVEL", unset = "1"))
sd$a_in_prior_sd <- TAU_SLOPE
sd$b_xi_prior_sd <- if (rung >= 1) TAU_LEVEL else PIN
sd$a_k0_prior_sd <- if (rung >= 2) TAU_SLOPE else PIN
sd$a_k1_prior_sd <- if (rung >= 3) TAU_SLOPE else PIN

tau_sfx <- if (TAU_SLOPE != 1) sprintf("_tau%g", TAU_SLOPE) else ""
b_sfx   <- if (grepl("runlvl", BUNDLE)) "_runlvl" else ""
TAG <- sprintf("joint_io_proc_mm_d%d%s%s", rung, b_sfx, tau_sfx)
cat(sprintf("===== %s =====\n  I=%d A=%d J=%d S=%d N=%d N_lwl=%d V_obs=%d n_instr=%d\n",
            TAG, sd$I, sd$A, sd$J, sd$S, sd$N, sd$N_lwl, sd$V_obs, sd$n_instr))
cat(sprintf("  input MM: sigma_r~N(%.2f,%.2f) [ESTIMATED]; mu_r_s_prior=%s\n",
            sd$sigma_r_prior_mean, sd$sigma_r_prior_sd, paste(sprintf("%.2f", sd$mu_r_s_prior_mean), collapse=",")))
cat(sprintf("  RT priors: tau~N(%.2f,%.2f) psi~N(%.2f,%.2f) srt0~N(%.3f,%.2f)\n",
            sd$mu_rt_prior_mean, sd$mu_rt_prior_sd, sd$psi_prior_mean, sd$mu_rtslope_prior_sd,
            sd$sigma_rt0_prior_mean, sd$sigma_rt0_prior_sd))
cat(sprintf("  ladder (per-SD effect prior SDs): a_in=%.3f b_xi=%.3f a_k0=%.3f a_k1=%.3f  [tau_slope=%.2f tau_level=%.2f]\n",
            sd$a_in_prior_sd, sd$b_xi_prior_sd, sd$a_k0_prior_sd, sd$a_k1_prior_sd, TAU_SLOPE, TAU_LEVEL))

chains  <- as.integer(Sys.getenv("STAN_CHAINS",        unset = "4"))
warmup  <- as.integer(Sys.getenv("STAN_WARMUP",        unset = "750"))
iter    <- as.integer(Sys.getenv("STAN_ITER",          unset = "750"))
adelta  <- as.numeric(Sys.getenv("STAN_ADAPT_DELTA",   unset = "0.95"))
mtd     <- as.integer(Sys.getenv("STAN_MAX_TREEDEPTH", unset = "10"))
threads <- as.integer(Sys.getenv("STAN_THREADS_PER_CHAIN", unset = "30"))  # reduce_sum
## grainsize ~ N/(2*threads): Stan's recommended chunk size for big reduce_sum
sd$grainsize <- max(1L, as.integer(floor(sd$N / (2 * threads))))
cat(sprintf("  reduce_sum: %d chains x %d threads = %d cores; grainsize=%d\n",
            chains, threads, chains * threads, sd$grainsize))

init_fun <- function() list(
  tau_s = rep(6.84, sd$S), psi_s = rep(-0.35, sd$S),
  delta = 7, s = 0.01,
  sigma_alpha = 1, sigma_zeta = 1, sigma_rt = c(0.143, 0.26), sigma_lwl = 0.24,
  mu_c = rep(8, sd$C), tau_c = rep(1, sd$C), sigma_lambda = 0.001,
  a_in = 0, b_xi = 0, a_k0 = 0, a_k1 = 0,
  # input measurement model
  sigma_r = 0.44, mu_r_s = sd$mu_r_s_prior_mean, sigma_meas = rep(0.3, sd$n_instr)
)

m <- cmdstan_model("model/stan/log_irt_long_proc_dp_joint_mm.stan",
                   cpp_options = list(stan_threads = TRUE))
fit <- m$sample(data = sd, chains = chains, parallel_chains = chains,
                threads_per_chain = threads,
                iter_warmup = warmup, iter_sampling = iter, init = init_fun,
                adapt_delta = adelta, max_treedepth = mtd, refresh = 100, seed = 1)

cat("\n=== diagnostics ===\n")
dg <- fit$diagnostic_summary()
cat(sprintf("divergences=%d  treedepth=%d\n", sum(dg$num_divergent), sum(dg$num_max_treedepth)))

SCALARS <- c("a_in","b_xi","a_k0","a_k1","delta","s",
             "sigma_r","sigma_alpha","sigma_zeta","sigma_rt0","sigma_rt1","rho_rt","sigma_lwl","sigma_lambda",
             "sigma_xi","sigma_kappa",
             "share_input_xi","share_proc_xi","share_resid_xi",
             "var_input_k","var_proc_k","var_resid_k",
             paste0("sigma_meas[", 1:sd$n_instr, "]"),
             paste0("tau_s[", 1:sd$S, "]"), paste0("psi_s[", 1:sd$S, "]"))
summ <- fit$summary(SCALARS)
SUMM_DIR <- file.path(FITS_DIR, "summaries")
dir.create(SUMM_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(summ, file.path(SUMM_DIR, paste0(TAG, ".summary.rds")))
saveRDS(fit$draws(SCALARS, format = "df"), file.path(SUMM_DIR, paste0(TAG, ".draws.rds")))

print(as.data.frame(summ |> select(variable, mean, q5, q95, rhat)), digits = 3)
cat(sprintf("\n%s DONE\n", TAG))
