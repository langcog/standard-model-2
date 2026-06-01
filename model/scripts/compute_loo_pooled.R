## PSIS-LOO comparison across the pooled IO ladder
## (baseline / gamma-additive / gamma-multiplicative).
##
## With N=404k the full log_lik matrix is too big to materialise, so we
## compute log_lik in R for a fixed random subsample of observations
## (same indices across models) directly from the saved posterior draws.
## This gives an unbiased relative elpd_diff for loo_compare.
##
## CDI likelihood only — the input channel is identical across the three
## models so wouldn't change the pairwise comparison.
##
## Output: fits/loo_pooled_results.rds + console table.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(posterior); library(loo)
})

b <- readRDS(file.path(PATHS$fits_dir, "io_pooled_subset_data.rds"))
sd <- b$stan_data
fits <- list(
  baseline = readRDS(file.path(PATHS$fits_dir, "io_pooled.rds")),
  add      = readRDS(file.path(PATHS$fits_dir, "io_pooled_gamma_add.rds")),
  mult     = readRDS(file.path(PATHS$fits_dir, "io_pooled_gamma_mult.rds"))
)

# ---- fixed subsample (same idx for all models) -------------------------
N_loo <- 20000L
set.seed(2026)
loo_idx <- sample.int(sd$N, N_loo)

# Per-obs info from the bundle, indexed by subsample
sub_aa  <- sd$aa[loo_idx]                       # admin
sub_jj  <- sd$jj[loo_idx]                       # item
sub_ch  <- sd$admin_to_child[sub_aa]            # child
sub_st  <- sd$study_of_child[sub_ch]            # study
sub_cc  <- sd$cc[sub_jj]                        # class
sub_logp <- sd$log_p[sub_jj]
sub_age <- sd$admin_age[sub_aa]
sub_y   <- sd$y[loo_idx]
log_H   <- sd$log_H; a0 <- sd$a0

# ---- compute log_lik per model -----------------------------------------
log_lik_for <- function(fit, model_type) {
  d <- fit$draws(format = "df")
  pull_mat <- function(pat) as.matrix(d[, grep(pat, colnames(d), value = TRUE), drop = FALSE])
  xi          <- pull_mat("^xi\\[")
  log_r_dev   <- pull_mat("^log_r_dev\\[")
  zeta        <- pull_mat("^zeta\\[")
  study_delta <- pull_mat("^study_delta\\[")
  lambda      <- pull_mat("^lambda\\[")
  delta_j     <- pull_mat("^delta_j\\[")
  beta_c      <- pull_mat("^beta_c\\[")
  delta   <- d$delta
  s_drw   <- d$s
  gamma_d <- if (model_type == "baseline") rep(0, nrow(d)) else d$gamma

  S <- nrow(d); M <- length(loo_idx)
  ll <- matrix(NA_real_, S, M)
  for (si in seq_len(S)) {
    base_slope <- 1 + delta[si] + study_delta[si, sub_st] + zeta[si, sub_ch]
    slope_obs <- switch(model_type,
      baseline = base_slope,
      add      = base_slope + gamma_d[si] * log_r_dev[si, sub_ch],
      mult     = (1 + delta[si] + study_delta[si, sub_st] + zeta[si, sub_ch]) *
                 (1 + gamma_d[si] * log_r_dev[si, sub_ch])
    )
    ae <- pmax(sub_age - s_drw[si], 0.01)
    la <- log(ae / a0)
    item_off <- beta_c[si, sub_cc] * sub_logp - delta_j[si, sub_jj]
    eta <- lambda[si, sub_jj] *
           (xi[si, sub_ch] + log_H + slope_obs * la + item_off)
    # bernoulli_logit_lpmf: y*eta - log1p(exp(eta))
    ll[si, ] <- sub_y * eta - log1p(exp(eta))
  }
  ll
}

# ---- run, collect, compare ---------------------------------------------
cat(sprintf("Subsample N_loo=%d. Computing log_lik for 3 models ...\n", N_loo))
loos <- list()
for (m in names(fits)) {
  cat(sprintf("  %s ...\n", m))
  t0 <- Sys.time()
  ll <- log_lik_for(fits[[m]], m)
  # r_eff per obs from posterior draws (length-S column ess / S)
  r_eff <- loo::relative_eff(exp(ll), chain_id = rep(1:4, each = nrow(ll)/4))
  loos[[m]] <- loo(ll, r_eff = r_eff, cores = 4)
  cat(sprintf("    elpd_loo = %.1f  (SE %.1f)   wall=%.1fs\n",
              loos[[m]]$estimates["elpd_loo","Estimate"],
              loos[[m]]$estimates["elpd_loo","SE"],
              as.numeric(difftime(Sys.time(), t0, units="secs"))))
}

cat("\n=== loo_compare (best on top, elpd_diff vs best) ===\n")
print(loo_compare(loos))

cat("\n=== Pareto-k diagnostics ===\n")
for (m in names(loos)) {
  k <- loos[[m]]$diagnostics$pareto_k
  cat(sprintf("  %-9s : ok %4d  | bad %3d  | very bad %3d\n",
              m, sum(k <= 0.7), sum(k > 0.7 & k <= 1), sum(k > 1)))
}

out_path <- file.path(PATHS$fits_dir, "loo_pooled_results.rds")
saveRDS(list(loos = loos, loo_idx = loo_idx, N_loo = N_loo), out_path)
cat(sprintf("\nSaved %s\n", out_path))
