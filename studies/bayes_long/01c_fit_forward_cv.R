## 01c_fit_forward_cv.R -- forward (prospective) cross-validation for one (dataset, model).
##
## Fit on each child's EARLIER administrations; score their held-out LAST one. This is
## the test the item-level LOO cannot do: item-level LOO leaves ~all of a child's other
## responses in the training set, so it rewards a kappa_i that merely interpolates that
## child's own sittings. Here the held-out administration is at an age the model never
## saw for that child, so M3's extra per-child slope has to actually predict forward.
##
## Scoring is a proper posterior predictive density, NOT a plug-in at the medians:
##   elpd_test = sum_obs log( (1/S) sum_s p(y_obs | draw_s) )
## averaged over posterior draws of (xi_i, kappa_i, delta_j). A plug-in estimate would
## ignore posterior spread and flatter the more flexible model.
##
## Usage:  Rscript studies/bayes_long/01c_fit_forward_cv.R <slug> <m2|m3>
## Env:    STAN_CHAINS/WARMUP/ITER/THREADS/ADAPT_DELTA, CV_DRAWS (default 400)
## Output: fits/bayes_long/summaries/<slug>_fcv_<model>.rds

suppressPackageStartupMessages({library(cmdstanr); library(posterior)})
args  <- commandArgs(trailingOnly = TRUE)
slug  <- args[1]; model <- args[2]
stopifnot(model %in% c("m2", "m3"))

CHAINS  <- as.integer(Sys.getenv("STAN_CHAINS",      "4"))
WARMUP  <- as.integer(Sys.getenv("STAN_WARMUP",      "1000"))
ITER    <- as.integer(Sys.getenv("STAN_ITER",        "1000"))
THREADS <- as.integer(Sys.getenv("STAN_THREADS",     "4"))
ADELTA  <- as.numeric(Sys.getenv("STAN_ADAPT_DELTA", "0.9"))
CV_DRAWS <- as.integer(Sys.getenv("CV_DRAWS", "400"))

PRI <- list(mu_xi_prior_mean=-6, mu_xi_prior_sd=5,
            delta_prior_mean=0,  delta_prior_sd=5,
            sigma_a_prior_sd=3,  sigma_b_prior_sd=5,
            tau_delta_prior_sd=5)
STAN_FILE <- c(m2 = "m2_efficiency", m3 = "m3_full")[model]

b  <- readRDS(file.path("fits", "bayes_long", sprintf("bundle_%s_fcv.rds", slug)))
tr <- b$stan_data; te <- b$test
grainsize <- max(1L, tr$N %/% (2L * THREADS))

dat <- list(N=tr$N, grainsize=grainsize, A=tr$A, I=tr$I, J=tr$J,
            aa=tr$aa, jj=tr$jj, y=tr$y,
            admin_to_child=tr$admin_to_child, admin_age=tr$admin_age,
            log_H=tr$log_H, a0=tr$a0,
            mu_xi_prior_mean=PRI$mu_xi_prior_mean, mu_xi_prior_sd=PRI$mu_xi_prior_sd,
            delta_prior_mean=PRI$delta_prior_mean, delta_prior_sd=PRI$delta_prior_sd,
            sigma_a_prior_sd=PRI$sigma_a_prior_sd,
            tau_delta_prior_sd=PRI$tau_delta_prior_sd)
if (model == "m3") dat$sigma_b_prior_sd <- PRI$sigma_b_prior_sd

SEED <- as.integer(sum(utf8ToInt(paste0(slug, "_fcv_", model))) %% 2147483647L)
cat(sprintf("=== %s fcv %s ===  train N=%d A=%d I=%d J=%d | test N=%d A=%d  (%d x %d+%d, %d thr)\n",
            slug, model, tr$N, tr$A, tr$I, tr$J, te$N, te$A, CHAINS, WARMUP, ITER, THREADS))

mod <- cmdstan_model(file.path("studies", "bayes_long", "stan", paste0(STAN_FILE, ".stan")),
                     cpp_options = list(stan_threads = TRUE))
fit <- mod$sample(data = dat, seed = SEED, chains = CHAINS, parallel_chains = CHAINS,
                  threads_per_chain = THREADS, iter_warmup = WARMUP, iter_sampling = ITER,
                  adapt_delta = ADELTA, refresh = 100)
dg <- fit$diagnostic_summary(quiet = TRUE)
cat(sprintf("divergences: %d | max_treedepth: %d | max rhat: %.3f\n",
            sum(dg$num_divergent), sum(dg$num_max_treedepth),
            max(fit$summary()$rhat, na.rm = TRUE)))

## ---- posterior draws needed for forward scoring -------------------------
## Child indices were preserved by the split, so xi[i]/kappa[i] line up with the test
## administrations' admin_to_child, and delta_j[j] with the test jj.
gd <- function(v) posterior::as_draws_matrix(fit$draws(v))
xi_d <- gd("xi"); dj_d <- gd("delta_j")
kap_d <- if (model == "m3") gd("kappa") else {
  ## m2 has a single population kappa: broadcast it to every child
  kp <- as.numeric(posterior::as_draws_matrix(fit$draws("kappa_pop")))
  matrix(kp, nrow = length(kp), ncol = tr$I)
}
S  <- nrow(xi_d)
sel <- if (S > CV_DRAWS) round(seq(1, S, length.out = CV_DRAWS)) else seq_len(S)
xi_d <- xi_d[sel, , drop = FALSE]; kap_d <- kap_d[sel, , drop = FALSE]; dj_d <- dj_d[sel, , drop = FALSE]
S <- length(sel)

## per-test-administration linear predictor offset: xi_i + kappa_i * log(age/a0) + log_H
child_of_te <- te$admin_to_child
logt_te     <- log(pmax(te$admin_age, 0.01) / tr$a0)
adm_base_d  <- xi_d[, child_of_te, drop = FALSE] +
               kap_d[, child_of_te, drop = FALSE] * rep(logt_te, each = S) + tr$log_H

## ---- chunked posterior predictive density -------------------------------
## Full S x N_test matrix would be ~2 GB; accumulate log-mean-exp in chunks instead.
logmeanexp_cols <- if (requireNamespace("matrixStats", quietly = TRUE)) {
  function(M) { m <- matrixStats::colMaxs(M)
                m + log(matrixStats::colMeans2(exp(sweep(M, 2, m, "-")))) }
} else {
  function(M) apply(M, 2, function(v) { m <- max(v); m + log(mean(exp(v - m))) })
}

CH <- 40000L; ll <- numeric(te$N)
for (start in seq(1L, te$N, by = CH)) {
  ix <- start:min(start + CH - 1L, te$N)
  eta <- adm_base_d[, te$aa[ix], drop = FALSE] - dj_d[, te$jj[ix], drop = FALSE]
  yv  <- rep(te$y[ix], each = S)
  ## log Bernoulli-logit: y*eta - log1p(exp(eta)), numerically stable
  lp  <- yv * eta - (pmax(eta, 0) + log1p(exp(-abs(eta))))
  ll[ix] <- logmeanexp_cols(lp)
  rm(eta, lp, yv)
}

## ---- report per observation, per administration, and per child ----------
## Raw totals are dominated by N; the interpretable units are per held-out
## administration and per child, which is also what the reviewer asked for.
adm_of_obs <- te$aa
elpd_by_adm   <- tapply(ll, adm_of_obs, sum)
child_of_adm  <- te$admin_to_child
elpd_by_child <- tapply(as.numeric(elpd_by_adm), child_of_adm[as.integer(names(elpd_by_adm))], sum)

res <- list(
  slug = slug, model = model, n_test_obs = te$N, n_test_adm = te$A,
  n_child = length(elpd_by_child),
  elpd_total   = sum(ll),
  elpd_per_obs = sum(ll) / te$N,
  elpd_per_adm = sum(ll) / te$A,
  elpd_by_adm   = as.numeric(elpd_by_adm),
  elpd_by_child = as.numeric(elpd_by_child),
  child_of_adm  = as.integer(child_of_adm[as.integer(names(elpd_by_adm))]),
  test_age      = te$admin_age,
  max_rhat = max(fit$summary()$rhat, na.rm = TRUE),
  divergences = sum(dg$num_divergent), draws_used = S, meta = b$meta)

OUT <- file.path("fits", "bayes_long", "summaries")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(OUT, sprintf("%s_fcv_%s.rds", slug, model)))
cat(sprintf("elpd_test=%.1f | per obs=%.5f | per held-out admin=%.2f | %d admins, %d children\n",
            res$elpd_total, res$elpd_per_obs, res$elpd_per_adm, te$A, res$n_child))
cat("wrote", file.path(OUT, sprintf("%s_fcv_%s.rds", slug, model)), "\n")
