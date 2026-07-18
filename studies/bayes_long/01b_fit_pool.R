## 01b_fit_pool.R -- fit the mega-model (m_pool.stan) on the combined cross-dataset
## bundle. Partial-pooled child variance across datasets; per-dataset means/items.
## Usage:  Rscript studies/bayes_long/01b_fit_pool.R [suffix]   (default "_a3")
suppressPackageStartupMessages({ library(cmdstanr); library(posterior); library(loo) })
SFX <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "_a3"

CHAINS  <- as.integer(Sys.getenv("STAN_CHAINS",      "4"))
WARMUP  <- as.integer(Sys.getenv("STAN_WARMUP",      "1000"))
ITER    <- as.integer(Sys.getenv("STAN_ITER",        "1000"))
THREADS <- as.integer(Sys.getenv("STAN_THREADS",     "6"))
ADELTA  <- as.numeric(Sys.getenv("STAN_ADAPT_DELTA", "0.9"))
INIT    <- Sys.getenv("STAN_INIT","")
LOO_DRAWS  <- as.integer(Sys.getenv("LOO_DRAWS",  "400"))
LOO_MAXOBS <- as.integer(Sys.getenv("LOO_MAXOBS", "500000"))

b   <- readRDS(file.path("fits","bayes_long", sprintf("bundle_pool%s.rds", SFX)))
sd0 <- b$stan_data
grainsize <- max(1L, sd0$N %/% (2L*THREADS))

dat <- list(
  N=sd0$N, grainsize=grainsize, A=sd0$A, I=sd0$I, J=sd0$J, D=sd0$D,
  aa=sd0$aa, jj=sd0$jj, y=sd0$y, admin_to_child=sd0$admin_to_child, admin_age=sd0$admin_age,
  child_ds=sd0$child_ds, item_ds=sd0$item_ds, log_H=sd0$log_H, a0=sd0$a0,
  mu_xi_prior_mean=-6, mu_xi_prior_sd=5, delta_prior_mean=0, delta_prior_sd=5,
  log_sa_prior_mean=log(1.5), log_sa_prior_sd=1,       # typical sigma_a ~ 1.5
  log_sb_prior_mean=log(5),   log_sb_prior_sd=1,       # typical sigma_b ~ 5
  s_scale_prior_sd=0.5,                                # how much datasets vary (log scale)
  tau_delta_prior_sd=5)

cat(sprintf("=== m_pool%s ===  N=%d A=%d I=%d J=%d D=%d  grainsize=%d  (%d x %d+%d, %d thr)\n",
            SFX, sd0$N, sd0$A, sd0$I, sd0$J, sd0$D, grainsize, CHAINS, WARMUP, ITER, THREADS))
mod <- cmdstan_model(file.path("studies","bayes_long","stan","m_pool.stan"),
                     cpp_options=list(stan_threads=TRUE))
SEED <- as.integer(sum(utf8ToInt(paste0("pool",SFX))) %% 2147483647L)
sargs <- list(data=dat, seed=SEED, chains=CHAINS, parallel_chains=CHAINS, threads_per_chain=THREADS,
              iter_warmup=WARMUP, iter_sampling=ITER, adapt_delta=ADELTA, refresh=100)
if (INIT != "") sargs$init <- as.numeric(INIT)
fit <- do.call(mod$sample, sargs)

dg <- fit$diagnostic_summary()
cat(sprintf("divergences: %d | max_treedepth: %d\n", sum(dg$num_divergent), sum(dg$num_max_treedepth)))

SCALARS <- c("mu_xi","delta","kappa_pop","sigma_a","sigma_b","rho_ab","tau_delta",
             "m_a","m_b","s_a","s_b","sigma_a_pop","sigma_b_pop")
summ <- fit$summary(intersect(SCALARS, fit$metadata()$stan_variables)); print(summ, n=100)

## LOO: reconstruct per-obs log_lik in R from admin_base + item_offset (same as 01_fit.R)
loo_res <- tryCatch({
  ab <- posterior::as_draws_matrix(fit$draws("admin_base"))
  io <- posterior::as_draws_matrix(fit$draws("item_offset"))
  di <- round(seq(1, nrow(ab), length.out=min(LOO_DRAWS, nrow(ab))))
  ab <- ab[di,,drop=FALSE]; io <- io[di,,drop=FALSE]
  set.seed(sum(utf8ToInt(paste0("pool",SFX))))
  oi <- if (sd0$N > LOO_MAXOBS) sort(sample.int(sd0$N, LOO_MAXOBS)) else seq_len(sd0$N)
  aa <- sd0$aa[oi]; jj <- sd0$jj[oi]; yy <- sd0$y[oi]
  ll <- matrix(0, nrow(ab), length(oi))
  for (d in seq_len(nrow(ab))) { eta <- ab[d, aa] + io[d, jj]
    ll[d, ] <- yy*eta - (pmax(eta,0) + log1p(exp(-abs(eta)))) }
  res <- loo::loo(ll, r_eff = loo::relative_eff(exp(ll), chain_id=rep(1L, nrow(ll))))
  attr(res,"n_obs_used") <- length(oi); print(res); res
}, error=function(e){ cat("LOO failed:", conditionMessage(e), "\n"); NULL })

OUT <- file.path("fits","bayes_long","summaries"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
tag <- sprintf("pool%s", SFX)
saveRDS(summ, file.path(OUT, paste0(tag,".summary.rds")))
saveRDS(fit$draws(intersect(SCALARS, fit$metadata()$stan_variables), format="df"), file.path(OUT, paste0(tag,".draws.rds")))
if (!is.null(loo_res)) saveRDS(loo_res, file.path(OUT, paste0(tag,".loo.rds")))
## per-child (efficiency/acceleration) + per-item (difficulty) exports, w/ dataset label
xi <- fit$summary("xi","median")$median; kap <- fit$summary("kappa","median")$median
saveRDS(data.frame(ii=seq_len(sd0$I), dataset=sd0$child_ds, xi=xi, kappa=kap),
        file.path(OUT, paste0(tag,"_child.rds")))
cat("saved", tag, "\n")
