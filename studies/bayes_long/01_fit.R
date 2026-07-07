## 01_fit.R -- fit one (dataset, model) of the bayes_long ladder.
##
## Usage:  Rscript studies/bayes_long/01_fit.R <slug> <m0|m1|m2|m3>
## Env:    STAN_CHAINS/WARMUP/ITER/THREADS/ADAPT_DELTA, LOO_DRAWS, LOO_MAXOBS
## Output: fits/bayes_long/summaries/<slug>_<model>.{summary,draws,loo}.rds

suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(loo)
})
# paths relative to project root (slurm cd's here; run locally from repo root too)

args  <- commandArgs(trailingOnly = TRUE)
slug  <- args[1]; model <- args[2]
stopifnot(model %in% c("m0","m1","m2","m3"))

CHAINS  <- as.integer(Sys.getenv("STAN_CHAINS",      "4"))
WARMUP  <- as.integer(Sys.getenv("STAN_WARMUP",      "1000"))
ITER    <- as.integer(Sys.getenv("STAN_ITER",        "1000"))
THREADS <- as.integer(Sys.getenv("STAN_THREADS",     "4"))
ADELTA  <- as.numeric(Sys.getenv("STAN_ADAPT_DELTA", "0.9"))
LOO_DRAWS  <- as.integer(Sys.getenv("LOO_DRAWS",  "400"))     # thin draws for LOO
LOO_MAXOBS <- as.integer(Sys.getenv("LOO_MAXOBS", "500000"))  # subsample obs for LOO if larger

PRI <- list(mu_xi_prior_mean=-6, mu_xi_prior_sd=5,
            delta_prior_mean=0,  delta_prior_sd=5,
            sigma_a_prior_sd=3,  sigma_b_prior_sd=5,
            tau_delta_prior_sd=5)

STAN_FILE <- c(m0="m0_accumulator", m1="m1_acceleration",
               m2="m2_efficiency",  m3="m3_full")[model]

b   <- readRDS(file.path("fits","bayes_long", sprintf("bundle_%s.rds", slug)))
sd0 <- b$stan_data
grainsize <- max(1L, sd0$N %/% (2L*THREADS))

## per-model data list (only the fields that model declares)
dat <- list(N=sd0$N, grainsize=grainsize, A=sd0$A, J=sd0$J,
            aa=sd0$aa, jj=sd0$jj, y=sd0$y, admin_age=sd0$admin_age,
            log_H=sd0$log_H, a0=sd0$a0,
            mu_xi_prior_mean=PRI$mu_xi_prior_mean, mu_xi_prior_sd=PRI$mu_xi_prior_sd,
            tau_delta_prior_sd=PRI$tau_delta_prior_sd)
if (model %in% c("m1","m2","m3")) { dat$delta_prior_mean<-PRI$delta_prior_mean; dat$delta_prior_sd<-PRI$delta_prior_sd }
if (model %in% c("m2","m3"))      { dat$I<-sd0$I; dat$admin_to_child<-sd0$admin_to_child; dat$sigma_a_prior_sd<-PRI$sigma_a_prior_sd }
if (model == "m3")                { dat$sigma_b_prior_sd<-PRI$sigma_b_prior_sd }

cat(sprintf("=== %s / %s ===  N=%d A=%d I=%d J=%d  grainsize=%d  (%d chains x %d+%d, %d threads)\n",
            slug, model, sd0$N, sd0$A, sd0$I, sd0$J, grainsize, CHAINS, WARMUP, ITER, THREADS))

mod <- cmdstan_model(file.path("studies","bayes_long","stan", paste0(STAN_FILE, ".stan")),
                     cpp_options=list(stan_threads=TRUE))
fit <- mod$sample(data=dat, chains=CHAINS, parallel_chains=CHAINS, threads_per_chain=THREADS,
                  iter_warmup=WARMUP, iter_sampling=ITER, adapt_delta=ADELTA, refresh=200)

dg <- fit$diagnostic_summary()
cat(sprintf("divergences: %d | max_treedepth hits: %d\n",
            sum(dg$num_divergent), sum(dg$num_max_treedepth)))

SCALARS <- intersect(c("mu_xi","delta","kappa_pop","sigma_a","sigma_b","rho_ab","tau_delta"),
                     fit$metadata()$stan_variables)
summ <- fit$summary(SCALARS); print(summ)

## ---- LOO: reconstruct per-obs log_lik in R from admin_base + item_offset ----
loo_res <- tryCatch({
  ab <- posterior::as_draws_matrix(fit$draws("admin_base"))   # draws x A
  io <- posterior::as_draws_matrix(fit$draws("item_offset"))  # draws x J
  di <- round(seq(1, nrow(ab), length.out=min(LOO_DRAWS, nrow(ab))))
  ab <- ab[di,,drop=FALSE]; io <- io[di,,drop=FALSE]
  # deterministic per-dataset obs subsample -> SAME obs across M0-M3 (comparable ELPD)
  set.seed(sum(utf8ToInt(slug)))
  oi <- if (sd0$N > LOO_MAXOBS) sort(sample.int(sd0$N, LOO_MAXOBS)) else seq_len(sd0$N)
  aa <- sd0$aa[oi]; jj <- sd0$jj[oi]; yy <- sd0$y[oi]
  ll <- matrix(0, nrow(ab), length(oi))
  for (d in seq_len(nrow(ab))) {
    eta <- ab[d, aa] + io[d, jj]
    ll[d, ] <- yy*eta - (pmax(eta,0) + log1p(exp(-abs(eta))))   # stable bernoulli_logit lpmf
  }
  res <- loo::loo(ll, r_eff = loo::relative_eff(exp(ll), chain_id=rep(1L, nrow(ll))))
  attr(res, "n_obs_used") <- length(oi); attr(res, "n_draws_used") <- nrow(ab)
  print(res); res
}, error=function(e){ cat("LOO failed:", conditionMessage(e), "\n"); NULL })

OUT <- file.path("fits","bayes_long","summaries"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
tag <- sprintf("%s_%s", slug, model)
saveRDS(summ, file.path(OUT, paste0(tag,".summary.rds")))
saveRDS(fit$draws(SCALARS, format="df"), file.path(OUT, paste0(tag,".draws.rds")))
if (!is.null(loo_res)) saveRDS(loo_res, file.path(OUT, paste0(tag,".loo.rds")))
cat("saved", tag, "\n")
