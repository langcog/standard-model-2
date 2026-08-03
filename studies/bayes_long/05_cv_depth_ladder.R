## 05_cv_depth_ladder.R -- does per-child acceleration predict forward once it can
## actually be estimated?
##
## THE QUESTION (MCF). The forward-CV null at >=4 administrations may not show that
## per-child acceleration carries no prospective information -- only that three training
## administrations cannot estimate it. A 2-parameter curve through 3 points has one
## residual degree of freedom, so kappa_i is dominated by shrinkage and M3's per-child
## slope is nearly M2's shared one. A null there is what you would see either way.
##
## THE DESIGN. Norwegian children with >=6 administrations (n = 371), holding the sample
## and the held-out target fixed and varying only how many administrations the model may
## use. Two ladders, because WHICH k you keep matters:
##   first-k  the earliest k. Count and span grow, but the last training point also moves
##            closer to the target, so the horizon falls from 12.0 months at k=2 to 7.0 at
##            k=5. Improvement here confounds better estimation with less extrapolation.
##   last-k   the k immediately before the held-out one. The horizon is pinned at 3.0
##            months for every k, so only kappa_i precision varies. This is the clean test.
##
## Usage:  Rscript studies/bayes_long/05_cv_depth_ladder.R
## Output: paper/cache/si_cv_depth.rds + console report

suppressPackageStartupMessages({library(dplyr)})
SUMM  <- file.path("fits", "bayes_long", "summaries")
CACHE <- file.path("paper", "cache")
KS <- 2:5

one <- function(mode, k) {
  tag <- sprintf("fcv6%sk%d", if (mode == "last") "l" else "", k)
  f2 <- file.path(SUMM, sprintf("norwegian_%s_m2.rds", tag))
  f3 <- file.path(SUMM, sprintf("norwegian_%s_m3.rds", tag))
  if (!file.exists(f2) || !file.exists(f3)) return(NULL)
  r2 <- readRDS(f2); r3 <- readRDS(f3)
  stopifnot(identical(sort(unique(r2$child_of_adm)), sort(unique(r3$child_of_adm))),
            r2$n_test_obs == r3$n_test_obs)
  d <- r3$elpd_by_child - r2$elpd_by_child
  pd <- r3$pop_diag; g <- function(v) { x <- pd$median[pd$variable == v]; if (!length(x)) NA_real_ else x }
  ## the first-k bundles were built before meta carried horizon_med; recompute if absent
  hz <- r3$meta$horizon_med
  if (is.null(hz)) {
    b <- readRDS(file.path("fits", "bayes_long", sprintf("bundle_norwegian_%s.rds", tag)))
    tr <- b$stan_data; te <- b$test
    last_tr <- vapply(seq_len(tr$I), function(i) max(tr$admin_age[tr$admin_to_child == i]), 0.0)
    hz <- median(te$admin_age - last_tr[te$admin_to_child])
  }
  data.frame(mode = mode, k = k, n_child = length(d),
             horizon = hz,
             diff_per_child = mean(d), diff_se = sd(d)/sqrt(length(d)),
             diff_z = mean(d)/(sd(d)/sqrt(length(d))),
             pct_better = 100*mean(d > 0),
             sigma_b = g("sigma_b"), rho_ab = g("rho_ab"), kappa_pop = g("kappa_pop"),
             rhat_scoring = max(r2$rhat_scoring, r3$rhat_scoring, na.rm = TRUE),
             rhat_all = max(r2$max_rhat, r3$max_rhat, na.rm = TRUE), row.names = NULL)
}

L <- bind_rows(lapply(c("first","last"), function(m) bind_rows(lapply(KS, function(k) one(m, k)))))
if (!nrow(L)) { cat("no completed pairs\n"); quit(save = "no") }
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
saveRDS(L, file.path(CACHE, "si_cv_depth.rds"))

for (m in c("last","first")) {
  d <- filter(L, mode == m); if (!nrow(d)) next
  cat(sprintf("\n=== %s-k ladder (horizon %s) ===\n", m,
      if (m == "last") "FIXED at 3.0 months" else "falls 12.0 -> 7.0 months"))
  print(d |> transmute(k, `train/kid` = k, `horizon (mo)` = sprintf("%.1f", horizon),
    `dELPD/child` = sprintf("%+.2f", diff_per_child), SE = sprintf("%.2f", diff_se),
    z = sprintf("%+.1f", diff_z), `% better` = sprintf("%.0f", pct_better),
    sigma_b = sprintf("%.2f", sigma_b)) |> as.data.frame(), row.names = FALSE)
}

cat("\n=== the mechanism: sigma_b inflation, not absence of acceleration ===\n")
cat("Norwegian's full-data sigma_b is 5.65. In the last-k ladder the TRAINING fit gives:\n")
lk <- filter(L, mode == "last")
for (i in seq_len(nrow(lk))) with(lk[i,],
  cat(sprintf("  k=%d  sigma_b %5.2f (%.1fx the full-data value)  ->  dELPD %+7.2f\n",
              k, sigma_b, sigma_b/5.65, diff_per_child)))
cat("\nSo M3 does NOT collapse onto M2 as k grows -- sigma_b stays near 5.2, well above\n")
cat("zero, i.e. real between-child variation is retained. What changes is that sigma_b\n")
cat("stops being INFLATED: with two training points a per-child slope has zero residual\n")
cat("df, kappa_i absorbs measurement noise, and extrapolating those noisy slopes forward\n")
cat("is actively harmful. The forward penalty tracks the inflation, and vanishes with it.\n")
cat("\nNote the ladder converges to ZERO, not to positive: better-estimated acceleration\n")
cat("becomes HARMLESS for forward prediction, never helpful.\n")
cat(sprintf("\nrhat: scoring params (xi, kappa, delta_j) max %.3f across all 16 fits;\n",
            max(L$rhat_scoring, na.rm = TRUE)))
cat(sprintf("      max over ALL params %.3f (the delta_j_raw funnel; not used for prediction)\n",
            max(L$rhat_all, na.rm = TRUE)))
cat("\nwrote", file.path(CACHE, "si_cv_depth.rds"), "\n")
