## 08_cv_headline_report.R -- does the HEADLINE claim survive prospective testing?
##
## THE QUESTION (MCF). The forward-CV campaign so far reports M2 vs M3, which asks whether
## per-child acceleration VARIANCE forecasts. That returned a null, and it sits awkwardly
## next to a headline claim (kappa >> 1) that had never been tested prospectively at all.
## A reader could easily conflate the two. This pairs M2k1 (per-child xi_i, kappa pinned
## to 1) against M2 (per-child xi_i, kappa free), which tests the headline directly.
##
## WHY NOT M0 vs M1, where the headline actually sits on the ladder. Neither M0 nor M1 has
## a per-child intercept, so both must predict every child's held-out administration from
## population values alone. That comparison is dominated by between-child level variance,
## and the kappa it rewards is fitted largely to CROSS-SECTIONAL age structure -- the very
## confound forward CV exists to escape. With xi_i in both models each child is anchored to
## their own level from their own training administrations, so the only thing left for the
## exponent to explain is WITHIN-child growth shape.
##
## Neither model estimates a per-child slope, so none of M3's identification problems
## apply: kappa is a population parameter and needs no depth to estimate.
##
## We expected that to make M2 vs M2k1 FLAT across the depth ladder. It is not -- it grows
## steeply, +92 per child at k=2 to +248 at k=5. The reason is visible in the bundles: the
## last-k design pins the forecast horizon at 3.0 months but lets each child's training
## window extend further BACK as k grows (earliest training age 12.0 months at k=2, 8.0 at
## k=5). With only two training administrations the free intercept xi_i can sit wherever
## best splits the difference locally and partly masks a wrong exponent; across five
## administrations spanning 8-35 months no single intercept can reconcile a kappa=1 curve
## with the observed within-child growth, so the extrapolation degrades.
##
## That is the signature of genuine functional-form misspecification rather than a fitting
## artifact: a merely less-flexible model does not get progressively WORSE as it is given
## more data. So both comparisons are depth-dependent, but they converge on opposite
## conclusions -- acceleration itself is overwhelmingly supported prospectively, while its
## between-child variation only climbs to break-even.
##
## Statistics match 02c/05: paired per-child difference in the predictive density of the
## held-out final administration, with child as the clustering unit.
##
## Usage:  Rscript studies/bayes_long/08_cv_headline_report.R
## Output: paper/cache/si_cv_headline.rds + console report

suppressPackageStartupMessages({library(dplyr)})
SUMM  <- file.path("fits", "bayes_long", "summaries")
CACHE <- file.path("paper", "cache")

## fcv4 is the main SI forward-CV setting; fcv6lk2..5 is the depth ladder's clean arm,
## where the forecast horizon is pinned at 3.0 months and only training depth varies.
TAGS <- c("fcv4", sprintf("fcv6lk%d", 2:5))

pair <- function(tag, hi, lo) {
  fh <- file.path(SUMM, sprintf("norwegian_%s_%s.rds", tag, hi))
  fl <- file.path(SUMM, sprintf("norwegian_%s_%s.rds", tag, lo))
  if (!file.exists(fh) || !file.exists(fl)) return(NULL)
  rh <- readRDS(fh); rl <- readRDS(fl)
  ## Both fits score the identical test bundle, so children pair by construction --
  ## asserted rather than assumed, as in 02c.
  idh <- sort(unique(rh$child_of_adm)); idl <- sort(unique(rl$child_of_adm))
  stopifnot(identical(idh, idl), rh$n_test_obs == rl$n_test_obs,
            length(idh) == length(rh$elpd_by_child))
  d <- rh$elpd_by_child - rl$elpd_by_child
  n <- length(d); se <- sd(d) / sqrt(n)
  hz <- rh$meta$horizon_med
  ## How far back each child's training window reaches, as the MEAN over children of
  ## their earliest training administration. Recorded because it is the mechanism behind
  ## the depth trend and is easy to get wrong: the sample-wide minimum runs 12.0 -> 8.0
  ## months across the ladder, but that is one extreme child, and the per-child mean --
  ## the quantity that describes the typical training window -- runs 26.5 -> 21.1.
  bf <- file.path("fits", "bayes_long", sprintf("bundle_norwegian_%s.rds", tag))
  first_age <- NA_real_
  if (file.exists(bf)) {
    sd0 <- readRDS(bf)$stan_data
    first_age <- mean(tapply(sd0$admin_age, sd0$admin_to_child, min))
  }
  data.frame(tag = tag, comparison = sprintf("%s - %s", hi, lo),
             n_child = n, n_test_obs = rh$n_test_obs,
             horizon_med = if (is.null(hz)) NA_real_ else hz,
             train_first_age_mean = first_age,
             diff_per_child = mean(d), diff_se = se, diff_z = mean(d) / se,
             diff_per_obs = sum(d) / rh$n_test_obs,
             pct_child_better = 100 * mean(d > 0),
             rhat_score_hi = rh$rhat_scoring, rhat_score_lo = rl$rhat_scoring,
             row.names = NULL)
}

## m2 - m2k1  = value of acceleration at all      (the HEADLINE claim)
## m3 - m2    = value of per-child acceleration   (the existing null, for contrast)
res <- bind_rows(
  bind_rows(lapply(TAGS, pair, hi = "m2", lo = "m2k1")),
  bind_rows(lapply(TAGS, pair, hi = "m3", lo = "m2")))
if (!nrow(res)) { cat("no completed pairs yet\n"); quit(save = "no") }

dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(CACHE, "si_cv_headline.rds"))

show <- function(cmp, title) {
  r <- filter(res, comparison == cmp)
  if (!nrow(r)) return(invisible(NULL))
  cat("\n===", title, "===\n")
  print(r |> transmute(
    setting = tag, n = n_child,
    horizon = ifelse(is.na(horizon_med), "-", sprintf("%.1f mo", horizon_med)),
    `dELPD/child` = sprintf("%+.2f", diff_per_child),
    SE = sprintf("%.2f", diff_se), z = sprintf("%+.1f", diff_z),
    `% kids better` = sprintf("%.1f", pct_child_better),
    `max rhat (scoring)` = sprintf("%.3f", pmax(rhat_score_hi, rhat_score_lo))) |>
    as.data.frame(), row.names = FALSE)
}
cat("Norwegian forward CV: predictive density of each child's held-out FINAL administration\n")
show("m2 - m2k1", "HEADLINE: does acceleration itself forecast?  (M2 vs M2k1, kappa free vs kappa=1)")
show("m3 - m2",   "for contrast: does PER-CHILD acceleration forecast?  (M3 vs M2)")

## Both comparisons rise with depth; what differs is where they end up.
h <- filter(res, comparison == "m2 - m2k1", grepl("^fcv6lk", tag)) |> arrange(tag)
v <- filter(res, comparison == "m3 - m2",   grepl("^fcv6lk", tag)) |> arrange(tag)
if (nrow(h) && nrow(v)) {
  cat("\nAcross the depth ladder (horizon pinned at 3.0 months; only training depth,\n")
  cat("and hence how far back the training window reaches, varies):\n")
  cat(sprintf("  M2 - M2k1 (acceleration exists):    %s\n",
              paste(sprintf("%+.0f", h$diff_per_child), collapse = " -> ")))
  cat(sprintf("  M3 - M2   (it varies between kids): %s\n",
              paste(sprintf("%+.0f", v$diff_per_child), collapse = " -> ")))
  cat("  -> both rise with depth, but the first rises AWAY from zero and the second only\n")
  cat("     up TO it. More within-child data is progressively harder to reconcile with\n")
  cat("     kappa = 1, while per-child kappa_i merely stops hurting once estimable.\n")
}
cat("\nwrote", file.path(CACHE, "si_cv_headline.rds"), "\n")
