## 02c_forward_cv_report.R -- summarise the forward-CV campaign (M3 vs M2).
##
## The reviewer's ask was not just "do grouped CV" but "report differences per held-out
## administration or child, rather than only very large total ELPD values". So the
## headline here is the PAIRED PER-CHILD difference: for each child, the predictive
## density of their held-out final administration under M3 minus under M2. Its mean and
## a clustered SE across children are what we report, alongside the share of children
## M3 predicts better. Totals are shown only for completeness.
##
## Pairing note: 01c saves elpd_by_child WITHOUT names, but it does save child_of_adm,
## and elpd_by_child is tapply()'d over exactly that vector -- so the child IDs are
## recoverable as sort(unique(child_of_adm)), and M2/M3 pair correctly because both fits
## score the identical test bundle. We assert the IDs match rather than assuming it.
##
## Usage:  Rscript studies/bayes_long/02c_forward_cv_report.R [tag]
##   tag defaults to "fcv4" (the >=4-administration bundles -- the FAIR test). Pass "fcv"
##   for the original >=3 run, in which holding out the last of three administrations
##   leaves only two in training and M3 is evaluated in the regime that cannot identify a
##   per-child slope. Reading the wrong tag silently reports the other experiment.
## Output: paper/cache/si_forward_cv[_<tag>].rds + console table

suppressPackageStartupMessages({library(dplyr)})
.args <- commandArgs(trailingOnly = TRUE)
TAG   <- if (length(.args)) .args[1] else "fcv4"
SUMM  <- file.path("fits", "bayes_long", "summaries")
CACHE <- file.path("paper", "cache")
DATASETS <- c(thal = "English (Thal)", smith = "English (Smith)",
              marchman = "English (Marchman)", norwegian = "Norwegian",
              japanese = "Japanese")

one <- function(slug, label) {
  f2 <- file.path(SUMM, sprintf("%s_%s_m2.rds", slug, TAG))
  f3 <- file.path(SUMM, sprintf("%s_%s_m3.rds", slug, TAG))
  if (!file.exists(f2) || !file.exists(f3)) { cat("  pending:", slug, "\n"); return(NULL) }
  r2 <- readRDS(f2); r3 <- readRDS(f3)

  id2 <- sort(unique(r2$child_of_adm)); id3 <- sort(unique(r3$child_of_adm))
  stopifnot(identical(id2, id3),
            length(id2) == length(r2$elpd_by_child),
            length(id3) == length(r3$elpd_by_child),
            r2$n_test_obs == r3$n_test_obs)

  d      <- r3$elpd_by_child - r2$elpd_by_child      # paired, per child
  n      <- length(d)
  se     <- sd(d) / sqrt(n)                          # child is the clustering unit
  data.frame(
    slug = slug, lang = label, n_child = n,
    n_test_obs = r3$n_test_obs, n_test_adm = r3$n_test_adm,
    test_age_lo = min(r3$test_age), test_age_hi = max(r3$test_age),
    elpd_m2 = r2$elpd_total, elpd_m3 = r3$elpd_total,
    diff_total = sum(d),
    diff_per_child = mean(d), diff_se = se, diff_z = mean(d) / se,
    diff_per_obs = sum(d) / r3$n_test_obs,
    pct_child_better = 100 * mean(d > 0),
    rhat_m2 = r2$max_rhat, rhat_m3 = r3$max_rhat,
    ## rhat over the parameters that actually enter the prediction (xi, kappa, delta_j),
    ## as opposed to max over everything -- the laggard is typically delta_j_raw, a raw
    ## non-centered nuisance whose product with tau_delta is what we use. NULL for runs
    ## made before 01c recorded it.
    rhat_score_m2 = if (is.null(r2$rhat_scoring)) NA_real_ else r2$rhat_scoring,
    rhat_score_m3 = if (is.null(r3$rhat_scoring)) NA_real_ else r3$rhat_scoring,
    row.names = NULL)
}

cat(sprintf("Forward CV [%s]: M3 vs M2 on each child's held-out FINAL administration\n", TAG))
res <- bind_rows(Map(one, names(DATASETS), DATASETS))
if (!nrow(res)) { cat("no completed pairs yet\n"); quit(save = "no") }
res <- res |> mutate(lang = factor(lang, levels = unname(DATASETS))) |> arrange(lang)

dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
OUTF <- if (TAG == "fcv4") "si_forward_cv.rds" else sprintf("si_forward_cv_%s.rds", TAG)
saveRDS(res, file.path(CACHE, OUTF))

cat("\n")
print(res |> transmute(
  Dataset = lang, n = n_child,
  `test ages` = sprintf("%.0f-%.0f", test_age_lo, test_age_hi),
  `dELPD/child` = sprintf("%+.2f", diff_per_child),
  SE = sprintf("%.2f", diff_se),
  z = sprintf("%+.1f", diff_z),
  `% kids M3 better` = sprintf("%.0f", pct_child_better),
  `dELPD total` = sprintf("%+.0f", diff_total)) |> as.data.frame(), row.names = FALSE)
cat(sprintf("\nrhat  all params: %.3f | SCORING params (xi, kappa, delta_j): %s\n",
            max(c(res$rhat_m2, res$rhat_m3), na.rm = TRUE),
            if (all(is.na(c(res$rhat_score_m2, res$rhat_score_m3)))) "not recorded for this run"
            else sprintf("%.3f", max(c(res$rhat_score_m2, res$rhat_score_m3), na.rm = TRUE))))
cat("wrote", file.path(CACHE, OUTF), "\n")
