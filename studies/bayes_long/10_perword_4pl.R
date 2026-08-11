## 10_perword_4pl.R -- one estimator, both systems.
##
## THE OBJECTION. Children are scored by a binary checklist and modelled with an IRT
## accumulator; LMs are scored by surprisal and modelled with a four-parameter logistic.
## Two different functional forms were chosen for the two systems, so part of the
## qualitative contrast may follow from that choice rather than from the learners. The
## companion worry is that the two kappas are not the same unit.
##
## THE TEST. Drop the IRT model on the child side and run the LMs' OWN estimator --
## the same four_pl_sc() applied to LM surprisal -- on child production curves. Same
## functional form, same fitting routine, same per-word unit, on both systems.
##
## WHY THIS IS A FAIR UNIT. The 4-PL s(x) = lo + (up-lo)/(1+exp((x-mid)/sc)) implies that
## p = (s-lo)/(up-lo), the fraction of the way still to go, has logit(p) = -(x-mid)/sc.
## So kappa = 0.434/sc is the slope of a LOGIT with respect to natural-log input on BOTH
## sides -- not, as is sometimes read, the slope of a bounded fraction. For the LM, p is
## the fraction of a word's surprisal reduction still outstanding; for the child, the
## fraction of mastery still outstanding, 1 - P(produce). Feeding 1 - P(produce) to the
## same routine therefore reproduces the LM computation exactly.
##
## THE ATTENUATION PROBLEM, AND WHY THE TEST IS STILL CONSERVATIVE. An LM ladder is one
## learner, so its per-word curve is a WITHIN-learner trajectory. A child's per-word curve
## can only be formed by pooling across children, and a mixture of logistic curves with
## different crossing points is flatter than any of its components. The pooled estimator
## must therefore UNDERSTATE the child's kappa. It cannot overstate it. So a child-side
## number that still greatly exceeds the LM side is a conservative result.
##
## We do not leave the size of that bias unknown. The calibration step below simulates
## responses from the fitted M3 -- known xi_i, kappa_i, delta_j -- and runs the identical
## pipeline over them. Whatever kappa comes back, against a known truth, is the estimator's
## attenuation factor, and it is what licenses reading the observed number.
##
## Usage:  Rscript studies/bayes_long/10_perword_4pl.R [slug_a3 ...]
## Output: paper/cache/si_perword_4pl.rds + console report

suppressPackageStartupMessages({library(dplyr); library(here)})
BL <- here("fits", "bayes_long"); SUMM <- file.path(BL, "summaries")
CACHE <- here("paper", "cache")
set.seed(20260811)

## ---- the LMs' estimator, copied verbatim from paper/build_cache_short.R ----------
four_pl_sc <- function(x, y) tryCatch(suppressWarnings({
  f <- nls(y ~ lo + (up - lo) / (1 + exp((x - mid) / sc)),
           start = list(up = max(y), lo = min(y), mid = mean(x), sc = 0.5),
           lower = c(up = min(y)-5, lo = min(y)-5, mid = min(x)-3, sc = 1e-3),
           upper = c(up = max(y)+5, lo = max(y)+5, mid = max(x)+3, sc = 50),
           algorithm = "port", control = nls.control(maxiter = 200, warnOnly = TRUE))
  c(sc = unname(coef(f)["sc"]), rng = unname(coef(f)["up"] - coef(f)["lo"]))
}), error = function(e) c(sc = NA_real_, rng = NA_real_))

## same retention filter the LM pipeline applies
keep_ok <- function(sc, rng, rng_min) is.finite(sc) & sc > 0.01 & sc < 10 & rng > rng_min

## Per-word kappa from a set of (x = log10 input, y = fraction still to go) curves.
## rng_min is 1 nat for LM surprisal; for a child curve y lives in [0,1], so the
## equivalent "did this word actually move" floor is a fraction of the unit interval.
perword_kappa <- function(df, rng_min) {
  df |> group_by(word) |>
    group_modify(~{ p <- four_pl_sc(.x$x, .x$y)
                    tibble(sc = p["sc"], rng = p["rng"]) }) |> ungroup() |>
    filter(keep_ok(sc, rng, rng_min)) |> mutate(kappa = 0.434 / sc)
}

## ---- child curves: P(not produced) by word and age, pooled over children ---------
## Ages are binned to whole months, the resolution CDI administrations actually carry.
child_curves <- function(slug) {
  b <- readRDS(file.path(BL, sprintf("bundle_%s.rds", slug)))
  sd0 <- b$stan_data
  data.frame(word = sd0$jj, age = round(sd0$admin_age[sd0$aa]), y01 = sd0$y) |>
    group_by(word, age) |>
    summarise(n = n(), p = mean(y01), .groups = "drop") |>
    filter(n >= 5) |>                               # a proportion from <5 kids is noise
    mutate(x = log10(age), y = 1 - p)               # fraction of mastery still outstanding
}

## ---- calibration: simulate from the fitted M3, then run the same pipeline --------
## If the pipeline returns the kappa that was simulated in, the pooled estimator is
## unbiased; whatever it returns instead is the attenuation the real estimate carries.
simulate_curves <- function(slug) {
  b <- readRDS(file.path(BL, sprintf("bundle_%s.rds", slug))); sd0 <- b$stan_data
  ch <- read.csv(file.path(SUMM, sprintf("%s_m3_child.csv", slug)))
  ps <- read.csv(file.path(SUMM, sprintf("%s_m3_psi.csv",  slug)))
  xi <- ch$xi_median[order(ch$ii)]; kap <- ch$kappa_median[order(ch$ii)]
  dj <- ps$delta_j[order(ps$jj)]
  ch_of_obs <- sd0$admin_to_child[sd0$aa]
  age_obs   <- sd0$admin_age[sd0$aa]
  eta <- xi[ch_of_obs] + sd0$log_H + kap[ch_of_obs] * log(age_obs / sd0$a0) - dj[sd0$jj]
  ysim <- rbinom(length(eta), 1, plogis(eta))
  list(truth = median(kap),
       curves = data.frame(word = sd0$jj, age = round(age_obs), y01 = ysim) |>
         group_by(word, age) |>
         summarise(n = n(), p = mean(y01), .groups = "drop") |>
         filter(n >= 5) |> mutate(x = log10(age), y = 1 - p))
}

slugs <- commandArgs(trailingOnly = TRUE)
if (!length(slugs)) slugs <- c("thal_a3", "marchman_a3")

RNG_MIN_CHILD <- 0.25    # the word's production fraction must move over >=1/4 of [0,1]
res <- list()
for (slug in slugs) {
  obs <- child_curves(slug)
  ko  <- perword_kappa(obs, RNG_MIN_CHILD)
  sim <- simulate_curves(slug)
  ks  <- perword_kappa(sim$curves, RNG_MIN_CHILD)
  ## the IRT kappa this dataset reports, for the head-to-head
  irt <- as.data.frame(readRDS(file.path(SUMM, sprintf("%s_m3.summary.rds", slug))))
  res[[slug]] <- list(observed = ko, simulated = ks, truth = sim$truth,
                      n_words_obs = n_distinct(obs$word), n_words_fit = nrow(ko),
                      kappa_obs = median(ko$kappa),
                      kappa_iqr = unname(quantile(ko$kappa, c(.25, .75))),
                      recovery  = median(ks$kappa) / sim$truth,
                      kappa_corrected = median(ko$kappa) / (median(ks$kappa) / sim$truth),
                      kappa_irt = irt$median[irt$variable == "kappa_pop"],
                      age_lo = min(obs$age), age_hi = max(obs$age))
  cat(sprintf("\n=== %s ===\n", slug))
  cat(sprintf("  observed   : median kappa %6.2f  [IQR %5.2f-%5.2f]  (%d/%d words fitted)\n",
              median(ko$kappa), quantile(ko$kappa,.25), quantile(ko$kappa,.75),
              nrow(ko), n_distinct(obs$word)))
  cat(sprintf("  simulated  : median kappa %6.2f   against a simulated truth of %.2f\n",
              median(ks$kappa), sim$truth))
  cat(sprintf("  -> pooled estimator recovers %.0f%% of the true kappa\n",
              100 * median(ks$kappa) / sim$truth))
}

## ---- the LM side, through the identical routine ---------------------------------
lad <- readr::read_csv(here("fits","llm","ladder_bestval_finer.csv"), show_col_types = FALSE) |>
  mutate(surprisal = as.numeric(surprisal))
nb <- n_distinct(lad$words)
lm_k <- lad |> group_by(seed, word) |> filter(n_distinct(words) == nb) |>
  group_modify(~{ p <- four_pl_sc(log10(.x$words), .x$surprisal)
                  tibble(sc = p["sc"], rng = p["rng"]) }) |> ungroup() |>
  filter(keep_ok(sc, rng, 1)) |> mutate(kappa = 0.434 / sc)
cat(sprintf("\n=== LMs (CHILDES ladder, same routine) ===\n"))
cat(sprintf("  per-word kappa: median %6.2f  [IQR %5.2f-%5.2f]  (%d fits)\n",
            median(lm_k$kappa), quantile(lm_k$kappa,.25), quantile(lm_k$kappa,.75), nrow(lm_k)))

cat("\n=== HEADLINE: one estimator, one unit, both systems ===\n")
for (slug in names(res)) cat(sprintf("  %-14s child per-word kappa %6.2f   (attenuated estimator; truth recovers at %.0f%%)\n",
    slug, median(res[[slug]]$observed$kappa),
    100 * median(res[[slug]]$simulated$kappa) / res[[slug]]$truth))
cat(sprintf("  %-14s LM per-word kappa    %6.2f\n", "CHILDES LMs", median(lm_k$kappa)))

## Does the independent estimator land where the IRT model says it should? This is the
## check that matters: two entirely different machineries on the same data.
cat("\n=== attenuation-corrected per-word kappa vs the IRT kappa ===\n")
cat(sprintf("  %-13s %8s %9s %10s %9s %7s\n","dataset","raw","recovery","corrected","IRT","ratio"))
for (slug in names(res)) with(res[[slug]],
  cat(sprintf("  %-13s %8.2f %8.0f%% %10.2f %9.2f %7.2f\n",
              slug, kappa_obs, 100*recovery, kappa_corrected, kappa_irt,
              kappa_corrected / kappa_irt)))

dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(child = res, lm = lm_k, rng_min_child = RNG_MIN_CHILD),
        file.path(CACHE, "si_perword_4pl.rds"))
cat("\nwrote", file.path(CACHE, "si_perword_4pl.rds"), "\n")
