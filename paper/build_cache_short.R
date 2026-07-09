## build_cache_short.R -- lean cache for the SHORT paper (3 main figures).
##
## Reads the by-dataset Bayesian M3 accumulator fits from studies/bayes_long
## (fits/bayes_long/) and writes only what the short draft needs:
##   paper/cache/fig1_fan.rds        Fig 1: M0-M3 schematic + per-dataset M3 fan
##   paper/cache/fig2_efficiency.rds Fig 2: per-dataset per-word exposures-to-learn
##                                   (PENDING a per-item delta_j export -- see below)
## Fig 3 (LLM) reuses the existing paper/cache/fig6_llm_slopes.rds unchanged.
##
## Model numbering (new bayes_long ladder): M0 = kappa=1 pure accumulator (LLM
## analog); M1 = +acceleration; M2 = +per-child efficiency; M3 = +per-child
## acceleration (the headline). Fits are the "_a3" (3+-admin) variant.
##
## Run from the repo root:  Rscript paper/build_cache_short.R
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(readr); library(tibble)
})
CACHE <- here("paper", "cache")
BL    <- here("fits", "bayes_long")
SUMM  <- file.path(BL, "summaries")
SFX   <- "_a3"                              # the 3+-admin fit variant
DATASETS <- c(thal = "English (Thal)", smith = "English (Smith)",
              marchman = "English (Marchman)", norwegian = "Norwegian",
              japanese = "Japanese")
QS <- c(.1, .25, .5, .75, .9)
QLAB <- c("10th", "25th", "50th", "75th", "90th")

## ============ 1. Fig 1: schematic (M0-M3) + per-dataset M3 fan ============

## ---- (A) illustrative schematic: the four rungs at hard-coded params ----
## M0 pure accumulator (kappa=1, no between-child var) -> M3 full (per-child
## efficiency + acceleration). Purely illustrative (not a fit).
ages_s <- seq(8, 30, length.out = 80); A0s <- 8; Ls <- log(ages_s / A0s)
zq <- qnorm(QS); XI0 <- -2.0; KAPPA_A <- 2.5
DELTA_S <- qnorm(ppoints(150), mean = 0, sd = 1.5)
vocab_of <- function(th) vapply(th, function(x) mean(plogis(x - DELTA_S)), numeric(1))
variants <- tribble(
  ~name,                       ~kappa,  ~sigma_a, ~sigma_k,
  "(M0) accumulator",          1,       0,        0,
  "(M1) + acceleration",       KAPPA_A, 0,        0,
  "(M2) + efficiency var.",    KAPPA_A, 1.1,      0,
  "(M3) + per-child accel.",   KAPPA_A, 1.1,      0.7)
schematic <- variants |>
  mutate(name = factor(name, levels = name)) |>
  rowwise() |>
  do({ v <- .
    bind_rows(lapply(seq_along(QS), function(i) {
      th <- XI0 + zq[i] * v$sigma_a + (v$kappa + zq[i] * v$sigma_k) * Ls
      tibble(name = v$name, qf = factor(QLAB[i], levels = QLAB),
             age = ages_s, vocab = vocab_of(th)) })) }) |>
  ungroup()

## ---- (B) per-dataset model-implied fan from the M3 scalar posteriors ----
## Synthetic population (xi, kappa) ~ MVN + synthetic item difficulties; overlay
## 10/25/50/75/90 vocab quantiles on empirical trajectories. Mirrors 03_fan.R.
set.seed(1); N_SIM <- 500; N_SPAG <- 150; M_ITEM <- 500
fan_one <- function(slug, label) {
  sf <- file.path(SUMM, paste0(slug, SFX, "_m3.summary.rds"))
  bf <- file.path(BL,   paste0("bundle_", slug, SFX, ".rds"))
  if (!file.exists(sf) || !file.exists(bf)) { cat("  skip", slug, "(no m3 fit)\n"); return(NULL) }
  s  <- as.data.frame(readRDS(sf)); g <- function(v) s$median[s$variable == v]
  sd <- readRDS(bf)$stan_data
  mu_xi <- g("mu_xi"); delta <- g("delta"); sa <- g("sigma_a"); sb <- g("sigma_b")
  rho <- g("rho_ab"); tau <- g("tau_delta"); log_H <- sd$log_H; a0 <- sd$a0

  Sig <- matrix(c(sa^2, rho * sa * sb, rho * sa * sb, sb^2), 2)
  Z <- matrix(rnorm(N_SIM * 2), N_SIM, 2) %*% chol(Sig)
  xi <- mu_xi + Z[, 1]; kappa <- 1 + delta + Z[, 2]
  base_j <- log_H - rnorm(M_ITEM, 0, tau)     # synthetic item difficulties

  emp <- tibble(aa = sd$aa, y = sd$y) |>
    group_by(aa) |> summarise(prop = mean(y), .groups = "drop") |>
    mutate(age = sd$admin_age[aa], child = sd$admin_to_child[aa])
  aw <- quantile(emp$age, c(.05, .95)); ages <- seq(floor(aw[1]), ceiling(aw[2]), by = 0.5)
  kids <- sample(unique(emp$child), min(N_SPAG, n_distinct(emp$child)))
  spag <- emp |> filter(child %in% kids) |> transmute(lang = label, child, age, prop)

  fan <- lapply(ages, function(t) {
    A <- log(t / a0); v <- rowMeans(plogis(outer(xi + kappa * A, base_j, "+")))
    tibble(lang = label, age = t, q = QS, vocab = quantile(v, QS, names = FALSE)) }) |>
    bind_rows() |> mutate(qf = factor(q, levels = QS, labels = QLAB))

  cat(sprintf("  %-10s kappa=%.2f sigma_a=%.2f sigma_b=%.2f rho=%.2f\n",
              slug, 1 + delta, sa, sb, rho))
  list(fan = fan, spag = spag,
       meta = tibble(slug = slug, lang = label, kappa_pop = 1 + delta,
                     sigma_a = sa, sigma_b = sb, rho = rho,
                     n_kids = length(unique(emp$child))))
}
cat("Fig 1: per-dataset M3 fan\n")
res <- Filter(Negate(is.null), Map(fan_one, names(DATASETS), DATASETS))
lvl <- unname(DATASETS)
fig1 <- list(
  schematic = schematic,
  fan  = bind_rows(lapply(res, `[[`, "fan"))  |> mutate(lang = factor(lang, levels = lvl)),
  spag = bind_rows(lapply(res, `[[`, "spag")) |> mutate(lang = factor(lang, levels = lvl)),
  meta = bind_rows(lapply(res, `[[`, "meta")))
saveRDS(fig1, file.path(CACHE, "fig1_fan.rds"))
cat(sprintf("Wrote %s (%d datasets)\n\n", file.path(CACHE, "fig1_fan.rds"), nrow(fig1$meta)))

## ============ 2. Fig 2: per-dataset exposures-to-learn ===================
## PENDING a per-item difficulty export from studies/bayes_long/01_fit.R.
## Expected per dataset: fits/bayes_long/summaries/<slug>_a3_m3_psi.csv
##   columns: jj, item, delta_j, emp_prod
##     delta_j  = posterior-median item difficulty (THE source of word AoA)
##     emp_prod = empirical per-item production rate (for the §34 validation)
##
## A word's age-of-acquisition comes from the MODEL's delta_j, NOT from raw
## corpus frequency. This is the journal §34 correction: delta_j must be
## validated against PRODUCTION (emp_prod), not frequency -- high-frequency
## function words ("in", "if") are late-learned, so cor(delta_j, freq) ~ 0 and a
## frequency-based check false-alarms. Guard: cor(delta_j, emp_prod) strongly
## negative (easy words are produced by more children). Show this validation.
##
## Per dataset, mirroring the retired build_cache.R section 5:
##   t_50   = a0 * exp((delta_j - log_H - mu_xi) / kappa_pop)   # AoA from delta_j
##   N_word = exp(mu_xi) * exp(log_H) * t_50 * prob             # cumulative exposures
## Word difficulty/AoA = delta_j. The `prob` on the y-axis (exposures) is per-item
## token frequency purely as the exposure COUNT (not as a difficulty proxy) --
## CONFIRM with the fitting session whether the new figure keeps that y-axis or
## switches to a delta_j/emp_prod-only quantity. Per-language freq: english_word_freq
## / norwegian_word_freq / Japanese TBD; lexical_class from CDI item metadata.
## Panels laid out per dataset, matching the Fig 1 fan facets.
## The _psi.csv exports now exist (delta_j + emp_prod), BUT the exposures figure
## also needs each item's word form -> token frequency and its lexical class, and
## the lean bayes_long bundles dropped that CDI item metadata: item ids are
## prefixed id:/ul: uni-lemmas (~25/681 join to english_word_freq), and there is
## no Japanese frequency table. So the per-dataset exposures figure is still
## blocked on an item-metadata source (word + lexical_class per item). The
## delta_j-vs-emp_prod §34 validation, by contrast, needs only the psi file.
cat("Fig 2: SKIPPED -- psi exports present, but item metadata (word->frequency,\n",
    "        lexical_class) is not; the exposures figure needs a metadata source.\n",
    "        fig-efficiency stays on fig4_exposure.rds for now.\n", sep = "")

## ============ 3. SI: LOO model-comparison ladder (M0-M3) =================
## Per dataset: loo_compare across the four rungs -> elpd_diff vs best (M3) + se.
suppressPackageStartupMessages(library(loo))
LADDER <- c(M0 = "m0", M1 = "m1", M2 = "m2", M3 = "m3")
loo_one <- function(slug, label) {
  fs <- file.path(SUMM, paste0(slug, SFX, "_", LADDER, ".loo.rds"))
  if (!all(file.exists(fs))) { cat("  skip", slug, "(missing loo rungs)\n"); return(NULL) }
  ls <- lapply(fs, readRDS); names(ls) <- names(LADDER)
  cmp <- loo::loo_compare(ls)
  data.frame(slug = slug, lang = label, model = rownames(cmp),
             elpd_diff = cmp[, "elpd_diff"], se_diff = cmp[, "se_diff"],
             row.names = NULL)
}
cat("SI: LOO ladder (M0-M3)\n")
si_loo <- bind_rows(Map(loo_one, names(DATASETS), DATASETS)) |>
  mutate(lang = factor(lang, levels = unname(DATASETS)),
         model = factor(model, levels = names(LADDER)))
saveRDS(si_loo, file.path(CACHE, "si_loo.rds"))
cat(sprintf("Wrote %s (%d datasets)\n\n",
            file.path(CACHE, "si_loo.rds"), dplyr::n_distinct(si_loo$slug)))

## ============ 4. SI: log-age vs linear-age (M3 vs M3-linear) =============
## Per dataset: loo_compare(M3-log, M3-linear). elpd_diff is the loser's deficit
## relative to the winner (log wins everywhere so far; Norwegian m3lin pending).
loglin_one <- function(slug, label) {
  f3 <- file.path(SUMM, paste0(slug, SFX, "_m3.loo.rds"))
  fl <- file.path(SUMM, paste0(slug, SFX, "_m3lin.loo.rds"))
  if (!file.exists(f3) || !file.exists(fl)) { cat("  skip", slug, "(no m3lin)\n"); return(NULL) }
  cmp <- loo::loo_compare(list(log = readRDS(f3), linear = readRDS(fl)))
  loser <- rownames(cmp)[2]
  data.frame(slug = slug, lang = label, winner = rownames(cmp)[1],
             elpd_diff = cmp[loser, "elpd_diff"], se_diff = cmp[loser, "se_diff"],
             row.names = NULL)
}
cat("SI: log vs linear age (M3 vs M3-linear)\n")
si_loglin <- bind_rows(Map(loglin_one, names(DATASETS), DATASETS)) |>
  mutate(lang = factor(lang, levels = unname(DATASETS)))
saveRDS(si_loglin, file.path(CACHE, "si_loglin.rds"))
cat(sprintf("Wrote %s (%d datasets)\n",
            file.path(CACHE, "si_loglin.rds"), nrow(si_loglin)))
