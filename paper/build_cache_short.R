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
##   columns: jj, item, delta_j   (posterior-median item difficulty)
## Then per dataset, mirroring the retired build_cache.R section 5:
##   t_50   = a0 * exp((delta_j - log_H - mu_xi) / kappa_pop)
##   N_word = exp(mu_xi) * exp(log_H) * t_50 * prob
## with prob = per-item token frequency (per-language table: english_word_freq /
## norwegian_word_freq / a Japanese source TBD) and lexical_class from the CDI
## item metadata. Panels laid out per dataset, matching the Fig 1 fan facets.
PSI <- function(slug) file.path(SUMM, paste0(slug, SFX, "_m3_psi.csv"))
if (!any(vapply(names(DATASETS), function(s) file.exists(PSI(s)), logical(1)))) {
  cat("Fig 2: SKIPPED -- no per-item delta_j exports yet (<slug>_a3_m3_psi.csv).\n",
      "        fig-efficiency stays on the existing fig4_exposure.rds until they land.\n", sep = "")
} else {
  stop("Fig 2 inputs detected but the per-dataset efficiency build is not wired yet ",
       "-- implement here once the psi export format + per-language frequency are settled.")
}
