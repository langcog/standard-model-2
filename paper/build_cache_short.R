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
## For a typical child, the number of exposures a word needs before it is learned.
## Word difficulty / age-of-acquisition comes from the MODEL's delta_j (NOT from
## frequency -- the journal §34 correction: delta_j is validated against PRODUCTION
## via emp_prod, since high-frequency function words are late-learned). Frequency
## enters only as the per-word exposure COUNT on the y-axis:
##   t_50   = a0 * exp((delta_j - log_H - mu_xi) / kappa_pop)   # AoA from delta_j
##   N_word = exp(mu_xi + log_H) * t_50 * prob                  # cumulative exposures
##
## The lean bayes_long bundles dropped CDI item metadata, so we recover the word +
## lexical_class per item via wordbankr (REQUIRES NETWORK at cache-build time): item
## ids are `ul:<uni_lemma>` (unambiguous cross-linked items) else `id:<item_def>`,
## joined back to get_item_data() by uni_lemma / item_definition. Per-language token
## frequency: english_word_freq / norwegian_word_freq. Japanese is EXCLUDED -- no
## CHILDES frequency pull exists for it (no pull_japanese_freq.R); addable later.
suppressPackageStartupMessages(library(wordbankr))
FIG2_CFG <- list(
  thal      = list(wl = "English (American)", form = "WS", freq = "english"),
  smith     = list(wl = "English (American)", form = "WS", freq = "english"),
  marchman  = list(wl = "English (American)", form = "WS", freq = "english"),
  norwegian = list(wl = "Norwegian",          form = "WS", freq = "norwegian"))
FREQ <- list(english   = readRDS(here("fits", "english_word_freq.rds")),
             norwegian = readRDS(here("fits", "norwegian_word_freq.rds")))
## Exposure-COUNT anchor. The descriptive M3 carries no input rate, so cumulative
## exposures are anchored to the paper's own population input-rate table
## (input_rate_table.rds): each per-source estimate is tokens/month. We take the
## span across sources as the plausible range and its median as the central
## anchor; the figure uses the central value, the text reports LO-HI.
## [MCF: confirm range = spread of the 6 per-source monthly means (167k-561k),
##  central = their median.]
.rates <- sort(unique(readRDS(here("paper", "cache", "input_rate_table.rds"))$mean_mo))
.rates <- .rates[is.finite(.rates)]
R_LO <- min(.rates); R_HI <- max(.rates); R_MID <- median(.rates)
clean_word <- function(x) tolower(trimws(sub("[ ]*\\(.*\\)$", "", x)))
KEEP_CLASS <- c("nouns", "predicates", "function_words")   # drop "other" (heterogeneous)
PSI <- function(slug) file.path(SUMM, paste0(slug, SFX, "_m3_psi.csv"))

eff_one <- function(slug, label) {
  cfg <- FIG2_CFG[[slug]]; if (is.null(cfg)) { cat("  skip", slug, "(no freq table)\n"); return(NULL) }
  pf <- PSI(slug); sf <- file.path(SUMM, paste0(slug, SFX, "_m3.summary.rds"))
  bf <- file.path(BL, paste0("bundle_", slug, SFX, ".rds"))
  if (!all(file.exists(c(pf, sf, bf)))) { cat("  skip", slug, "(missing inputs)\n"); return(NULL) }
  psi <- read.csv(pf) |> mutate(kind = sub(":.*", "", item), key = sub("^[^:]+:", "", item))
  s <- as.data.frame(readRDS(sf)); g <- function(v) s$median[s$variable == v]
  sd <- readRDS(bf)$stan_data
  mu_xi <- g("mu_xi"); kappa <- g("kappa_pop"); log_H <- sd$log_H; a0 <- sd$a0

  it   <- wordbankr::get_item_data(language = cfg$wl, form = cfg$form) |> filter(item_kind == "word")
  byul <- it |> filter(!is.na(uni_lemma)) |> distinct(uni_lemma, item_definition, lexical_category)
  byid <- it |> distinct(item_definition, lexical_category)
  m <- bind_rows(
    psi |> filter(kind == "ul") |> left_join(byul, by = c("key" = "uni_lemma")),
    psi |> filter(kind == "id") |> left_join(byid, by = c("key" = "item_definition")) |>
      mutate(item_definition = key))
  fr <- FREQ[[cfg$freq]] |> transmute(w = tolower(w), prob)

  m <- m |>
    mutate(word = coalesce(item_definition, key), w = clean_word(word)) |>
    left_join(fr, by = "w") |>
    filter(!is.na(lexical_category), lexical_category %in% KEEP_CLASS, !is.na(prob), prob > 0) |>
    mutate(lexical_class = factor(lexical_category, levels = KEEP_CLASS),
           t_50   = a0 * exp((delta_j - log_H - mu_xi) / kappa),  # AoA from delta_j (model)
           # cumulative exposures = (input rate) x (months to learn) x (word share);
           # input rate anchored externally (R_MID), lo/hi span the source range.
           N_word    = R_MID * t_50 * prob,
           N_word_lo = R_LO  * t_50 * prob,
           N_word_hi = R_HI  * t_50 * prob,
           slug = slug, lang = label)
  r <- cor(m$delta_j, m$emp_prod)                        # §34 validation
  cat(sprintf("  %-10s items=%d  cor(delta_j,emp_prod)=%.2f\n", slug, nrow(m), r))
  if (is.na(r) || r > -0.5)
    stop(sprintf("Fig 2 %s: cor(delta_j, emp_prod)=%.2f -- expected strongly negative (§34)", slug, r))
  list(items = m |> select(slug, lang, item, word, lexical_class, delta_j, emp_prod,
                           prob, t_50, N_word, N_word_lo, N_word_hi),
       meta  = data.frame(slug = slug, lang = label, cor_dj_prod = r, n_items = nrow(m)))
}
cat("Fig 2: per-dataset exposures-to-learn (wordbankr metadata join)\n")
e2 <- Filter(Negate(is.null), Map(eff_one, names(DATASETS), DATASETS))
fig2 <- list(items = bind_rows(lapply(e2, `[[`, "items")) |> mutate(lang = factor(lang, levels = unname(DATASETS))),
             meta  = bind_rows(lapply(e2, `[[`, "meta")),
             anchor = list(lo = R_LO, mid = R_MID, hi = R_HI))  # tokens/month input-rate anchor
saveRDS(fig2, file.path(CACHE, "fig2_efficiency.rds"))
cat(sprintf("Wrote %s (%d datasets, Japanese excluded -- no frequency)\n\n",
            file.path(CACHE, "fig2_efficiency.rds"), nrow(fig2$meta)))

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

## ============ 5. SI: per-child BLUPs (efficiency xi, acceleration kappa) ==
## From the M3 per-child exports (<slug>_a3_m3_child.csv: xi_median, kappa_median).
## Replaces the retired glmer blups_demographics.rds for the dip-test / histogram
## in "Characterizing Variation".
cat("SI: per-child BLUPs (M3)\n")
child_one <- function(slug, label) {
  f <- file.path(SUMM, paste0(slug, SFX, "_m3_child.csv"))
  if (!file.exists(f)) { cat("  skip", slug, "(no child csv)\n"); return(NULL) }
  read.csv(f) |> transmute(lang = label, ckey, n_admins,
                           xi = xi_median, kappa = kappa_median)
}
si_blups <- bind_rows(Map(child_one, names(DATASETS), DATASETS)) |>
  mutate(lang = factor(lang, levels = unname(DATASETS)))
saveRDS(si_blups, file.path(CACHE, "si_blups.rds"))
cat(sprintf("Wrote %s (%d children, %d datasets)\n",
            file.path(CACHE, "si_blups.rds"), nrow(si_blups), dplyr::n_distinct(si_blups$lang)))
