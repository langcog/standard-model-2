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

## ============ 2. Fig 2: English exposures-to-learn (pooled) ==============
## A single English figure, so the item cloud can be shown and words labeled. We
## take a sample-weighted (by n_kids) average of the three English M3 fits
## (thal/smith/marchman): pooled population params (mu_xi, kappa_pop) and pooled
## per-item difficulty delta_j / production emp_prod, over items present in all
## three. Word AoA comes from the MODEL's delta_j, NOT frequency (§34: delta_j is
## validated against PRODUCTION via emp_prod, since high-frequency function words
## are late-learned). Frequency enters only as the per-word exposure COUNT:
##   t_50   = a0 * exp((delta_j - log_H - mu_xi) / kappa_pop)   # AoA from delta_j
##   N_word = R * t_50 * prob            # R = external input-rate anchor (below)
## word + lexical_class recovered via wordbankr (get_item_data, English WS;
## REQUIRES NETWORK at build time). Frequency: english_word_freq.
suppressPackageStartupMessages(library(wordbankr))
EN <- c("thal", "smith", "marchman")
## Exposure-count anchor: the descriptive M3 has no input rate, so anchor
## cumulative exposures to the paper's population input-rate table
## (input_rate_table.rds, tokens/month). Figure uses the central (median-of-
## sources) value; lo/hi span the range for in-text reporting.
## [MCF: confirm range = spread of the 6 per-source monthly means, central = median.]
.rates <- sort(unique(readRDS(here("paper", "cache", "input_rate_table.rds"))$mean_mo))
.rates <- .rates[is.finite(.rates)]
R_LO <- min(.rates); R_HI <- max(.rates); R_MID <- median(.rates)
clean_word <- function(x) tolower(trimws(sub("[ ]*\\(.*\\)$", "", x)))
KEEP_CLASS <- c("nouns", "predicates", "function_words")   # drop "other" (heterogeneous)
PSI <- function(slug) file.path(SUMM, paste0(slug, SFX, "_m3_psi.csv"))

cat("Fig 2: English exposures-to-learn (sample-weighted pool of thal/smith/marchman)\n")
parts <- lapply(EN, function(slug) {
  b <- readRDS(file.path(BL, paste0("bundle_", slug, SFX, ".rds")))
  s <- as.data.frame(readRDS(file.path(SUMM, paste0(slug, SFX, "_m3.summary.rds"))))
  g <- function(v) s$median[s$variable == v]
  list(w = b$meta$n_kids, mu_xi = g("mu_xi"), kappa = g("kappa_pop"),
       log_H = b$stan_data$log_H, a0 = b$stan_data$a0,
       psi = read.csv(PSI(slug)) |> transmute(item, delta_j, emp_prod))
})
W <- vapply(parts, `[[`, numeric(1), "w"); wsum <- sum(W)
wavg <- function(f) sum(W * vapply(parts, `[[`, numeric(1), f)) / wsum
mu_xi <- wavg("mu_xi"); kappa <- wavg("kappa"); log_H <- wavg("log_H"); a0 <- wavg("a0")
## sample-weighted per-item delta_j / emp_prod, over items present in all three
psi <- bind_rows(Map(function(p, w) mutate(p$psi, w = w), parts, W)) |>
  group_by(item) |> filter(n() == length(EN)) |>
  summarise(delta_j = weighted.mean(delta_j, w), emp_prod = weighted.mean(emp_prod, w),
            .groups = "drop") |>
  mutate(kind = sub(":.*", "", item), key = sub("^[^:]+:", "", item))
cat(sprintf("  pooled: n_kids=%.0f  mu_xi=%.2f kappa=%.2f  items=%d\n",
            wsum, mu_xi, kappa, nrow(psi)))

## wordbankr metadata (word + lexical_class) + CHILDES frequency
it   <- wordbankr::get_item_data(language = "English (American)", form = "WS") |> filter(item_kind == "word")
byul <- it |> filter(!is.na(uni_lemma)) |> distinct(uni_lemma, item_definition, lexical_category)
byid <- it |> distinct(item_definition, lexical_category)
m <- bind_rows(
  psi |> filter(kind == "ul") |> left_join(byul, by = c("key" = "uni_lemma")),
  psi |> filter(kind == "id") |> left_join(byid, by = c("key" = "item_definition")) |>
    mutate(item_definition = key))
fr <- readRDS(here("fits", "english_word_freq.rds")) |> transmute(w = tolower(w), prob)
items <- m |>
  mutate(word = coalesce(item_definition, key), w = clean_word(word)) |>
  left_join(fr, by = "w") |>
  filter(!is.na(lexical_category), lexical_category %in% KEEP_CLASS, !is.na(prob), prob > 0) |>
  mutate(lang = factor("English"), lexical_class = factor(lexical_category, levels = KEEP_CLASS),
         t_50   = a0 * exp((delta_j - log_H - mu_xi) / kappa),
         N_word    = R_MID * t_50 * prob,
         N_word_lo = R_LO  * t_50 * prob,
         N_word_hi = R_HI  * t_50 * prob) |>
  select(lang, item, word, lexical_class, delta_j, emp_prod, prob, t_50,
         N_word, N_word_lo, N_word_hi)
r <- cor(items$delta_j, items$emp_prod)                     # §34 validation
cat(sprintf("  items=%d  cor(delta_j,emp_prod)=%.2f\n", nrow(items), r))
if (is.na(r) || r > -0.5)
  stop(sprintf("Fig 2: cor(delta_j, emp_prod)=%.2f -- expected strongly negative (§34)", r))
fig2 <- list(items = items,
             meta  = data.frame(lang = "English", n_kids = wsum, mu_xi = mu_xi,
                                kappa = kappa, cor_dj_prod = r, n_items = nrow(items)),
             anchor = list(lo = R_LO, mid = R_MID, hi = R_HI))  # tokens/month input-rate anchor
saveRDS(fig2, file.path(CACHE, "fig2_efficiency.rds"))
cat(sprintf("Wrote %s (English pooled, %d items)\n\n",
            file.path(CACHE, "fig2_efficiency.rds"), nrow(items)))

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

## ============ 6. SI: datasets table (from the 5 M3 bundles) ==============
## Single-source replacement for the retired build_table1.R / table1_datasets.csv
## (which pooled glmer + io/proc + cross-sectional). Just the 5 longitudinal
## bundles now. Citations are plain text (kable can't render [@key] markers);
## Thal/Smith/Marchman years are placeholders. [MCF: confirm citations/years.]
cat("SI: datasets table\n")
DS_CITE <- c(thal = "Thal et al. (20XX)", smith = "Smith et al. (20XX)",
             marchman = "Marchman et al. (20XX)", norwegian = "Simonsen et al. (2014)",
             japanese = "Hagihara et al. (2023)")
DS_LANG <- c(thal = "English (American)", smith = "English (American)",
             marchman = "English (American)", norwegian = "Norwegian",
             japanese = "Japanese")
ds_one <- function(slug, label) {
  bf <- file.path(BL, paste0("bundle_", slug, SFX, ".rds"))
  if (!file.exists(bf)) return(NULL)
  b <- readRDS(bf); m <- b$meta; ag <- b$stan_data$admin_age
  data.frame(citation = DS_CITE[[slug]], language = DS_LANG[[slug]],
             n_kids = m$n_kids, n_admins = m$n_admins,
             min_age = min(ag), max_age = max(ag), mean_age = mean(ag),
             med_admins = m$med_admins_per_kid, stringsAsFactors = FALSE)
}
si_datasets <- bind_rows(Map(ds_one, names(DATASETS), DATASETS))
saveRDS(si_datasets, file.path(CACHE, "si_datasets.rds"))
cat(sprintf("Wrote %s (%d datasets, %s children total)\n",
            file.path(CACHE, "si_datasets.rds"), nrow(si_datasets),
            format(sum(si_datasets$n_kids), big.mark = ",")))

## ============ 7. Inline text values (rendered as `r ...` in the paper) ====
## One reproducible bundle of the numbers that appear inline in
## standard_model_short.qmd, so they are cache-derived rather than hand-typed.
## Raw values here; the qmd preamble formats them to the intended precision.
cat("Inline text values\n")
kap <- bind_rows(lapply(names(DATASETS), function(slug) {
  r <- as.data.frame(readRDS(file.path(SUMM, paste0(slug, SFX, "_m3.summary.rds"))))
  r <- r[r$variable == "kappa_pop", ]
  data.frame(slug = slug, med = r$median, q5 = r$q5, q95 = r$q95)
}))
klo <- kap[which.min(kap$med), ]; khi <- kap[which.max(kap$med), ]  # min / max kappa dataset
qc_pct <- function(slug) {                                          # % children QC-excluded
  m <- readRDS(file.path(BL, paste0("bundle_", slug, SFX, ".rds")))$meta
  100 * m$n_qc_dropped / (m$n_kids + m$n_qc_dropped)
}
llm <- readRDS(file.path(CACHE, "fig6_llm_slopes.rds"))             # children EN/NO kappa
kap_grp <- function(g) {
  v <- llm$slopes$slope_natural[llm$slopes$group == g]
  c(median = median(v), sd = sd(v))
}
en <- kap_grp("Children (English)"); no <- kap_grp("Children (Norwegian)")
inline <- list(
  age_lo = min(si_datasets$min_age), age_hi = max(si_datasets$max_age),
  loo_min = min(abs(si_loo$elpd_diff[si_loo$model == "M2"])),  # smallest M3-vs-next-best gap
  kappa_lo = klo$med, kappa_lo_q5 = klo$q5, kappa_lo_q95 = klo$q95,
  kappa_hi = khi$med, kappa_hi_q5 = khi$q5, kappa_hi_q95 = khi$q95,
  en_kappa = unname(en["median"]), no_kappa = unname(no["median"]),
  en_sd = unname(en["sd"]), no_sd = unname(no["sd"]),
  qc_marchman = qc_pct("marchman"), qc_norwegian = qc_pct("norwegian"))
saveRDS(inline, file.path(CACHE, "si_inline.rds"))
cat(sprintf("Wrote %s\n", file.path(CACHE, "si_inline.rds")))
