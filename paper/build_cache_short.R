## build_cache_short.R -- lean cache for the SHORT paper (3 main figures).
##
## Reads the by-dataset Bayesian M3 accumulator fits from studies/bayes_long
## (fits/bayes_long/) and writes only what the short draft needs:
##   paper/cache/fig1_fan.rds        Fig 1: M0-M3 schematic + per-dataset M3 fan
##   paper/cache/fig2_efficiency.rds Fig 2: per-dataset per-word exposures-to-learn
##                                   (PENDING a per-item delta_j export -- see below)
## Fig 3A (LLM slope density) reuses the existing paper/cache/fig6_llm_slopes.rds
## (slopes/summary built by build_cache.R); section 8 below AUGMENTS that rds with
## the Fig 3B scaling-ladder data (scaling_bud/scaling_par). Run order: build_cache.R
## (builds slopes/summary) THEN build_cache_short.R (augments + fig1/fig2 + inline).
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
## MAIN TEXT = 3+ (_a3): the ladder, the log-vs-linear comparison, and the mega-model
## are all 3+, so the main-text figures/numbers use the same threshold. This matters most
## for sigma_b: at 2+ the two-wave children inflate it (Smith 5.1->8.1, Marchman 6.6->6.9,
## Norwegian 5.7->7.6) because a two-point slope can't separate acceleration from
## trajectory shape. The SI reports BOTH ladders (see section 3).
SFX   <- "_a3"
SFX2  <- ""                                 # the 2+ (all-data) variant, for the SI dual ladder
                                            # NB: sample/QC numbers live in the self-contained
                                            # paper/cache/bayes_long_sample.rds (built by
                                            # studies/bayes_long/qc_exclusion_report.R), which
                                            # only needs the bundles -- so the methods numbers
                                            # render before these fit-dependent caches rebuild.
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

## ---- (B) per-dataset model-implied fan from the FITTED per-child posteriors ----
## Quantiles over each child's fitted (xi_i, kappa_i) with the fitted item difficulties
## -- NOT an MVN resample of the population variance. The fitted kappa distribution is
## right-skewed (see supplement, fig-kappa-skew), so a Gaussian resample invents
## low-acceleration children and over-disperses the lower tail at older ages.
set.seed(1); N_SPAG <- 150
fan_one <- function(slug, label) {
  sf <- file.path(SUMM, paste0(slug, SFX, "_m3.summary.rds"))
  bf <- file.path(BL,   paste0("bundle_", slug, SFX, ".rds"))
  cf <- file.path(SUMM, paste0(slug, SFX, "_m3_child.csv"))
  pf <- file.path(SUMM, paste0(slug, SFX, "_m3_psi.csv"))
  if (!all(file.exists(c(sf, bf, cf, pf)))) { cat("  skip", slug, "(missing m3 fit/exports)\n"); return(NULL) }
  s  <- as.data.frame(readRDS(sf)); g <- function(v) s$median[s$variable == v]
  sd <- readRDS(bf)$stan_data; log_H <- sd$log_H; a0 <- sd$a0
  sa <- g("sigma_a"); sb <- g("sigma_b"); rho <- g("rho_ab"); delta <- g("delta")
  ch <- read.csv(cf); dj <- read.csv(pf)$delta_j      # fitted per-child (xi,kappa) + item difficulties
  xi <- ch$xi_median; kappa <- ch$kappa_median; base_j <- log_H - dj

  emp <- tibble(aa = sd$aa, y = sd$y) |>
    group_by(aa) |> summarise(prop = mean(y), .groups = "drop") |>
    mutate(age = sd$admin_age[aa], child = sd$admin_to_child[aa])
  ages <- seq(floor(min(emp$age)), ceiling(max(emp$age)), by = 0.5)   # full observed range (bands cover the data)
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

## ---- (A2) conceptual theta -> CDI mechanism (short-paper Fig 1, block B) ----
## Replaces the old M0-M3 schematic in the short paper: shows the pure vs
## accelerating accumulator and the two dimensions of individual variation
## (efficiency xi = a level shift; acceleration kappa = a fan) in latent-ability
## (theta) space and projected into words-produced (CDI) space. Illustrative.
## Tuned so the CDI (right) facet shows a FULL sigmoid inside the plotted window:
## the Baseline accelerator sits at ~0 until it takes off at ~14 months and saturates
## (~0.94) by 36, which makes it visually distinct from the concave log curve in the
## theta facet. Narrowing the difficulty spread (sd 1.5 -> 0.8) sharpens the sweep
## through the item distribution.
## Plotted range starts at 8 months to match the DATA range in block C (Marchman and
## Norwegian both begin at 8), so panels B and C share an x-origin. Note a0c stays 9:
## it is only the internal anchor where log(t/a0c)=0, and moving it would shift the
## whole curve up by kappa*log(9/8) ~ 0.77 and pull the takeoff earlier, undoing the
## tuning. Extending further left is not useful -- the curve is already at 0.002 by 8mo
## (visually zero, so the bottom of the S is fully shown), and t=0 is a log singularity
## (theta -> -Inf) so it cannot be plotted at all.
a0c <- 9; ages_c <- seq(8, 36, length.out = 200); Lc <- log(ages_c / a0c)
XI_C <- -6; KA_C <- 6.5; SD_D <- 0.8
DELTA_C <- qnorm(ppoints(400), 0, SD_D)
vocab_c <- function(th) vapply(th, function(x) mean(plogis(x - DELTA_C)), numeric(1))
qlab_c <- c(theta = "Latent Ability (θ)", cdi = "Words Produced (CDI)")
## Pure (kappa=1) gets its own REFIT intercept rather than sharing the accelerator's
## xi. Sharing xi buries it at ~0.01 (invisible); refitting by least squares against
## the Baseline curve is the M0-vs-M3 ladder comparison in miniature -- the best a
## pure accumulator can do -- and it makes the point harder: even at its optimum the
## unit accumulator starts too high (~0.30) and ends too low (~0.59), never tracing
## the S. Illustrative, not estimated.
base_cdi_c <- vocab_c(XI_C + KA_C * Lc)
XI_PURE <- optimize(function(xi) sum((vocab_c(xi + Lc) - base_cdi_c)^2), c(-8, 6))$minimum
conceptual <- tribble(
  ~scenario,        ~xi,          ~kappa,      ~kind,
  "Pure (κ=1)",     XI_PURE,      1.0,         "pure",
  "Baseline",       XI_C,         KA_C,        "accel",
  "↑ Efficiency",   XI_C + 1.6,   KA_C,        "accel",
  "↑ Acceleration", XI_C,         KA_C + 1.8,  "accel") |>
  mutate(scenario = factor(scenario, levels = scenario)) |>
  rowwise() |>
  mutate(d = list(tibble(age = ages_c, theta = xi + kappa * Lc, cdi = vocab_c(xi + kappa * Lc)))) |>
  ungroup() |> unnest(d) |>
  pivot_longer(c(theta, cdi), names_to = "q", values_to = "value") |>
  mutate(quantity = factor(qlab_c[q], levels = qlab_c))
conceptual_lab <- conceptual |> group_by(scenario, quantity) |> slice_max(age, n = 1) |> ungroup()

fig1 <- list(
  schematic  = schematic,          # retained (unused by short-paper Fig 1; kept for compatibility)
  conceptual = conceptual, conceptual_lab = conceptual_lab,
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
CLASS4 <- c("nouns", "action words", "descriptive words", "function words")  # predicates split; "other" dropped
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

## wordbankr metadata (word + lexical_class) + CHILDES frequency.
## The WS item dictionary is a static instrument definition, but fetching it needs a
## live Wordbank connection -- and an outage used to abort the whole script, blocking
## every LATER section (LOO ladder, BLUPs, datasets, random effects) that needs no
## network at all. So: fetch when reachable and mirror to paper/cache; fall back to
## the mirror when not. Delete the mirror to force a refresh.
WB_MIRROR <- file.path(CACHE, "wordbank_items_en_ws.rds")
## NB: on a failed connection wordbankr::get_item_data() prints a message and returns
## NULL rather than throwing, so tryCatch alone is not enough -- test the value.
it <- tryCatch({
  x <- wordbankr::get_item_data(language = "English (American)", form = "WS")
  if (is.null(x) || !nrow(x)) stop("wordbank returned no rows")
  saveRDS(x, WB_MIRROR); cat("  wordbank: live fetch (mirrored to cache)\n"); x
}, error = function(e) {
  if (file.exists(WB_MIRROR)) {
    cat("  wordbank: UNREACHABLE -- using local mirror", basename(WB_MIRROR), "\n")
    return(readRDS(WB_MIRROR))
  }
  cat("  wordbank: UNREACHABLE and no local mirror -- SKIPPING Fig 2,\n",
      "    keeping the committed fig2_efficiency.rds. Rerun when Wordbank is up.\n", sep = "")
  NULL
})
if (is.null(it)) fig2_skipped <- TRUE else {
fig2_skipped <- FALSE
it <- it |> filter(item_kind == "word")
byul <- it |> filter(!is.na(uni_lemma)) |> distinct(uni_lemma, item_definition, lexical_category, category)
byid <- it |> distinct(item_definition, lexical_category, category)
m <- bind_rows(
  psi |> filter(kind == "ul") |> left_join(byul, by = c("key" = "uni_lemma")),
  psi |> filter(kind == "id") |> left_join(byid, by = c("key" = "item_definition")) |>
    mutate(item_definition = key))
fr <- readRDS(here("fits", "english_word_freq.rds")) |> transmute(w = tolower(w), prob)
items <- m |>
  mutate(word = coalesce(item_definition, key), w = clean_word(word)) |>
  left_join(fr, by = "w") |>
  mutate(lex4 = dplyr::case_when(                       # split predicates; drop "other"
           lexical_category == "nouns"          ~ "nouns",
           category == "action_words"           ~ "action words",
           category == "descriptive_words"      ~ "descriptive words",
           lexical_category == "function_words" ~ "function words",
           TRUE ~ NA_character_)) |>
  filter(!is.na(lex4), !grepl(" ", w), !is.na(prob), prob > 0) |>  # drop multi-word items (bad unigram freq)
  mutate(lang = factor("English"), lexical_class = factor(lex4, levels = CLASS4),
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
## Plotted age range. Stored rather than hard-coded in both places, because the M0
## overlay lines below are evaluated over it and would silently disagree with the
## figure if the two drifted.
FIG2_XLIM <- c(10, 42)
fig2 <- list(items = items,
             meta  = data.frame(lang = "English", n_kids = wsum, mu_xi = mu_xi,
                                kappa = kappa, cor_dj_prod = r, n_items = nrow(items)),
             anchor = list(lo = R_LO, mid = R_MID, hi = R_HI),  # tokens/month input-rate anchor
             xlim = FIG2_XLIM)

## ---- M0 overlay for Fig 2 ------------------------------------------------
## Fig 2 draws M3's per-class fits; the accompanying text claims the fitted PURE
## accumulator predicts a far shallower decline. These are the M0 fits that back that
## sentence, so a reader can see the contrast rather than take it on faith.
##
## The coefficients come from the SAME M0 counterfactual the SI reports
## (studies/bayes_long/04_fig2_m0_counterfactual.R -> si_fig2_m0.rds), so the main-text
## overlay and SI: Age of Acquisition in the Pure Accumulator Model cannot drift apart.
##
## GEOMETRY. Fig 2 has LINEAR age on x and log10 exposures on y, so the line to draw is
## the semi-log fit lm(log10(N) ~ t_50) -- NOT the log-log elasticity the SI table
## reports, which belongs to that section's log-x figure. Fitted over all of M0's points
## (its AoAs run to 254 months) and evaluated over each class's plotted span, so each
## dashed line covers the same ages as the solid one it is being compared against.
m0f <- file.path(CACHE, "si_fig2_m0.rds")
fig2$m0_lines <- if (!file.exists(m0f)) {
  cat("  note: si_fig2_m0.rds absent -- Fig 2 M0 overlay skipped\n"); NULL
} else {
  span <- items |> filter(t_50 >= FIG2_XLIM[1], t_50 <= FIG2_XLIM[2]) |>
    group_by(lexical_class) |> summarise(x1 = min(t_50), x2 = max(t_50), .groups = "drop")
  readRDS(m0f)$m0 |> group_by(lexical_class) |>
    group_modify(~ tibble::tibble(b0 = coef(lm(log10(N_word) ~ t_50, data = .x))[1],
                                  b1 = coef(lm(log10(N_word) ~ t_50, data = .x))[2])) |>
    ungroup() |> left_join(span, by = "lexical_class") |> rowwise() |>
    reframe(lexical_class = lexical_class, t_50 = c(x1, x2),
            N_word = 10^(b0 + b1 * c(x1, x2)), slope_per_month = b1)
}
if (!is.null(fig2$m0_lines)) {
  cat("  M0 overlay slopes (log10 exposures per month), M3 for comparison:\n")
  m3s <- items |> group_by(lexical_class) |>
    summarise(m3 = coef(lm(log10(N_word) ~ t_50))[2], .groups = "drop")
  print(fig2$m0_lines |> distinct(lexical_class, slope_per_month) |>
        left_join(m3s, by = "lexical_class") |>
        transmute(class = lexical_class, M3 = round(m3, 4),
                  M0 = round(slope_per_month, 4)) |> as.data.frame(), row.names = FALSE)
}
saveRDS(fig2, file.path(CACHE, "fig2_efficiency.rds"))
cat(sprintf("Wrote %s (English pooled, %d items)\n\n",
            file.path(CACHE, "fig2_efficiency.rds"), nrow(items)))
}  # end if (!fig2_skipped)

## ============ 3. SI: LOO model-comparison ladder (M0-M3) =================
## Per dataset: loo_compare across the four rungs -> elpd_diff vs best (M3) + se.
suppressPackageStartupMessages(library(loo))
## Computed at BOTH thresholds: 3+ (the main-text analysis, where per-child slopes are
## identifiable) and 2+ (the full longitudinal sample). The SI reports both; the main
## text quotes the 3+ rows. `threshold` distinguishes them.
LADDER <- c(M0 = "m0", M1 = "m1", M2 = "m2", M3 = "m3")
loo_one <- function(slug, label, sfx, thr) {
  fs <- file.path(SUMM, paste0(slug, sfx, "_", LADDER, ".loo.rds"))
  if (!all(file.exists(fs))) { cat("  skip", slug, " (", thr, ": missing loo rungs)\n", sep=""); return(NULL) }
  ls <- lapply(fs, readRDS); names(ls) <- names(LADDER)
  cmp <- loo::loo_compare(ls)
  data.frame(slug = slug, lang = label, threshold = thr, model = rownames(cmp),
             elpd_diff = cmp[, "elpd_diff"], se_diff = cmp[, "se_diff"],
             row.names = NULL)
}
cat("SI: LOO ladder (M0-M3), both thresholds\n")
## NB: unname() the Map results -- bind_rows() with two *named* lists binds them as
## columns (thal/smith/...) instead of stacking rows.
si_loo <- bind_rows(
    unname(Map(function(s, l) loo_one(s, l, SFX,  "3+"), names(DATASETS), DATASETS)),
    unname(Map(function(s, l) loo_one(s, l, SFX2, "2+"), names(DATASETS), DATASETS))) |>
  mutate(lang = factor(lang, levels = unname(DATASETS)),
         model = factor(model, levels = names(LADDER)),
         threshold = factor(threshold, levels = c("3+", "2+")))
saveRDS(si_loo, file.path(CACHE, "si_loo.rds"))
cat(sprintf("Wrote %s (%d datasets)\n\n",
            file.path(CACHE, "si_loo.rds"), dplyr::n_distinct(si_loo$slug)))

## ============ 4. SI: log-age vs linear-age (M3 vs M3-linear) =============
## Per dataset: loo_compare(M3-log, M3-linear). elpd_diff is the loser's deficit
## relative to the winner (log wins everywhere so far; Norwegian m3lin pending).
## Log-vs-linear is a 3+ (_a3) analysis (m3lin was only fit there; the two-point kids
## can't identify per-child slope curvature). So force the 3+ suffix here regardless
## of the global SFX (which is "" / 2+ for the main ladder + Fig 1-2).
loglin_one <- function(slug, label) {
  f3 <- file.path(SUMM, paste0(slug, "_a3_m3.loo.rds"))
  fl <- file.path(SUMM, paste0(slug, "_a3_m3lin.loo.rds"))
  if (!file.exists(f3) || !file.exists(fl)) { cat("  skip", slug, "(no _a3 m3lin)\n"); return(NULL) }
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
  ## Pre-exclusion N: the full longitudinal sample (>=2 administrations, SFX2="").
  ## The main-text analysis is the >=3-admin subset (n_kids); reporting both makes
  ## the exclusion read as a modeling requirement (a two-admin child gives one slope
  ## with zero residual df, so kappa_i is unidentified) rather than a small dataset.
  bf2 <- file.path(BL, paste0("bundle_", slug, SFX2, ".rds"))
  n_all <- if (file.exists(bf2)) readRDS(bf2)$stan_data$I else NA_integer_
  data.frame(citation = DS_CITE[[slug]], language = DS_LANG[[slug]],
             n_kids = m$n_kids, n_kids_all = n_all, n_admins = m$n_admins,
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
## (QC exclusion % moved to paper/cache/bayes_long_sample.rds -- built from bundles only,
##  so it doesn't depend on these fits.)
llm <- readRDS(file.path(CACHE, "fig6_llm_slopes.rds"))             # children EN/NO kappa
kap_grp <- function(g) {
  v <- llm$slopes$slope_natural[llm$slopes$group == g]
  c(median = median(v), sd = sd(v))
}
en <- kap_grp("Children (English)"); no <- kap_grp("Children (Norwegian)")
## ---- random-effect parameters, for main-text reporting ------------------
## Reviewers reasonably ask for the between-child SD of kappa (sigma_b), the
## intercept-slope correlation (rho), their uncertainty, and the share of children
## above the kappa = 1 null -- these are more informative than distributional-shape
## claims. Assembled per dataset (si_ranef.rds) plus pooled scalars (inline).
ranef <- bind_rows(lapply(names(DATASETS), function(slug) {
  s <- as.data.frame(readRDS(file.path(SUMM, paste0(slug, SFX, "_m3.summary.rds"))))
  g <- function(v, q) s[[q]][s$variable == v]
  ch <- read.csv(file.path(SUMM, paste0(slug, SFX, "_m3_child.csv")))
  data.frame(slug = slug, lang = DATASETS[[slug]], n_kids = nrow(ch),
             kappa = g("kappa_pop","median"), kappa_q5 = g("kappa_pop","q5"), kappa_q95 = g("kappa_pop","q95"),
             sigma_a = g("sigma_a","median"), sigma_a_q5 = g("sigma_a","q5"), sigma_a_q95 = g("sigma_a","q95"),
             sigma_b = g("sigma_b","median"), sigma_b_q5 = g("sigma_b","q5"), sigma_b_q95 = g("sigma_b","q95"),
             rho = g("rho_ab","median"), rho_q5 = g("rho_ab","q5"), rho_q95 = g("rho_ab","q95"),
             pct_gt1 = 100 * mean(ch$kappa_median > 1),
             n_neg = sum(ch$kappa_median < 0), row.names = NULL)
})) |> mutate(lang = factor(lang, levels = unname(DATASETS)))
saveRDS(ranef, file.path(CACHE, "si_ranef.rds"))
cat(sprintf("  random effects: sigma_b %.2f-%.2f, rho %+.2f to %+.2f, %.1f%% of children kappa>1\n",
            min(ranef$sigma_b), max(ranef$sigma_b), min(ranef$rho), max(ranef$rho),
            100 * sum(ranef$pct_gt1/100 * ranef$n_kids) / sum(ranef$n_kids)))

inline <- list(
  age_lo = min(si_datasets$min_age), age_hi = max(si_datasets$max_age),
  # smallest M3-vs-next-best gap, from the MAIN-TEXT (3+) ladder only
  loo_min = min(abs(si_loo$elpd_diff[si_loo$model == "M2" & si_loo$threshold == "3+"])),
  kappa_lo = klo$med, kappa_lo_q5 = klo$q5, kappa_lo_q95 = klo$q95,
  kappa_hi = khi$med, kappa_hi_q5 = khi$q5, kappa_hi_q95 = khi$q95,
  en_kappa = unname(en["median"]), no_kappa = unname(no["median"]),
  en_sd = unname(en["sd"]), no_sd = unname(no["sd"]),
  ## random-effect summaries (pooled across the five 3+ fits)
  sb_lo = min(ranef$sigma_b), sb_hi = max(ranef$sigma_b),
  rho_lo = min(ranef$rho), rho_hi = max(ranef$rho),
  pct_kappa_gt1 = 100 * sum(ranef$pct_gt1/100 * ranef$n_kids) / sum(ranef$n_kids),
  n_kappa_neg = sum(ranef$n_neg), n_kids_tot = sum(ranef$n_kids),
  ## between- vs within-dataset uncertainty in kappa: the honest scale of uncertainty
  ## is the spread ACROSS samples, not the (model-conditional) within-sample interval,
  ## because CDI items are not locally independent.
  kappa_between_sd = sd(ranef$kappa),
  kappa_ci_width_min = min(ranef$kappa_q95 - ranef$kappa_q5))
saveRDS(inline, file.path(CACHE, "si_inline.rds"))
cat(sprintf("Wrote %s\n", file.path(CACHE, "si_inline.rds")))

## ============ 8. Fig 3 (LM panels B & C): scaling law + matched-input return ====
## Augment fig6_llm_slopes.rds (panel A = EN/NO child + LM slope densities, from
## build_cache.R) with panels B and C. On the distinct-input ladder, aggregate loss
## (mean CDI-word surprisal per data budget) follows the Chinchilla form
## L = E + B*D^-beta.
##   Panel B: excess loss (L - E) falls as a power law (slope -beta), shown with the
##            per-seed points so the model-to-model variability is visible.
##   Panel C: children and LMs on ONE dimensionless axis -- the fractional return
##            gamma = (reduction of excess loss per e-fold of experience) / (excess
##            loss remaining) -- at matched ABSOLUTE input (cumulative word tokens).
##            For a power law gamma = beta (LMs, constant). For the child accumulator
##            gamma_i(t) = kappa_i * mean_j(1 - p_ij) / mean_j(-log p_ij), whose
##            ceiling is the child's kappa_i. Child AGE is mapped to cumulative tokens
##            via the pooled input rate from studies/input_estimation (age = experience).
fig6        <- readRDS(file.path(CACHE, "fig6_llm_slopes.rds"))
ladder      <- read_csv(here("fits", "llm", "ladder_bestval_finer.csv"), show_col_types = FALSE)
scaling_bud <- ladder |>
  group_by(words) |> summarise(L = mean(surprisal), .groups = "drop") |> arrange(words)
scaling_fit <- nls(L ~ E + B * words^(-beta), data = scaling_bud,
                   start = list(E = 3, B = 50, beta = 0.3),
                   control = nls.control(maxiter = 500, warnOnly = TRUE))
scaling_par <- as.list(coef(scaling_fit))          # E (entropy floor), B, beta
scaling_bud <- scaling_bud |>
  mutate(kappa = scaling_par$beta * (L - scaling_par$E))   # retained for SI derivation
fig6$scaling_bud <- scaling_bud
fig6$scaling_par <- scaling_par
E <- scaling_par$E

## (B) per-seed excess loss (the scaling-law points, with variability)
seed_excess <- ladder |>
  group_by(seed, words) |> summarise(L = mean(surprisal), .groups = "drop") |>
  mutate(excess = L - E) |> filter(excess > 0)
## (C, LM) per-seed smooth local exponent gamma(D) = -d log(excess)/d log D
Dgrid <- 10^seq(log10(3e6), log10(2.4e7), length = 50)
lm_gamma <- seed_excess |> group_by(seed) |> group_modify(~{
  f <- lm(log(excess) ~ poly(log(words), 2, raw = TRUE), data = .x)
  b <- coef(f); ld <- log(Dgrid)
  tibble(words = Dgrid, gamma = -(b[2] + 2 * b[3] * ld)) }) |> ungroup() |> filter(gamma > 0)

## input rate + between-child SD of log rate, from the input-estimation validation set
## (POOLED row); single-sourced here rather than duplicated in the figure chunk.
vpool <- read_csv(here("studies", "input_estimation", "validation_set.csv"), show_col_types = FALSE) |>
  filter(grepl("POOLED", source)) |> slice(1)
rate_hr <- vpool$tokens_per_hour_mean; sig_r <- vpool$log_r_sd
tok_per_mo <- rate_hr * 365; f1 <- exp(sig_r)      # x / division per +/-1 SD of log input rate
a0 <- 18; log_H <- log(365)                        # model constants (Methods)

## (C, children) per-dataset item difficulties (delta_j) and measured age spans
## (spans per si_datasets.rds: Thal 12-29, Smith 16-30, Marchman 8-30).
EN_INFO <- tibble(lang = c("English (Thal)", "English (Smith)", "English (Marchman)"),
                  slug = c("thal", "smith", "marchman"), amin = c(12, 16, 8), amax = c(29, 30, 30))
dj_of <- lapply(setNames(EN_INFO$slug, EN_INFO$lang), function(s) read.csv(PSI(s))$delta_j)
gamma_curve <- function(xi, kap, dj, ages)
  vapply(ages, function(t) { p <- plogis(xi + kap * log(t / a0) + log_H - dj)
                             kap * mean(1 - p) / mean(-log(p)) }, numeric(1))
elig <- readRDS(file.path(CACHE, "si_blups.rds")) |>
  filter(lang %in% EN_INFO$lang, kappa > 0, n_admins >= 3) |> left_join(EN_INFO, by = "lang")
## all eligible English children: gamma over each child's DATASET measured age span
child_bg <- elig |> rowwise() |>
  reframe(ckey = ckey, age = seq(amin, amax, by = 0.4),
          gamma = gamma_curve(xi, kappa, dj_of[[lang]], seq(amin, amax, by = 0.4))) |>
  ungroup() |> mutate(tokens = age * tok_per_mo)
## three exemplars (low/med/high kappa) from Marchman (widest span), with +/-1 SD input ribbon
mar  <- filter(elig, slug == "marchman")
exk  <- quantile(mar$kappa, c(.1, .5, .9))
exid <- vapply(exk, function(k) mar$ckey[which.min(abs(mar$kappa - k))], character(1))
child_ex <- elig |> filter(ckey %in% exid) |> rowwise() |>
  reframe(ckey = ckey, kappa = kappa, age = seq(amin, amax, by = 0.25),
          gamma = gamma_curve(xi, kappa, dj_of[[lang]], seq(amin, amax, by = 0.25))) |>
  ungroup() |> mutate(tokens = age * tok_per_mo, lo = tokens / f1, hi = tokens * f1)

fig6$fig3 <- list(seed_excess = seed_excess, lm_gamma = lm_gamma,
                  child_bg = child_bg, child_ex = child_ex,
                  consts = list(rate_hr = rate_hr, sig_r = sig_r, tok_per_mo = tok_per_mo, f1 = f1,
                                band_lo = 8 * tok_per_mo, band_hi = 30 * tok_per_mo))

## (C) Panel C: between-learner variability in developmental kappa.
## Per-learner acceleration on one axis: children = per-child kappa_i BLUPs (the
## same _a3 M3 exports as panel A); LMs = per-LADDER median per-word 4-PL kappa
## on the distinct-input axis -- CHILDES (10 seeds x 18 rungs,
## ladder_bestval_finer.csv) and the composition control (register_bestval.csv;
## experiments log L7): BabyLM-minus-CHILDES and ClimbMix, each 3 disjoint
## subsets x 8 seeds x 12 rungs. ~35k nls fits, adds a few minutes to the build.
##
## CONVERGENCE WARNINGS. Many of these fits fail to converge and say so, which used to
## bury the build log under thousands of "singular convergence (7)" messages. The
## individual warnings are suppressed here, but NOT silently: `warnOnly = TRUE` already
## meant a failed fit returned its current parameter values rather than erroring, and the
## sc/rng filter below discards the degenerate ones, so the number that actually matters
## is how many words survive to contribute to a ladder's median. ladder_kappa() reports
## that retention rate per population instead.
##
## The rate is not incidental. CHILDES retains ~95% of words, but the register-control
## ladders retain only ~60%: about a third of their fits have sc pinned at the lower
## bound, meaning the word's surprisal drops between two adjacent rungs and a 4-PL scale
## is not identified over 12 budgets. Those ladder medians are therefore taken over a
## thinner set of words than the CHILDES ones, which is worth seeing in the log.
four_pl_sc <- function(x, y) tryCatch(suppressWarnings({
  f <- nls(y ~ lo + (up - lo) / (1 + exp((x - mid) / sc)),
           start = list(up = max(y), lo = min(y), mid = mean(x), sc = 0.5),
           lower = c(up = min(y)-5, lo = min(y)-5, mid = min(x)-3, sc = 1e-3),
           upper = c(up = max(y)+5, lo = max(y)+5, mid = max(x)+3, sc = 50),
           algorithm = "port", control = nls.control(maxiter = 200, warnOnly = TRUE))
  c(sc = unname(coef(f)["sc"]), rng = unname(coef(f)["up"] - coef(f)["lo"]))
}), error = function(e) c(sc = NA_real_, rng = NA_real_))
ladder_kappa <- function(d, idcols, label = "") {  # per-word 4-PL over budgets -> per-ladder median
  # one row per (ladder, word, budget): collapses duplicate evals logged at the
  # same final step (on_train_end + last scheduled step can coincide)
  d <- d |> distinct(across(all_of(c(idcols, "word", "words"))), .keep_all = TRUE)
  nb <- n_distinct(d$words)                       # full-ladder budget count
  raw <- d |> group_by(across(all_of(c(idcols, "word")))) |> filter(n_distinct(words) == nb) |>
    group_modify(~{ p <- four_pl_sc(log10(.x$words), .x$surprisal)
                    tibble(sc = p["sc"], rng = p["rng"]) }) |> ungroup()
  kept <- raw |> filter(is.finite(sc), sc > 0.01, sc < 10, rng > 1)
  ## Stands in for the thousands of suppressed nls warnings: how many words actually
  ## reach each ladder's median, and why the rest did not. Reasons are assigned in the
  ## order the filter applies, so they are disjoint and sum to the number dropped -- a
  ## word with both an unidentified scale and a flat curve is counted once, as the former.
  RSN <- c("unidentified scale", "flat curve", "error")   # factor levels, so a zero-count
  why <- factor(with(raw,                                 # reason still prints as 0
           ifelse(!is.finite(sc), "error",
           ifelse(sc <= 0.01 | sc >= 10, "unidentified scale",
           ifelse(rng <= 1, "flat curve", "kept")))), levels = c(RSN, "kept"))
  cat(sprintf("  %-22s %d/%d words fitted (%.1f%%)  [dropped: %s]\n",
              label, nrow(kept), nrow(raw), 100 * nrow(kept) / max(nrow(raw), 1),
              paste(sprintf("%d %s", table(why)[RSN], RSN), collapse = ", ")))
  kept |> mutate(kappa = 0.434 / sc) |>
    group_by(across(all_of(idcols))) |> summarise(kappa = median(kappa), .groups = "drop")
}
cat("\n  per-word 4-PL fits for Fig 3C (nls warnings suppressed; retention reported):\n")
reg <- read_csv(here("fits", "llm", "register_bestval.csv"), show_col_types = FALSE)
k_chi <- ladder_kappa(ladder |> mutate(surprisal = as.numeric(surprisal)), "seed",
                      "CHILDES") |>
  transmute(population = "LMs: CHILDES", id = paste0("s", seed), kappa)
k_reg <- ladder_kappa(reg, c("corpus", "subset", "seed"), "register control") |>
  transmute(population = ifelse(corpus == "babylm", "LMs: BabyLM", "LMs: ClimbMix"),
            id = paste(subset, seed), kappa)
k_kid <- fig6$slopes |>
  filter(group %in% c("Children (English)", "Children (Norwegian)")) |>
  transmute(population = as.character(group), id = paste0("c", row_number()),
            kappa = slope_natural)
fig6$vary <- bind_rows(k_kid, k_chi, k_reg) |>
  mutate(population = factor(population, levels = c(
    "Children (English)", "Children (Norwegian)",
    "LMs: CHILDES", "LMs: BabyLM", "LMs: ClimbMix")))

## ---- bootstrap CIs on the coefficient of variation (Fig 3C) --------------
## A reviewer noted that a CV computed from 8-10 seeds -- with BOTH the mean and the
## variance estimated from fitted kappa values -- can be unstable, and asked for
## intervals. We resample learners within population with replacement. The child
## populations have hundreds of learners so their CIs are tight; the LM populations have
## 8-30 and theirs are correspondingly wide, which is exactly the point to display
## honestly rather than to hide behind a point estimate.
set.seed(20260802); NBOOT <- 4000
fig6$vary_cv <- fig6$vary |>
  group_by(population) |>
  group_modify(~ {
    k <- .x$kappa; n <- length(k)
    cv <- function(v) 100 * sd(v) / mean(v)
    bs <- replicate(NBOOT, cv(sample(k, n, replace = TRUE)))
    bs <- bs[is.finite(bs)]
    tibble(n_learners = n, cv = cv(k),
           cv_lo = unname(quantile(bs, 0.025)), cv_hi = unname(quantile(bs, 0.975)))
  }) |> ungroup()
cat("Fig 3C: bootstrap CIs on the CV of kappa (2.5-97.5%)\n")
for (i in seq_len(nrow(fig6$vary_cv))) with(fig6$vary_cv[i, ],
  cat(sprintf("  %-22s n=%4d  CV %5.1f%% [%.1f, %.1f]\n", population, n_learners, cv, cv_lo, cv_hi)))

saveRDS(fig6, file.path(CACHE, "fig6_llm_slopes.rds"))
cat(sprintf("Augmented fig6_llm_slopes.rds: scaling (beta=%.3f, E=%.2f) + fig3 (n_child=%d, exemplars kappa=%s)\n",
            scaling_par$beta, scaling_par$E, dplyr::n_distinct(child_bg$ckey),
            paste(round(exk, 1), collapse = "/")))
