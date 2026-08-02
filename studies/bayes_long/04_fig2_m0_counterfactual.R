## 04_fig2_m0_counterfactual.R -- is Fig 2's decline evidence for acceleration?
##
## THE OBJECTION. Fig 2 plots estimated cumulative exposures at acquisition,
##   N_j = R * t50_j * p_j,   t50_j = a0 * exp((delta_j - log_H - mu_xi)/kappa),
## and the text says a pure accumulator "would show a constant pattern, with the same
## average number of instances heard before learning each word". That does not follow
## from our M0. M0 fixes kappa = 1 but leaves delta_j FREE, and free delta_j is precisely
## a free per-word bucket size. N_j is constant only if exp(delta_j) is proportional to
## 1/p_j -- i.e. only if difficulty is entirely frequency-determined, which M0 does not
## assume. So the decline could instead reflect later-acquired words simply being rarer.
##
## THE TEST. Rebuild the same figure from the FITTED M0 (its own delta_j and mu_xi, with
## kappa = 1) and compare the slope to M3's. Both use the identical frequency vector, so
## any difference is attributable to the item model rather than to the corpus.
## If M0 reproduces the decline, Fig 2 is not evidence for acceleration and should be
## reframed descriptively. If M0 is flat or much shallower, Fig 2 is a real model
## implication. No new fits are needed: the M0 rung of the ladder already exports psi.
##
## Usage:  Rscript studies/bayes_long/04_fig2_m0_counterfactual.R
## Output: paper/cache/si_fig2_m0.rds + console report

suppressPackageStartupMessages({library(dplyr); library(here)})
SUMM <- here("fits", "bayes_long", "summaries"); BL <- here("fits", "bayes_long")
EN   <- c("thal", "smith", "marchman")          # Fig 2 is the English pool
CACHE <- here("paper", "cache")

## sample-weighted pooling, mirroring paper/build_cache_short.R section 2
grab <- function(model) {
  parts <- lapply(EN, function(slug) {
    b <- readRDS(file.path(BL, sprintf("bundle_%s_a3.rds", slug)))
    s <- as.data.frame(readRDS(file.path(SUMM, sprintf("%s_a3_%s.summary.rds", slug, model))))
    g <- function(v) { x <- s$median[s$variable == v]; if (!length(x)) NA_real_ else x }
    ## M0 has no kappa_pop: acceleration is fixed at 1 by construction.
    kap <- if (model == "m0") 1 else g("kappa_pop")
    list(w = b$meta$n_kids, mu_xi = g("mu_xi"), kappa = kap,
         log_H = b$stan_data$log_H, a0 = b$stan_data$a0,
         psi = read.csv(file.path(SUMM, sprintf("%s_a3_%s_psi.csv", slug, model))) |>
                 transmute(item, delta_j))
  })
  W <- vapply(parts, `[[`, numeric(1), "w"); wsum <- sum(W)
  wavg <- function(f) sum(W * vapply(parts, `[[`, numeric(1), f)) / wsum
  psi <- bind_rows(Map(function(p, w) mutate(p$psi, w = w), parts, W)) |>
    group_by(item) |> filter(n() == length(EN)) |>
    summarise(delta_j = weighted.mean(delta_j, w), .groups = "drop")
  list(mu_xi = wavg("mu_xi"), kappa = wavg("kappa"),
       log_H = wavg("log_H"), a0 = wavg("a0"), psi = psi, n_kids = wsum)
}

m3 <- grab("m3"); m0 <- grab("m0")
cat(sprintf("M3: kappa=%.2f mu_xi=%.2f | M0: kappa=%.2f mu_xi=%.2f | items %d/%d\n",
            m3$kappa, m3$mu_xi, m0$kappa, m0$mu_xi, nrow(m3$psi), nrow(m0$psi)))

## frequency + lexical class, exactly as Fig 2 builds them
fig2 <- readRDS(file.path(CACHE, "fig2_efficiency.rds"))$items |>
  select(item, word, lexical_class, prob)
R_MID <- readRDS(file.path(CACHE, "fig2_efficiency.rds"))$anchor$mid

mk <- function(m, label) {
  inner_join(m$psi, fig2, by = "item") |>
    mutate(model = label,
           t_50 = m$a0 * exp((delta_j - m$log_H - m$mu_xi) / m$kappa),
           N_word = R_MID * t_50 * prob) |>
    filter(is.finite(t_50), is.finite(N_word), N_word > 0)
}
d3 <- mk(m3, "M3 (accelerating)"); d0 <- mk(m0, "M0 (pure accumulator)")

## TWO THINGS MUST BE MEASURED SEPARATELY.
##
## (a) The decline in Fig 2 is a WITHIN-LEXICAL-CLASS effect. Pooled across all 548 items
##     the slope is essentially flat (+0.007 log10/month, R2 = 0.001); it is only within
##     nouns, action words, descriptive words and function words that N falls with AoA,
##     with the class offsets cancelling in the pool. Any comparison to M0 therefore has
##     to be made class by class, which is also how the figure draws its fit lines.
##
## (b) A raw per-month slope is NOT comparable across the two item models, because they
##     imply very different AoA ranges: kappa divides the exponent, so M0 stretches the
##     same delta_j spread over a far wider span of months and mechanically flattens any
##     per-month slope. The scale-free quantity is the ELASTICITY d log10(N)/d log10(t50),
##     which equals 1 + d log10(p)/d log10(t50) under either model. That is what we compare.
##     Note what elasticity = 1 means: N rising in direct proportion to AoA, i.e. exposures
##     needed INCREASING with age -- the opposite of the paper's claim.
CLASSES <- levels(fig2$lexical_class)
slopes <- function(d, label) {
  bind_rows(lapply(CLASSES, function(cl) {
    dd <- filter(d, lexical_class == cl, is.finite(N_word), N_word > 0, t_50 > 0)
    if (nrow(dd) < 15) return(NULL)
    f_lin <- lm(log10(N_word) ~ t_50, data = dd)          # per-month, as drawn
    f_ela <- lm(log10(N_word) ~ log10(t_50), data = dd)   # scale-free
    data.frame(model = label, lexical_class = cl, n = nrow(dd),
               slope_per_month = unname(coef(f_lin)[2]), r2_lin = summary(f_lin)$r.squared,
               elasticity = unname(coef(f_ela)[2]), r2_ela = summary(f_ela)$r.squared)
  }))
}
S <- bind_rows(slopes(d3, "M3 (accelerating)"), slopes(d0, "M0 (pure accumulator)"))

cat("\n=== AoA range each item model implies (all items) ===\n")
for (d in list(d3, d0)) cat(sprintf("  %-22s t50 %.1f - %.1f months (median %.1f)\n",
    d$model[1], min(d$t_50), max(d$t_50), median(d$t_50)))
cat("  [the data span roughly 8-36 months; an item model implying acquisition at 254\n")
cat("   months is not reproducing the observed range]\n")

cat("\n=== pooled across all items (what the figure does NOT show) ===\n")
for (d in list(d3, d0)) {
  f <- lm(log10(N_word) ~ log10(t_50), data = d)
  cat(sprintf("  %-22s elasticity %+.2f  (R2 %.3f)\n", d$model[1], coef(f)[2], summary(f)$r.squared))
}

cat("\n=== BY LEXICAL CLASS: elasticity d log10(N) / d log10(AoA) ===\n")
cat(sprintf("  %-20s %22s %22s\n", "", "M3 (accelerating)", "M0 (pure accumulator)"))
for (cl in CLASSES) {
  a <- S[S$model == "M3 (accelerating)"     & S$lexical_class == cl, ]
  b <- S[S$model == "M0 (pure accumulator)" & S$lexical_class == cl, ]
  if (!nrow(a) || !nrow(b)) next
  cat(sprintf("  %-20s %+8.2f (R2 %.2f, n=%3d) %+8.2f (R2 %.2f, n=%3d)\n",
              cl, a$elasticity, a$r2_ela, a$n, b$elasticity, b$r2_ela, b$n))
}

cat("\n=== how to read this ===\n")
cat("  Negative elasticity  = exposures needed FALL with age (the paper's claim).\n")
cat("  Elasticity near +1   = exposures needed RISE in proportion to age.\n")
cat("  If M0 reproduces M3's negative elasticities, Fig 2 does not discriminate the\n")
cat("  models and must be reframed descriptively; if M0 is far less negative, it does.\n")

saveRDS(list(m3 = d3, m0 = d0, slopes = S), file.path(CACHE, "si_fig2_m0.rds"))
cat("\nwrote", file.path(CACHE, "si_fig2_m0.rds"), "\n")
