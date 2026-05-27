## 7-language × 7-model mega-figure: empirical vocab(age) data with
## glmer-predicted quantile fans for each ladder model.
##
## Approach: for each (lang, model) fit, simulate 200 fake kids from
## the model's implied kid-RE distribution, compute each one's vocab
## curve (sum of inv_logit linear predictor over all items at each age),
## then plot 10/25/50/75/90 quantiles across fake kids.
##
## For A and B models (no kid RE), the implied "spread" is zero — we
## plot a single population line.
##
## Inputs: fits/glmer_ladder/{data,fit}_<lang>_<model>.rds  (need fit RDS
##         not just summary; ranef(fit)$item is required and isn't in
##         the per-kid ranef CSV).
## Output: outputs/figs/longitudinal/glmer_ladder_mega.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(lme4); library(MASS); library(patchwork)
})

LADDER_DIR <- file.path(PATHS$fits_dir, "glmer_ladder")
OUT_PNG    <- file.path(PATHS$figs_dir, "longitudinal",
                        "glmer_ladder_mega.png")
OUT_CSV    <- file.path(PATHS$outputs_dir, "glmer_ladder",
                        "predictions_quantiles.csv")
dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

N_SIM_KIDS <- 500            # smooths top-quantile sampling noise
QUANTILES  <- c(0.1, 0.25, 0.5, 0.75, 0.9)
## AGE_GRID is now per-language (set inside the loop) so we don't
## extrapolate predictions beyond the observed data window.

## ---- Core simulator -------------------------------------------------
## Given a fitted glmer, simulate N_SIM_KIDS "fake kids" from the
## model's kid-RE distribution. For each (sim, age), compute predicted
## vocab = sum_j inv_logit(eta_ij). Return long data with quantiles.
predict_curves <- function(fit, model_id, ages, n_sim = 200, a0 = NA_real_) {
  fe <- fixef(fit)
  vc <- VarCorr(fit)

  intercept <- if ("(Intercept)" %in% names(fe)) fe[["(Intercept)"]] else 0
  age_coef  <- if ("log_age" %in% names(fe)) fe[["log_age"]]
                else if ("age_c" %in% names(fe)) fe[["age_c"]] else NA_real_

  # Determine age scaling used by the formula
  uses_log <- model_id == "A" || grepl("_log$", model_id)
  age_term <- if (uses_log) log(ages / a0) else ages - a0
  # Model A: pure unit accumulator (κ ≡ 1, no fitted slope). For A,
  # there's no age_coef in fixef, but the offset pins the log_age
  # coefficient at 1, so effective slope = 1.
  if (model_id == "A") age_coef <- 1

  # Item BLUPs
  re_item   <- ranef(fit)$item
  item_blups <- if (!is.null(re_item)) re_item[, 1] else 0

  # Child RE structure. Bootstrap from the actual BLUPs (the shrunken
  # empirical-Bayes per-kid estimates) rather than sampling from
  # MVN(0, Sigma_hat). Bootstrapping BLUPs gives the model's
  # predictions conditional on the *observed* kid distribution, which
  # is what we want when comparing model curves to data. Parametric
  # MVN sampling implies a wider population (its variance components
  # are unshrunken) and produces extreme (intercept, slope) combos
  # that the data never showed — those distort the upper-quantile
  # lines, especially at thin-data ages.
  has_int_re   <- "child" %in% names(vc)
  re_child_blup <- if (has_int_re) ranef(fit)$child else NULL
  if (has_int_re && nrow(re_child_blup) > 0) {
    set.seed(1)
    idx <- sample.int(nrow(re_child_blup), n_sim, replace = TRUE)
    int_col   <- re_child_blup[idx, 1]
    slope_col <- if (ncol(re_child_blup) >= 2) re_child_blup[idx, 2]
                  else rep(0, n_sim)
    kid_re <- cbind(int_col, slope_col)
  } else {
    kid_re <- matrix(0, nrow = 1, ncol = 2)  # collapse to single curve
  }
  n_use <- nrow(kid_re)

  # For each sim, compute vocab per age
  J <- length(item_blups)
  out <- tibble(sim = integer(), age = numeric(), vocab = numeric())
  for (s in seq_len(n_use)) {
    kid_int   <- kid_re[s, 1]
    kid_slope <- kid_re[s, 2]
    # eta_base over ages: vector length(ages)
    eta_kid_age <- intercept + (age_coef + kid_slope) * age_term + kid_int
    # vocab per age = sum_j inv_logit(eta_kid_age[a] + item_blups[j])
    # vectorized: outer addition then plogis then row-sum
    M <- outer(eta_kid_age, item_blups, FUN = "+")
    vocab_age <- rowSums(plogis(M))
    out <- bind_rows(out,
                      tibble(sim = s, age = ages, vocab = vocab_age))
  }
  out
}

## ---- Build empirical scatter -----------------------------------------
## For each language, load the bundle's df and compute per-(kid, age,
## form) vocab counts.
##
## Note on de-duplication: a few kids in wordbank have multiple admin_id
## records at the same (kid, age, form) — they get pulled as duplicate
## rows in get_instrument_data. We dedupe by (child_id, age, form, item)
## first so a kid's vocab is bounded by n_items. (The model fits don't
## dedupe — flagged for v2 cleanup.)
empirical_for <- function(lang_slug) {
  path <- file.path(LADDER_DIR, sprintf("data_%s.rds", lang_slug))
  if (!file.exists(path)) return(NULL)
  d <- readRDS(path)
  d$df |>
    distinct(child_id, age, form, item, .keep_all = TRUE) |>
    group_by(child_id, age, form) |>
    summarise(vocab = sum(produces, na.rm = TRUE), .groups = "drop") |>
    mutate(language = d$language, lang_slug = lang_slug,
            n_items   = d$n_items)
}

## ---- Loop over (lang, model) cells -----------------------------------
## Spanish (Mexican) excluded: WS-only, narrow age window (17-30 mo),
## n=119 — fixed-age effects identification is degenerate (B/D slopes
## go negative because all variation is between-kid), and D_log
## optimizer blows σ_slope up to 43. Not comparable to the others.
LANGS  <- c("english_american", "norwegian", "finnish",
            "french_quebecois", "japanese", "french_french")
LANG_LABELS <- c(
  english_american = "English (American)",
  norwegian        = "Norwegian",
  finnish          = "Finnish",
  french_quebecois = "French (Quebecois)",
  japanese         = "Japanese",
  french_french    = "French (French)"
)
MODELS <- c("A", "B_lin", "B_log", "C_lin", "C_log", "D_lin", "D_log")

cat("Building empirical scatters...\n")
emp <- bind_rows(lapply(LANGS, empirical_for))
cat(sprintf("Empirical rows: %d\n", nrow(emp)))

## summary CSVs (for AIC labels and a0 pivot) — filter to our LANGS
summ <- list.files(LADDER_DIR, pattern = "^summary_.*\\.csv$",
                    full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |>
  bind_rows() |>
  filter(lang_slug %in% LANGS)

## Per-language age range: use 5th–95th percentile of observed admin
## ages, not min/max. The thin tails (admins outside this range) are
## sparse enough that extrapolating BLUPs to them produces extreme
## young-age predictions for high-intercept / low-slope kids that the
## data doesn't actually contain.
age_ranges <- emp |>
  group_by(lang_slug) |>
  summarise(age_min = floor(quantile(age, 0.05, na.rm = TRUE)),
            age_max = ceiling(quantile(age, 0.95, na.rm = TRUE)),
            .groups = "drop")
cat("\nAge ranges per language (5th–95th percentile):\n"); print(age_ranges)

## predictions
cat("\nLoading fits and simulating...\n")
all_preds <- list()
for (lang in LANGS) {
  for (mid in MODELS) {
    fit_path <- file.path(LADDER_DIR, sprintf("fit_%s_%s.rds", lang, mid))
    if (!file.exists(fit_path)) {
      cat(sprintf("  [skip] %s/%s — fit not found\n", lang, mid))
      next
    }
    cat(sprintf("  %s / %s ... ", lang, mid)); t0 <- Sys.time()
    fit <- readRDS(fit_path)
    a0 <- summ$a0[summ$lang_slug == lang & summ$model == mid][1]
    if (is.na(a0)) a0 <- median(emp$age[emp$lang_slug == lang],
                                  na.rm = TRUE)
    ar  <- age_ranges[age_ranges$lang_slug == lang, ]
    ages_lang <- seq(ar$age_min, ar$age_max, by = 1)
    pr <- predict_curves(fit, mid, ages_lang, n_sim = N_SIM_KIDS, a0 = a0)
    pr$lang_slug <- lang
    pr$model     <- mid
    all_preds[[paste(lang, mid, sep = "/")]] <- pr
    rm(fit); gc(verbose = FALSE)
    cat(sprintf("done (%.1f s)\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

preds <- bind_rows(all_preds)

## Clamp predicted vocab to (0, n_items) per language — the model
## itself can't predict more items than exist, and pathological D fits
## (Spanish-MX) can produce simulated curves that briefly explode.
n_items_lang <- emp |> distinct(lang_slug, n_items)
preds <- preds |>
  left_join(n_items_lang, by = "lang_slug") |>
  mutate(vocab = pmin(pmax(vocab, 0), n_items))

## quantile summarize per (lang, model, age) — long form for plotting
qtiles <- preds |>
  group_by(lang_slug, model, age) |>
  reframe(qprob = QUANTILES,
          vocab = quantile(vocab, QUANTILES))

write_csv(qtiles, OUT_CSV)
cat(sprintf("\nWrote %s\n", OUT_CSV))

## ---- Plot ------------------------------------------------------------
qtiles <- qtiles |>
  mutate(language = factor(LANG_LABELS[lang_slug], levels = LANG_LABELS[LANGS]),
          model    = factor(model, levels = MODELS))
emp_lab <- emp |>
  mutate(language = factor(LANG_LABELS[lang_slug], levels = LANG_LABELS[LANGS]))

# AIC labels per cell, relative to within-lang best
aic_labs <- summ |>
  group_by(language, lang_slug) |>
  mutate(AIC_best = min(AIC), dAIC = AIC - AIC_best) |>
  ungroup() |>
  mutate(language  = factor(language, levels = LANG_LABELS[LANGS]),
          model     = factor(model, levels = MODELS),
          label     = ifelse(dAIC == 0, "best",
                              sprintf("Δ%+.0fk", dAIC / 1000)))

# With Spanish-MX excluded, no remaining cells have pathological
# σ_slope. bad_fits is empty; we drop the annotation layers entirely.

## quantile lines: factor with wordbank palette
qtiles <- qtiles |>
  mutate(qprob_f = factor(qprob, levels = c(0.1, 0.25, 0.5, 0.75, 0.9),
                           labels = c("0.1", "0.25", "0.5", "0.75", "0.9")))
WORDBANK_PALETTE <- c("0.1"  = "#1f78b4",   # dark blue
                       "0.25" = "#a6cee3",   # light blue
                       "0.5"  = "#33a02c",   # green
                       "0.75" = "#fdbf6f",   # gold
                       "0.9"  = "#e31a1c")   # red

p <- ggplot() +
  # empirical data
  geom_point(data = emp_lab,
              aes(x = age, y = vocab),
              colour = "grey25", alpha = 0.3, size = 0.4) +
  # quantile lines for models with kid RE (5 lines per panel)
  geom_line(data = qtiles |> filter(model %in% c("C_lin", "C_log",
                                                  "D_lin", "D_log")),
             aes(x = age, y = vocab,
                 colour = qprob_f, group = qprob_f,
                 linewidth = ifelse(as.character(qprob_f) == "0.5", 0.9, 0.55))) +
  # A and B collapse to a single curve — show only the median in green
  # so we don't stack 5 identical lines and end up showing only the
  # last-drawn one.
  geom_line(data = qtiles |> filter(!model %in% c("C_lin", "C_log",
                                                   "D_lin", "D_log"),
                                     qprob_f == "0.5"),
             aes(x = age, y = vocab, colour = qprob_f),
             linewidth = 0.95) +
  facet_grid(language ~ model, scales = "free_y", space = "fixed") +
  geom_text(data = aic_labs,
             aes(x = -Inf, y = Inf, label = label),
             hjust = -0.15, vjust = 1.5, size = 2.6, colour = "grey25") +
  scale_colour_manual(values = WORDBANK_PALETTE, name = "Percentile") +
  scale_linewidth_identity() +
  labs(x = "Age (months)",
       y = "Productive vocabulary",
       title = "glmer ladder: predictions vs empirical vocab(age) — 6 languages × 7 models",
       subtitle = sprintf("Lines = 10/25/50/75/90 quantiles across %d simulated kids drawn from each model's child-RE distribution. Dots = data. Corner label = ΔAIC vs best within language.",
                          N_SIM_KIDS)) +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        strip.text    = element_text(size = 8, face = "bold"),
        legend.position = "top",
        panel.spacing = unit(0.4, "lines"))

ggsave(OUT_PNG, p, width = 18, height = 16, dpi = 180)
cat(sprintf("\nWrote %s\n", OUT_PNG))
