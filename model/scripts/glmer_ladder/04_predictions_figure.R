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
predict_curves <- function(fit, model_id, ages, n_sim = 200, a0 = NA_real_,
                            keep_items = NULL) {
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

  # Item BLUPs — optionally subset to a given set of items so that we
  # sum predicted vocab only over items present on a specific form
  # (matches the form the empirical scatter shows).
  re_item   <- ranef(fit)$item
  if (!is.null(re_item) && !is.null(keep_items)) {
    re_item <- re_item[rownames(re_item) %in% keep_items, , drop = FALSE]
  }
  item_blups <- if (!is.null(re_item) && nrow(re_item) > 0) re_item[, 1] else 0

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

## ---- Pick a single "main" form per language --------------------------
## When languages have multiple forms (WG/WS, plus Short variants), each
## admin's empirical vocab is capped by the items on whichever form the
## kid took (e.g., Finnish WSShort has only 100 items). To compare
## predictions to data on the same scale, restrict the figure to one
## form per language. Choose the form with the most items (most
## informative; usually WS or the long WS variant).
main_form_for <- function(lang_slug) {
  d <- readRDS(file.path(LADDER_DIR, sprintf("data_%s.rds", lang_slug)))
  d$df |>
    distinct(form, item) |>
    count(form, name = "n_items") |>
    slice_max(n_items, n = 1, with_ties = FALSE) |>
    pull(form)
}

## items on the main form for each language (used to subset prediction
## sum and to know the ceiling). MAIN_FORM is initialized after LANGS
## is defined below.
main_form_items <- function(lang_slug) {
  d <- readRDS(file.path(LADDER_DIR, sprintf("data_%s.rds", lang_slug)))
  d$df |>
    filter(form == MAIN_FORM[[lang_slug]]) |>
    distinct(item) |>
    pull(item)
}

## ---- Build empirical scatter (main form only) ------------------------
empirical_for <- function(lang_slug) {
  path <- file.path(LADDER_DIR, sprintf("data_%s.rds", lang_slug))
  if (!file.exists(path)) return(NULL)
  d <- readRDS(path)
  main_form <- MAIN_FORM[[lang_slug]]
  n_main_items <- length(main_form_items(lang_slug))
  d$df |>
    filter(form == main_form) |>
    distinct(child_id, age, item, .keep_all = TRUE) |>
    group_by(child_id, age) |>
    summarise(vocab = sum(produces, na.rm = TRUE), .groups = "drop") |>
    mutate(language = d$language, lang_slug = lang_slug,
            form      = main_form,
            n_items   = n_main_items)
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

MAIN_FORM <- setNames(sapply(LANGS, main_form_for), LANGS)
cat("\nMain form per language:\n"); print(MAIN_FORM)

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
    items_main <- main_form_items(lang)
    pr <- predict_curves(fit, mid, ages_lang, n_sim = N_SIM_KIDS, a0 = a0,
                          keep_items = items_main)
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
  mutate(qprob_f = factor(qprob, levels = c(0.1, 0.25, 0.5, 0.75, 0.9),
                           labels = c("0.1", "0.25", "0.5", "0.75", "0.9")))
WORDBANK_PALETTE <- c("0.1"  = "#1f78b4",   # dark blue
                       "0.25" = "#a6cee3",   # light blue
                       "0.5"  = "#33a02c",   # green
                       "0.75" = "#fdbf6f",   # gold
                       "0.9"  = "#e31a1c")   # red

## Build the plot for a given language subset. Used twice: once for all
## 6 languages (supplementary) and once for the 4 well-powered languages
## (main text).
build_ladder_plot <- function(langs_subset, out_png, title_n_langs,
                                width = 18, height = 16) {
  qt <- qtiles |>
    filter(lang_slug %in% langs_subset) |>
    mutate(language = factor(LANG_LABELS[lang_slug],
                              levels = LANG_LABELS[langs_subset]),
            model    = factor(model, levels = MODELS))
  el <- emp |>
    filter(lang_slug %in% langs_subset) |>
    mutate(language = factor(LANG_LABELS[lang_slug],
                              levels = LANG_LABELS[langs_subset]))
  al <- summ |>
    filter(lang_slug %in% langs_subset) |>
    group_by(language, lang_slug) |>
    mutate(AIC_best = min(AIC), dAIC = AIC - AIC_best) |>
    ungroup() |>
    mutate(language = factor(language,
                              levels = LANG_LABELS[langs_subset]),
            model    = factor(model, levels = MODELS),
            label    = ifelse(dAIC == 0, "best",
                              sprintf("Δ%+.0fk", dAIC / 1000)))

  p <- ggplot() +
    geom_line(data = el,
               aes(x = age, y = vocab, group = child_id),
               colour = "grey25", alpha = 0.15, linewidth = 0.25) +
    geom_point(data = el,
                aes(x = age, y = vocab),
                colour = "grey25", alpha = 0.2, size = 0.25) +
    geom_line(data = qt |> filter(model %in% c("C_lin", "C_log",
                                                 "D_lin", "D_log")),
               aes(x = age, y = vocab,
                   colour = qprob_f, group = qprob_f,
                   linewidth = ifelse(as.character(qprob_f) == "0.5",
                                       0.9, 0.55))) +
    geom_line(data = qt |> filter(!model %in% c("C_lin", "C_log",
                                                  "D_lin", "D_log"),
                                    qprob_f == "0.5"),
               aes(x = age, y = vocab, colour = qprob_f),
               linewidth = 0.95) +
    facet_grid(language ~ model, scales = "free_y", space = "fixed") +
    geom_text(data = al,
               aes(x = -Inf, y = Inf, label = label),
               hjust = -0.15, vjust = 1.5, size = 2.6, colour = "grey25") +
    scale_colour_manual(values = WORDBANK_PALETTE, name = "Percentile") +
    scale_linewidth_identity() +
    labs(x = "Age (months)",
         y = "Productive vocabulary",
         title = sprintf("glmer ladder: predictions vs empirical vocab(age) — %d languages × 7 models",
                          title_n_langs),
         subtitle = sprintf("Lines = 10/25/50/75/90 quantiles across %d simulated kids (BLUP bootstrap from each model's child-RE distribution). Grey lines = per-kid longitudinal trajectories. Empirical and predictions both restricted to the largest form per language. Corner label = ΔAIC vs best within language.",
                            N_SIM_KIDS)) +
    theme_minimal(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey25"),
          strip.text    = element_text(size = 8, face = "bold"),
          legend.position = "top",
          panel.spacing = unit(0.4, "lines"))

  ggsave(out_png, p, width = width, height = height, dpi = 180)
  cat(sprintf("Wrote %s\n", out_png))
}

## Main-text version: 4 well-powered languages
LANGS_MAIN <- c("english_american", "norwegian",
                "french_quebecois", "japanese")
build_ladder_plot(
  LANGS_MAIN,
  file.path(PATHS$figs_dir, "longitudinal", "glmer_ladder_main.png"),
  title_n_langs = length(LANGS_MAIN),
  width = 18, height = 11
)

## Supplementary version: all 6 languages
build_ladder_plot(
  LANGS,
  OUT_PNG,
  title_n_langs = length(LANGS),
  width = 18, height = 16
)
