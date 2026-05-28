## glmer ladder — SIMULATION step (slow).
##
## For each (lang, model) fit, bootstrap N_SIM_KIDS "fake kids" from the
## model's per-kid BLUP distribution, compute each one's vocab curve
## (sum of inv_logit linear predictor over the main-form items at each
## age), and summarise 10/25/50/75/90 quantiles per (lang, model, age).
##
## Also collects the empirical spaghetti data and AIC summaries.
##
## Everything needed for plotting is saved to a cache RDS so the plot
## step (04b_plot.R) can iterate on aesthetics without re-simulating.
##
## Inputs:  fits/glmer_ladder/{data,fit}_<lang>_<model>.rds
## Output:  outputs/glmer_ladder/sim_cache.rds
##          outputs/glmer_ladder/predictions_quantiles.csv

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(lme4)
})

LADDER_DIR <- file.path(PATHS$fits_dir, "glmer_ladder")
OUT_CACHE  <- file.path(PATHS$outputs_dir, "glmer_ladder", "sim_cache.rds")
OUT_CSV    <- file.path(PATHS$outputs_dir, "glmer_ladder",
                        "predictions_quantiles.csv")
dir.create(dirname(OUT_CACHE), recursive = TRUE, showWarnings = FALSE)

N_SIM_KIDS <- 500            # smooths top-quantile sampling noise
QUANTILES  <- c(0.1, 0.25, 0.5, 0.75, 0.9)
## Age grid is per-language (5th–95th pct of observed ages) to avoid
## extrapolating predictions beyond the data window.

## ---- Core simulator -------------------------------------------------
## Bootstrap from the actual BLUPs (shrunken empirical-Bayes per-kid
## estimates) rather than sampling MVN(0, Sigma_hat). BLUP bootstrap
## gives predictions conditional on the *observed* kid distribution;
## parametric MVN sampling implies a wider (unshrunken) population and
## produces extreme (intercept, slope) combos the data never showed,
## which distort upper quantiles at thin-data ages.
predict_curves <- function(fit, model_id, ages, n_sim = 500, a0 = NA_real_,
                            keep_items = NULL) {
  fe <- fixef(fit)
  vc <- VarCorr(fit)

  intercept <- if ("(Intercept)" %in% names(fe)) fe[["(Intercept)"]] else 0
  age_coef  <- if ("log_age" %in% names(fe)) fe[["log_age"]]
                else if ("age_c" %in% names(fe)) fe[["age_c"]] else NA_real_

  uses_log <- model_id == "A" || grepl("_log$", model_id)
  age_term <- if (uses_log) log(ages / a0) else ages - a0
  if (model_id == "A") age_coef <- 1   # κ ≡ 1 (offset pins slope at 1)

  re_item   <- ranef(fit)$item
  if (!is.null(re_item) && !is.null(keep_items)) {
    re_item <- re_item[rownames(re_item) %in% keep_items, , drop = FALSE]
  }
  item_blups <- if (!is.null(re_item) && nrow(re_item) > 0) re_item[, 1] else 0

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

  out <- tibble(sim = integer(), age = numeric(), vocab = numeric())
  for (s in seq_len(n_use)) {
    kid_int   <- kid_re[s, 1]
    kid_slope <- kid_re[s, 2]
    eta_kid_age <- intercept + (age_coef + kid_slope) * age_term + kid_int
    M <- outer(eta_kid_age, item_blups, FUN = "+")
    vocab_age <- rowSums(plogis(M))
    out <- bind_rows(out, tibble(sim = s, age = ages, vocab = vocab_age))
  }
  out
}

## ---- Main form per language (largest item set) ----------------------
main_form_for <- function(lang_slug) {
  d <- readRDS(file.path(LADDER_DIR, sprintf("data_%s.rds", lang_slug)))
  d$df |>
    distinct(form, item) |>
    count(form, name = "n_items") |>
    slice_max(n_items, n = 1, with_ties = FALSE) |>
    pull(form)
}
main_form_items <- function(lang_slug) {
  d <- readRDS(file.path(LADDER_DIR, sprintf("data_%s.rds", lang_slug)))
  d$df |>
    filter(form == MAIN_FORM[[lang_slug]]) |>
    distinct(item) |>
    pull(item)
}

## ---- Empirical scatter (main form only) -----------------------------
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
            form = main_form, n_items = n_main_items)
}

## ---- Languages + models ---------------------------------------------
## Spanish (Mexican) excluded: WS-only, narrow age window, n=119 →
## degenerate D fits. We still simulate the 6 comparable languages.
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
cat("Main form per language:\n"); print(MAIN_FORM)

cat("\nBuilding empirical spaghetti data...\n")
emp <- bind_rows(lapply(LANGS, empirical_for))
cat(sprintf("Empirical rows: %d\n", nrow(emp)))

summ <- list.files(LADDER_DIR, pattern = "^summary_.*\\.csv$",
                    full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |>
  bind_rows() |>
  filter(lang_slug %in% LANGS)

age_ranges <- emp |>
  group_by(lang_slug) |>
  summarise(age_min = floor(quantile(age, 0.05, na.rm = TRUE)),
            age_max = ceiling(quantile(age, 0.95, na.rm = TRUE)),
            .groups = "drop")
cat("\nAge ranges per language (5th–95th pct):\n"); print(age_ranges)

## ---- Simulate -------------------------------------------------------
cat("\nLoading fits and simulating...\n")
all_preds <- list()
for (lang in LANGS) {
  for (mid in MODELS) {
    fit_path <- file.path(LADDER_DIR, sprintf("fit_%s_%s.rds", lang, mid))
    if (!file.exists(fit_path)) {
      cat(sprintf("  [skip] %s/%s — fit not found\n", lang, mid)); next
    }
    cat(sprintf("  %s / %s ... ", lang, mid)); t0 <- Sys.time()
    fit <- readRDS(fit_path)
    a0 <- summ$a0[summ$lang_slug == lang & summ$model == mid][1]
    if (is.na(a0)) a0 <- median(emp$age[emp$lang_slug == lang], na.rm = TRUE)
    ar  <- age_ranges[age_ranges$lang_slug == lang, ]
    ages_lang <- seq(ar$age_min, ar$age_max, by = 1)
    items_main <- main_form_items(lang)
    pr <- predict_curves(fit, mid, ages_lang, n_sim = N_SIM_KIDS, a0 = a0,
                          keep_items = items_main)
    pr$lang_slug <- lang; pr$model <- mid
    all_preds[[paste(lang, mid, sep = "/")]] <- pr
    rm(fit); gc(verbose = FALSE)
    cat(sprintf("done (%.1f s)\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
preds <- bind_rows(all_preds)

## clamp to (0, n_items) per language
n_items_lang <- emp |> distinct(lang_slug, n_items)
preds <- preds |>
  left_join(n_items_lang, by = "lang_slug") |>
  mutate(vocab = pmin(pmax(vocab, 0), n_items))

## quantiles per (lang, model, age)
qtiles <- preds |>
  group_by(lang_slug, model, age) |>
  reframe(qprob = QUANTILES, vocab = quantile(vocab, QUANTILES))

write_csv(qtiles, OUT_CSV)
cat(sprintf("\nWrote %s\n", OUT_CSV))

## ---- Save cache for plotting ----------------------------------------
saveRDS(list(
  qtiles      = qtiles,
  emp         = emp,
  summ        = summ,
  LANGS       = LANGS,
  LANG_LABELS = LANG_LABELS,
  MODELS      = MODELS,
  N_SIM_KIDS  = N_SIM_KIDS,
  age_ranges  = age_ranges,
  MAIN_FORM   = MAIN_FORM
), OUT_CACHE)
cat(sprintf("Wrote %s\n", OUT_CACHE))
cat("\nNow run: Rscript glmer_ladder/04b_plot.R\n")
