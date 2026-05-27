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

N_SIM_KIDS <- 200            # Mike's choice
AGE_GRID   <- seq(8, 36, by = 1)
QUANTILES  <- c(0.1, 0.25, 0.5, 0.75, 0.9)

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

  # Child RE structure
  has_int_re   <- "child" %in% names(vc)
  has_slope_re <- has_int_re &&
                  length(attr(vc$child, "stddev")) >= 2
  if (has_int_re) {
    sds <- attr(vc$child, "stddev")
    if (has_slope_re) {
      cor <- attr(vc$child, "correlation")[1, 2]
      Sigma <- diag(sds) %*% matrix(c(1, cor, cor, 1), 2, 2) %*% diag(sds)
      set.seed(1)
      kid_re <- MASS::mvrnorm(n_sim, mu = c(0, 0), Sigma = Sigma)
    } else {
      set.seed(1)
      kid_re <- cbind(rnorm(n_sim, 0, sds[1]), rep(0, n_sim))
    }
  } else {
    kid_re <- matrix(0, nrow = n_sim, ncol = 2)  # zero spread
    n_sim_use <- 1L                              # collapse to single curve
    kid_re <- kid_re[1, , drop = FALSE]
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
## form) vocab counts (sum of produces=1 over items they were scored
## on).
empirical_for <- function(lang_slug) {
  path <- file.path(LADDER_DIR, sprintf("data_%s.rds", lang_slug))
  if (!file.exists(path)) return(NULL)
  d <- readRDS(path)
  d$df |>
    group_by(child_id, age, form) |>
    summarise(vocab = sum(produces, na.rm = TRUE), .groups = "drop") |>
    mutate(language = d$language, lang_slug = lang_slug)
}

## ---- Loop over (lang, model) cells -----------------------------------
LANGS  <- c("english_american", "norwegian", "finnish",
            "french_quebecois", "japanese", "spanish_mexican",
            "french_french")
LANG_LABELS <- c(
  english_american = "English (American)",
  norwegian        = "Norwegian",
  finnish          = "Finnish",
  french_quebecois = "French (Quebecois)",
  japanese         = "Japanese",
  spanish_mexican  = "Spanish (Mexican)",
  french_french    = "French (French)"
)
MODELS <- c("A", "B_lin", "B_log", "C_lin", "C_log", "D_lin", "D_log")

cat("Building empirical scatters...\n")
emp <- bind_rows(lapply(LANGS, empirical_for))
cat(sprintf("Empirical rows: %d\n", nrow(emp)))

## summary CSVs (for AIC labels and a0 pivot)
summ <- list.files(LADDER_DIR, pattern = "^summary_.*\\.csv$",
                    full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |>
  bind_rows()

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
    pr <- predict_curves(fit, mid, AGE_GRID, n_sim = N_SIM_KIDS, a0 = a0)
    pr$lang_slug <- lang
    pr$model     <- mid
    all_preds[[paste(lang, mid, sep = "/")]] <- pr
    rm(fit); gc(verbose = FALSE)
    cat(sprintf("done (%.1f s)\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

preds <- bind_rows(all_preds)

## quantile summarize per (lang, model, age)
qtiles <- preds |>
  group_by(lang_slug, model, age) |>
  reframe(qprob = QUANTILES,
          vocab = quantile(vocab, QUANTILES))

write_csv(qtiles, OUT_CSV)
cat(sprintf("\nWrote %s\n", OUT_CSV))

## ---- Plot ------------------------------------------------------------
qtiles <- qtiles |>
  mutate(language  = factor(LANG_LABELS[lang_slug],
                              levels = LANG_LABELS[LANGS]),
          model     = factor(model, levels = MODELS),
          family    = case_when(model == "A" ~ "Pure accumulator",
                                grepl("lin$", model) ~ "Linear age",
                                grepl("log$", model) ~ "Log age"))
emp_lab <- emp |>
  mutate(language = factor(LANG_LABELS[lang_slug],
                            levels = LANG_LABELS[LANGS]))

# AIC labels per cell, relative to within-lang best
aic_labs <- summ |>
  group_by(language, lang_slug) |>
  mutate(AIC_best = min(AIC), dAIC = AIC - AIC_best) |>
  ungroup() |>
  mutate(language  = factor(language, levels = LANG_LABELS[LANGS]),
          model     = factor(model, levels = MODELS),
          label     = ifelse(dAIC == 0, "BEST",
                              sprintf("+%.0fK", dAIC / 1000)))

p <- ggplot() +
  # empirical data
  geom_point(data = emp_lab,
              aes(x = age, y = vocab),
              colour = "grey40", alpha = 0.25, size = 0.5) +
  # model quantile lines
  geom_line(data = qtiles,
             aes(x = age, y = vocab,
                 group = interaction(qprob, family),
                 colour = family,
                 alpha = ifelse(qprob == 0.5, 1.0, 0.6),
                 linewidth = ifelse(qprob == 0.5, 0.8, 0.4))) +
  facet_grid(language ~ model,
              scales = "free_y", space = "fixed") +
  geom_text(data = aic_labs,
             aes(x = -Inf, y = Inf, label = label),
             hjust = -0.1, vjust = 1.3, size = 2.6,
             colour = "grey25") +
  scale_colour_manual(values = c(`Pure accumulator` = "grey50",
                                  `Linear age`       = "#d7191c",
                                  `Log age`          = "#1f78b4")) +
  scale_alpha_identity() +
  scale_linewidth_identity() +
  labs(x = "Age (months)",
       y = "Productive vocabulary",
       title = "glmer ladder predictions vs empirical vocab(age) — 7 languages × 7 models",
       subtitle = sprintf("Lines = 10/25/50/75/90 quantiles across %d simulated kids drawn from each model's child-RE distribution. Grey dots = data. Corner label = ΔAIC vs best in row.",
                          N_SIM_KIDS),
       colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        strip.text    = element_text(size = 8, face = "bold"),
        legend.position = "top",
        panel.spacing = unit(0.4, "lines"))

ggsave(OUT_PNG, p, width = 18, height = 18, dpi = 180)
cat(sprintf("\nWrote %s\n", OUT_PNG))
