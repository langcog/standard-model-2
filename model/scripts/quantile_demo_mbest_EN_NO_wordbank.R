## EN + NO M_best side-by-side quantile-fan plot, with empirical
## reference drawn from the WHOLE WORDBANK corpus (not just the
## bundle's I=500 stratified sample). This is the generalization
## check -- "does the model fit the broader population?" -- as
## opposed to the bundle-internal version in
## model/scripts/quantile_demo_mbest_EN_NO.R.
##
## Same simulation logic, same GAMLSS BEINF smoother for the empirical
## quantile fan (lambda = 10000, per MB-CDI manual approach).
##
## Output: outputs/figs/longitudinal/m_best_quantile_EN_NO_wordbank.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/empirical_xsec_helper.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(patchwork)
})

OUT_DIR    <- file.path(PATHS$figs_dir, "longitudinal")
SUMMARIES  <- file.path(PATHS$fits_dir, "summaries")
LONG_ITEMS <- file.path(PATHS$fits_dir, "long_items.rds")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

AGE_RANGE      <- c(16, 30)
TAUS           <- c(0.10, 0.25, 0.50, 0.75, 0.90)
SEED           <- 20260523
N_DRAWS_USE    <- 200
N_KIDS_PER_AGE <- 200
LAMBDA         <- 10000

PALETTE <- c("0.1" = "#1f78b4", "0.25" = "#a6cee3",
             "0.5" = "#33a02c", "0.75" = "#fdbf6f",
             "0.9" = "#e31a1c")

## ---- Generic per-language panel builder --------------------------
build_panel <- function(language_label,
                        bundle_path,
                        tag,
                        cache_path,
                        language_filter,
                        title) {
  bundle <- readRDS(bundle_path)
  sd_b   <- bundle$stan_data
  item_col <- if ("item" %in% names(bundle$word_info)) "item" else "item_definition"
  bundle_items <- bundle$word_info[[item_col]]
  J_ITEMS <- sd_b$J

  cat(sprintf("[%s] I=%d, J=%d\n", language_label, sd_b$I, J_ITEMS))

  ## Empirical: wordbank-wide cross-sectional, one admin per child.
  emp_df <- build_xsec_empirical_wordbank(
    item_definitions = bundle_items,
    language_filter  = language_filter,
    form_filter      = "WS",
    age_range        = AGE_RANGE,
    long_items_path  = LONG_ITEMS,
    cache_path       = cache_path
  )
  cat(sprintf("[%s] empirical (wordbank x-sec): %d kids, median %d items/kid\n",
              language_label, nrow(emp_df), median(emp_df$n_items)))

  ## Load draws
  draws <- readRDS(file.path(SUMMARIES, paste0(tag, ".draws.rds")))
  if ("draws_df" %in% class(draws)) draws <- as.data.frame(draws)
  if (!"sigma_s" %in% names(draws)) draws$sigma_s <- 0
  if (!"s" %in% names(draws)) draws$s <- 0
  delta_j_df <- read_csv(file.path(SUMMARIES, paste0(tag, "_psi.csv")),
                          show_col_types = FALSE)
  med_col <- if ("delta_j_median" %in% names(delta_j_df)) "delta_j_median" else "psi_median"
  delta_j <- delta_j_df[[med_col]]

  ## Simulate
  set.seed(SEED)
  draw_idx <- sort(sample(seq_len(nrow(draws)), N_DRAWS_USE))
  AGE_GRID <- seq(AGE_RANGE[1], AGE_RANGE[2], by = 0.5)
  log_H <- sd_b$log_H; a0 <- sd_b$a0
  mu_r  <- sd_b$mu_r;  sigma_r <- sd_b$sigma_r

  rows <- vector("list", length(draw_idx) * length(AGE_GRID))
  idx <- 0
  for (k in draw_idx) {
    sa <- as.numeric(draws$sigma_alpha[k])
    sz <- as.numeric(draws$sigma_zeta[k])
    dd <- as.numeric(draws$delta[k])
    sd_d <- as.numeric(draws$s[k])
    sigma_xi <- sqrt(sa^2 + sigma_r^2)
    for (a in AGE_GRID) {
      idx <- idx + 1
      set.seed(SEED + k * 1000 + as.integer(a * 10))
      xi   <- rnorm(N_KIDS_PER_AGE, mu_r, sigma_xi)
      zeta <- rnorm(N_KIDS_PER_AGE, 0, sz)
      kappa <- 1 + dd + zeta
      log_age <- log(pmax(a - sd_d, 0.01) / a0)
      theta <- xi + kappa * log_age
      eta <- outer(theta + log_H, delta_j, "-")
      rows[[idx]] <- data.frame(age = a, vocab = rowSums(plogis(eta)))
    }
  }
  sim_df <- bind_rows(rows)
  sim_preds <- sim_df |>
    group_by(age) |>
    summarise(
      `0.1`  = quantile(vocab, 0.10, na.rm = TRUE),
      `0.25` = quantile(vocab, 0.25, na.rm = TRUE),
      `0.5`  = quantile(vocab, 0.50, na.rm = TRUE),
      `0.75` = quantile(vocab, 0.75, na.rm = TRUE),
      `0.9`  = quantile(vocab, 0.90, na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_longer(-age, names_to = "tau", values_to = "vocab")

  ## Empirical quantile fan via GAMLSS BEINF, lambda = 10000.
  emp_pred <- fit_xsec_quantile_fan(
    emp_df, taus = TAUS,
    age_grid = seq(AGE_RANGE[1], AGE_RANGE[2], by = 0.25),
    J_total = J_ITEMS,
    lambda = LAMBDA
  )

  ## Scatter (stratified-by-age)
  set.seed(SEED)
  emp_scatter <- emp_df |>
    mutate(age_int = round(age)) |>
    group_by(age_int) |>
    group_modify(~ slice_sample(.x, n = min(nrow(.x), 80))) |>
    ungroup() |>
    select(-age_int)

  ## Plot
  ggplot(emp_scatter, aes(age, vocab)) +
    geom_jitter(width = 0.3, alpha = 0.20, size = 0.5, colour = "grey30") +
    geom_line(data = emp_pred,
              aes(age, vocab, group = tau),
              colour = "grey45", linetype = "dashed",
              linewidth = 0.5) +
    geom_line(data = sim_preds,
              aes(age, vocab, group = tau, colour = tau),
              linewidth = 1.1) +
    scale_colour_manual(values = PALETTE, name = "Percentile") +
    scale_x_continuous(breaks = c(16, 20, 24, 28)) +
    coord_cartesian(xlim = c(AGE_RANGE[1] - 0.4, AGE_RANGE[2] + 0.4),
                    ylim = c(0, J_ITEMS)) +
    labs(x = "Age (months)",
         y = sprintf("Productive vocabulary (of %d items)", J_ITEMS),
         title = title,
         subtitle = sprintf("%d kids (1 admin/kid; wordbank-wide x-sec)", nrow(emp_df))) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey25"),
          legend.position = "bottom")
}

p_en <- build_panel(
  language_label  = "EN",
  bundle_path     = file.path(PATHS$fits_dir, "long_subset_data.rds"),
  tag             = "long_no_freq_slopes",
  cache_path      = file.path(PATHS$fits_dir, "demo_xsec_empirical_wordbank_EN.rds"),
  language_filter = "English (American)",
  title           = "English"
)
p_no <- build_panel(
  language_label  = "NO",
  bundle_path     = file.path(PATHS$fits_dir, "long_subset_data_nor.rds"),
  tag             = "long_no_freq_slopes_norwegian",
  cache_path      = file.path(PATHS$fits_dir, "demo_xsec_empirical_wordbank_NO.rds"),
  language_filter = "Norwegian",
  title           = "Norwegian"
)

combined <- (p_en | p_no) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "M_best vs. wordbank-wide cross-sectional empirical",
    subtitle = "Model: longitudinal I=500 fit (α + ζ + δ, no s, no s_i). Empirical: full wordbank WS, one admin/kid (random), GAMLSS BEINF λ=10000."
  ) &
  theme(legend.position = "bottom")

out_png <- file.path(OUT_DIR, "m_best_quantile_EN_NO_wordbank.png")
ggsave(out_png, combined, width = 11, height = 5.5, dpi = 200)
cat(sprintf("\nWrote: %s\n", out_png))
