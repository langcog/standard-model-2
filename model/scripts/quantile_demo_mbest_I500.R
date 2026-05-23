## Single-panel quantile-fan plot for the I=500 M_best refit.
##
## Matches the slopes-only panel of model/scripts/quantile_demo.R so it
## can be compared visually to the I=200 pilot's 5-panel figure
## (outputs/figs/longitudinal/quantile_demo.png), without sharing axes
## or being added as a panel to it. Same simulation logic, same TAUS,
## same empirical reference (Wordbank WS, restricted to bundle items).
##
## Output: outputs/figs/longitudinal/m_best_quantile_I500.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/empirical_xsec_helper.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
})

OUT_DIR    <- file.path(PATHS$figs_dir, "longitudinal")
SUMMARIES  <- file.path(PATHS$fits_dir, "summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BUNDLE_PATH    <- file.path(PATHS$fits_dir, "long_subset_data.rds")
CACHE_FULL_EMP <- file.path(PATHS$fits_dir, "demo_xsec_empirical_long.rds")
LONG_ITEMS     <- file.path(PATHS$fits_dir, "long_items.rds")
AGE_RANGE      <- c(16, 30)
TAUS           <- c(0.10, 0.25, 0.50, 0.75, 0.90)
TAG            <- "long_no_freq_slopes"
SEED           <- 20260522
N_DRAWS_USE    <- 200      # posterior draws used for simulation
N_KIDS_PER_AGE <- 200      # simulated kids per posterior draw per age

bundle  <- readRDS(BUNDLE_PATH)
sd_b    <- bundle$stan_data
cat(sprintf("Bundle: I=%d, A=%d, J=%d, N=%d\n",
            sd_b$I, sd_b$A, sd_b$J, sd_b$N))
J_ITEMS <- sd_b$J
# Bundle's word_info uses `item` (the older `item_definition` column was
# renamed in the bundle pipeline). Fall back to `item_definition` if a
# legacy bundle is loaded.
item_col <- if ("item" %in% names(bundle$word_info)) {
  "item"
} else if ("item_definition" %in% names(bundle$word_info)) {
  "item_definition"
} else {
  stop("word_info has no item/item_definition column")
}
bundle_items <- bundle$word_info[[item_col]]
cat(sprintf("Bundle items: %d (column '%s')\n", length(bundle_items), item_col))

## ---- Empirical reference: cross-sectional, one admin per kid -----
## Drawn from the bundle's own admin-level data (bundle$df) so we have
## proper admin_key disambiguation. long_items.rds was previously used
## but doesn't track admin_id -- multiple admins at the same age for
## the same kid got summed together, inflating vocab beyond J.
emp_df <- build_xsec_empirical(
  bundle_df  = bundle$df,
  age_range  = AGE_RANGE,
  cache_path = CACHE_FULL_EMP
)
cat(sprintf("Empirical (x-sec): %d kids, median %d items per kid\n",
            nrow(emp_df), median(emp_df$n_items)))

## ---- Load M_best draws ---------------------------------------------
draws_path   <- file.path(SUMMARIES, paste0(TAG, ".draws.rds"))
delta_j_path <- file.path(SUMMARIES, paste0(TAG, "_psi.csv"))
draws        <- readRDS(draws_path)
if ("draws_df" %in% class(draws)) draws <- as.data.frame(draws)
# M_best has no sigma_s -- fill with 0 for sim consistency
if (!"sigma_s" %in% names(draws)) draws$sigma_s <- 0
# s is essentially pinned at 0 in this variant
if (!"s" %in% names(draws)) draws$s <- 0

delta_j_df <- read_csv(delta_j_path, show_col_types = FALSE)
median_col <- if ("delta_j_median" %in% names(delta_j_df)) "delta_j_median" else "psi_median"
delta_j <- delta_j_df[[median_col]]
cat(sprintf("Draws: %d posterior draws, delta_j: %d items\n",
            nrow(draws), length(delta_j)))

## ---- Simulate ------------------------------------------------------
set.seed(SEED)
draw_idx <- sort(sample(seq_len(nrow(draws)), N_DRAWS_USE))
AGE_GRID <- seq(AGE_RANGE[1], AGE_RANGE[2], by = 0.5)
log_H    <- sd_b$log_H
a0       <- sd_b$a0
mu_r     <- sd_b$mu_r
sigma_r  <- sd_b$sigma_r

cat(sprintf("Simulating %d draws x %d ages x %d kids...\n",
            length(draw_idx), length(AGE_GRID), N_KIDS_PER_AGE))

rows <- vector("list", length(draw_idx) * length(AGE_GRID))
idx <- 0
for (k in draw_idx) {
  sa <- as.numeric(draws$sigma_alpha[k])
  sz <- as.numeric(draws$sigma_zeta[k])
  delta_d <- as.numeric(draws$delta[k])
  s_d  <- as.numeric(draws$s[k])
  sigma_xi <- sqrt(sa^2 + sigma_r^2)

  for (a in AGE_GRID) {
    idx <- idx + 1
    set.seed(SEED + k * 1000 + as.integer(a * 10))
    xi   <- rnorm(N_KIDS_PER_AGE, mu_r, sigma_xi)
    zeta <- rnorm(N_KIDS_PER_AGE, 0, sz)
    kappa <- 1 + delta_d + zeta
    log_age <- log(pmax(a - s_d, 0.01) / a0)
    theta   <- xi + kappa * log_age
    base_kid <- theta + log_H
    eta <- outer(base_kid, delta_j, "-")
    vocab <- rowSums(plogis(eta))
    rows[[idx]] <- data.frame(age = a, vocab = vocab)
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
  pivot_longer(cols = -age, names_to = "tau", values_to = "vocab")

## ---- Empirical quantile fan via GAMLSS (MB-CDI manual approach)
emp_pred <- fit_xsec_quantile_fan(
  emp_df,
  taus = TAUS,
  age_grid = seq(AGE_RANGE[1], AGE_RANGE[2], by = 0.25),
  J_total = J_ITEMS,
  lambda = 10000
)

## ---- Scatter (stratified-by-age, ~80 per bin) ----------------------
set.seed(SEED)
emp_scatter <- emp_df |>
  mutate(age_int = round(age)) |>
  group_by(age_int) |>
  group_modify(~ slice_sample(.x, n = min(nrow(.x), 80))) |>
  ungroup() |>
  select(-age_int)

## ---- Plot ---------------------------------------------------------
p <- ggplot(emp_scatter, aes(age, vocab)) +
  geom_jitter(width = 0.3, alpha = 0.20, size = 0.5, colour = "grey30") +
  geom_line(data = emp_pred,
            aes(age, vocab, group = tau),
            colour = "grey45", linetype = "dashed",
            linewidth = 0.5) +
  geom_line(data = sim_preds,
            aes(age, vocab, group = tau, colour = tau),
            linewidth = 1.1) +
  scale_colour_manual(
    values = c("0.1" = "#1f78b4", "0.25" = "#a6cee3",
               "0.5" = "#33a02c", "0.75" = "#fdbf6f",
               "0.9" = "#e31a1c"),
    name = "Percentile"
  ) +
  scale_x_continuous(breaks = c(16, 20, 24, 28)) +
  coord_cartesian(xlim = c(AGE_RANGE[1] - 0.4, AGE_RANGE[2] + 0.4),
                  ylim = c(0, J_ITEMS)) +
  labs(x = "Age (months)",
       y = sprintf("Productive vocabulary (of %d bundle items)", J_ITEMS),
       title = "M_best at I=500, J=671: alpha + zeta + delta (no s, no s_i)",
       subtitle = sprintf(
         "EN CDI:WS, %d kids (1 admin/kid x-sec). Empirical: grey dots + dashed lines. Predicted: coloured lines.",
         nrow(emp_df))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey25"),
        legend.position = "bottom")

out_png <- file.path(OUT_DIR, "m_best_quantile_I500.png")
ggsave(out_png, p, width = 7, height = 5, dpi = 200)
cat(sprintf("\nWrote: %s\n", out_png))
