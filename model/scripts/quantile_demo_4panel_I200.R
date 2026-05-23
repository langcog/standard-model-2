## Four-panel architecture demo at I=200 J=198 (English).
##
## Panels build up the M_best model component-by-component, all fit
## on the same I=200 J=198 N=145K bundle:
##   1. Pure accumulator       - long_demo_pure  (alpha=0, zeta=0, delta=0)
##   2. + per-child efficiency - long_demo_alpha (+ sigma_alpha free)
##   3. + age dynamics         - long_demo_kappa (+ sigma_zeta, delta free)
##   4. M_best (+ both)        - long_no_freq_slopes_english_I200
##
## Empirical: bundle-internal x-sec (1 admin/kid) via build_xsec_empirical,
## GAMLSS BEINF lambda=10000 for the quantile-fan smoother. Same machinery
## as the EN/NO side-by-side plots.
##
## Output: outputs/figs/longitudinal/quantile_demo_4panel_I200.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/empirical_xsec_helper.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(patchwork)
})

OUT_DIR    <- file.path(PATHS$figs_dir, "longitudinal")
SUMMARIES  <- file.path(PATHS$fits_dir, "summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BUNDLE_PATH    <- file.path(PATHS$fits_dir, "long_subset_data_I200.rds")
CACHE_PATH     <- file.path(PATHS$fits_dir, "demo_xsec_empirical_I200.rds")
AGE_RANGE      <- c(16, 30)
TAUS           <- c(0.10, 0.25, 0.50, 0.75, 0.90)
SEED           <- 20260523
N_DRAWS_USE    <- 200
N_KIDS_PER_AGE <- 200
LAMBDA         <- 10000

PALETTE <- c("0.1" = "#1f78b4", "0.25" = "#a6cee3",
             "0.5" = "#33a02c", "0.75" = "#fdbf6f",
             "0.9" = "#e31a1c")

VARIANTS <- list(
  pure   = list(tag = "long_demo_pure",
                title = "1. Pure accumulator"),
  alpha  = list(tag = "long_demo_alpha",
                title = "2. + child efficiency (alpha)"),
  kappa  = list(tag = "long_demo_kappa",
                title = "3. + age dynamics (kappa)"),
  mbest  = list(tag = "long_no_freq_slopes_english_I200",
                title = "4. M_best (alpha + zeta + delta)")
)

bundle  <- readRDS(BUNDLE_PATH)
sd_b    <- bundle$stan_data
J_ITEMS <- sd_b$J
cat(sprintf("I=200 bundle: I=%d, J=%d, N=%d\n",
            sd_b$I, J_ITEMS, sd_b$N))

## Build empirical (cross-sectional, one admin per kid).
emp_df <- build_xsec_empirical(
  bundle_df  = bundle$df,
  age_range  = AGE_RANGE,
  cache_path = CACHE_PATH
)
cat(sprintf("Empirical (x-sec): %d kids, median %d items/kid\n",
            nrow(emp_df), median(emp_df$n_items)))

emp_pred <- fit_xsec_quantile_fan(
  emp_df, taus = TAUS,
  age_grid = seq(AGE_RANGE[1], AGE_RANGE[2], by = 0.25),
  J_total = J_ITEMS,
  lambda = LAMBDA
)

set.seed(SEED)
emp_scatter <- emp_df |>
  mutate(age_int = round(age)) |>
  group_by(age_int) |>
  group_modify(~ slice_sample(.x, n = min(nrow(.x), 80))) |>
  ungroup() |>
  select(-age_int)

## Per-variant panel builder.
build_panel <- function(variant_key) {
  v <- VARIANTS[[variant_key]]
  draws <- readRDS(file.path(SUMMARIES, paste0(v$tag, ".draws.rds")))
  if ("draws_df" %in% class(draws)) draws <- as.data.frame(draws)
  for (p in c("sigma_alpha", "sigma_zeta", "delta", "s", "sigma_s")) {
    if (!p %in% names(draws)) draws[[p]] <- 0
  }
  delta_j_df <- read_csv(file.path(SUMMARIES, paste0(v$tag, "_psi.csv")),
                          show_col_types = FALSE)
  med_col <- if ("delta_j_median" %in% names(delta_j_df)) "delta_j_median" else "psi_median"
  delta_j <- delta_j_df[[med_col]]
  J_panel <- length(delta_j)
  cat(sprintf("[%s] tag=%s, %d draws, J=%d (bundle J=%d)\n",
              v$title, v$tag, nrow(draws), J_panel, J_ITEMS))

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

  ggplot(emp_scatter, aes(age, vocab)) +
    geom_jitter(width = 0.3, alpha = 0.20, size = 0.5, colour = "grey30") +
    geom_line(data = emp_pred,
              aes(age, vocab, group = tau),
              colour = "grey45", linetype = "dashed",
              linewidth = 0.5) +
    geom_line(data = sim_preds,
              aes(age, vocab, group = tau, colour = tau),
              linewidth = 1.0) +
    scale_colour_manual(values = PALETTE, name = "Percentile") +
    scale_x_continuous(breaks = c(16, 20, 24, 28)) +
    coord_cartesian(xlim = c(AGE_RANGE[1] - 0.4, AGE_RANGE[2] + 0.4),
                    ylim = c(0, J_ITEMS)) +
    labs(x = "Age (months)",
         y = sprintf("Productive vocabulary (of %d items)", J_ITEMS),
         title = v$title) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "bottom")
}

panels <- lapply(names(VARIANTS), build_panel)
combined <- (panels[[1]] | panels[[2]]) /
            (panels[[3]] | panels[[4]]) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Building up M_best (English I=200, J=198)",
    subtitle = "Each panel = one component added. Empirical: grey dots (1 admin/kid x-sec) + dashed lines (GAMLSS BEINF lambda=10000)."
  ) &
  theme(legend.position = "bottom")

out_png <- file.path(OUT_DIR, "quantile_demo_4panel_I200.png")
ggsave(out_png, combined, width = 11, height = 9, dpi = 200)
cat(sprintf("\nWrote: %s\n", out_png))
