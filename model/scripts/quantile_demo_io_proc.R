## Vocab-space quantile-fan plot for the IO + Peekbank fits.
## Three panels: BabyView | SEEDLingS | Stanford-linked (Peekbank).
## Same layout as the EN/NO M_best panels (bundle-internal empirical
## via the bundle's own df + GAMLSS BEINF quantile smoother).
##
## Each panel uses the base no_freq_slopes variant of its dataset.
## For SEEDLingS, the comp/std/comp_std variants add channels that
## constrain log_alpha jointly but produce similar vocab predictions;
## a separate side script can compare them if needed.
##
## Output: outputs/figs/longitudinal/m_best_quantile_io_proc.png

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

TAUS           <- c(0.10, 0.25, 0.50, 0.75, 0.90)
SEED           <- 20260523
N_DRAWS_USE    <- 200
N_KIDS_PER_AGE <- 200
LAMBDA         <- 10000

PALETTE <- c("0.1" = "#1f78b4", "0.25" = "#a6cee3",
             "0.5" = "#33a02c", "0.75" = "#fdbf6f",
             "0.9" = "#e31a1c")

## ---- Per-bundle panel builder ------------------------------------
build_panel <- function(label,
                        bundle_path,
                        tag,
                        cache_path,
                        title,
                        age_range = NULL) {
  bundle <- readRDS(bundle_path)
  sd_b   <- bundle$stan_data
  J_ITEMS <- sd_b$J
  if (is.null(age_range)) {
    age_range <- c(floor(min(sd_b$admin_age)), ceiling(max(sd_b$admin_age)))
  }
  cat(sprintf("[%s] I=%d, J=%d, age=[%g, %g]\n",
              label, sd_b$I, J_ITEMS, age_range[1], age_range[2]))

  emp_df <- build_xsec_empirical(
    bundle_df  = bundle$df,
    age_range  = age_range,
    cache_path = cache_path
  )
  cat(sprintf("[%s] empirical: %d kids, median %d items/kid\n",
              label, nrow(emp_df), median(emp_df$n_items)))

  draws <- readRDS(file.path(SUMMARIES, paste0(tag, ".draws.rds")))
  if ("draws_df" %in% class(draws)) draws <- as.data.frame(draws)
  if (!"sigma_s" %in% names(draws)) draws$sigma_s <- 0
  if (!"s" %in% names(draws)) draws$s <- 0
  delta_j_df <- read_csv(file.path(SUMMARIES, paste0(tag, "_psi.csv")),
                          show_col_types = FALSE)
  med_col <- if ("delta_j_median" %in% names(delta_j_df)) "delta_j_median" else "psi_median"
  delta_j <- delta_j_df[[med_col]]
  if (length(delta_j) != J_ITEMS)
    cat(sprintf("[%s] WARN: delta_j length %d != J %d (may mismatch)\n",
                label, length(delta_j), J_ITEMS))

  set.seed(SEED)
  draw_idx <- sort(sample(seq_len(nrow(draws)), N_DRAWS_USE))
  AGE_GRID <- seq(age_range[1], age_range[2], by = 0.5)
  log_H <- sd_b$log_H; a0 <- sd_b$a0
  # IO variants fit mu_r and sigma_r as parameters (since each kid
  # has observed log_r_obs); proc/EN/NO variants have them pinned in
  # stan_data. Prefer per-draw values if present, else stan_data.
  has_mu_r    <- "mu_r"    %in% names(draws)
  has_sigma_r <- "sigma_r" %in% names(draws)
  mu_r_fixed    <- sd_b$mu_r
  sigma_r_fixed <- sd_b$sigma_r

  rows <- vector("list", length(draw_idx) * length(AGE_GRID))
  idx <- 0
  for (k in draw_idx) {
    sa <- as.numeric(draws$sigma_alpha[k])
    sz <- as.numeric(draws$sigma_zeta[k])
    dd <- as.numeric(draws$delta[k])
    sd_d <- as.numeric(draws$s[k])
    mu_r    <- if (has_mu_r)    as.numeric(draws$mu_r[k])    else mu_r_fixed
    sigma_r <- if (has_sigma_r) as.numeric(draws$sigma_r[k]) else sigma_r_fixed
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

  # Skip the multi-quantile GAMLSS fan -- N is too small (BabyView=20,
  # SEEDLingS=44, Stanford=62) for reliable empirical q10/q90 lines.
  # Just show a single smoothed empirical median (LOESS).
  emp_scatter <- emp_df  # no down-sampling needed at this N

  ggplot(emp_scatter, aes(age, vocab)) +
    geom_jitter(width = 0.3, alpha = 0.40, size = 1.0, colour = "grey30") +
    geom_smooth(method = "loess", se = FALSE, span = 1.0,
                colour = "grey25", linewidth = 0.6,
                linetype = "dashed") +
    geom_line(data = sim_preds,
              aes(age, vocab, group = tau, colour = tau),
              linewidth = 1.1) +
    scale_colour_manual(values = PALETTE, name = "Percentile") +
    coord_cartesian(xlim = age_range + c(-0.4, 0.4),
                    ylim = c(0, J_ITEMS)) +
    labs(x = "Age (months)",
         y = sprintf("Productive vocabulary (of %d items)", J_ITEMS),
         title = title,
         subtitle = sprintf("I=%d, %d kids (1 admin/kid x-sec)",
                            sd_b$I, nrow(emp_df))) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey25"),
          legend.position = "bottom")
}

p_bv <- build_panel(
  label       = "BabyView",
  bundle_path = file.path(PATHS$fits_dir, "babyview_subset_data.rds"),
  tag         = "io_no_freq_slopes",
  cache_path  = file.path(PATHS$fits_dir, "demo_xsec_empirical_babyview.rds"),
  title       = "BabyView (English)"
)
p_sd <- build_panel(
  label       = "SEEDLingS",
  bundle_path = file.path(PATHS$fits_dir, "seedlings_subset_data.rds"),
  tag         = "io_no_freq_slopes_seedlings",
  cache_path  = file.path(PATHS$fits_dir, "demo_xsec_empirical_seedlings.rds"),
  title       = "SEEDLingS (English)"
)
p_sl <- build_panel(
  label       = "Stanford-linked",
  bundle_path = file.path(PATHS$fits_dir, "stanford_linked_subset_data.rds"),
  tag         = "long_proc_no_freq_slopes",
  cache_path  = file.path(PATHS$fits_dir, "demo_xsec_empirical_stanford.rds"),
  title       = "Stanford-linked (Peekbank LWL)"
)

combined <- (p_bv | p_sd | p_sl) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "M_best on IO and processing bundles (vocab-space quantile fans)",
    subtitle = "All variants fit with α + ζ + δ (no s, no s_i). Empirical: bundle's own x-sec scatter + LOESS median (dashed). N too small for stable q10/q90 empirical."
  ) &
  theme(legend.position = "bottom")

out_png <- file.path(OUT_DIR, "m_best_quantile_io_proc.png")
ggsave(out_png, combined, width = 14, height = 5, dpi = 200)
cat(sprintf("\nWrote: %s\n", out_png))
