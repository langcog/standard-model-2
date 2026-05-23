## Input-conditional vocab quantile fan for the IO fits (BabyView,
## SEEDLingS). Visualizes the three-way age × input × vocab relationship
## by grouping kids into input quartiles and showing model-predicted
## trajectories (conditional on input level) overlaid with empirical
## dots coloured by the kid's actual input quartile.
##
## Per-quartile model prediction: fix log_r to the quartile's median,
## sample alpha and zeta from their posteriors, propagate forward to
## get vocab(age) per simulated kid, then take the median across the
## (draws x simulated-kids) pool at each age. Same machinery as the
## marginal quantile-fan plot, except log_r is conditioned rather than
## marginalized.
##
## Output: outputs/figs/longitudinal/m_best_input_quartile_io.png

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

SEED           <- 20260523
N_DRAWS_USE    <- 200
N_KIDS_PER_AGE <- 200
Q_LABELS       <- c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)")

# Cool-to-warm palette for input level. Low input = blue, high = red.
INPUT_PALETTE <- c("Q1 (lowest)" = "#2c7bb6",
                   "Q2"          = "#abd9e9",
                   "Q3"          = "#fdae61",
                   "Q4 (highest)" = "#d7191c")

## ---- Per-bundle panel builder ------------------------------------
build_panel <- function(label,
                        bundle_path,
                        tag,
                        cache_path,
                        title) {
  bundle <- readRDS(bundle_path)
  sd_b   <- bundle$stan_data
  J_ITEMS <- sd_b$J
  age_range <- c(floor(min(sd_b$admin_age)), ceiling(max(sd_b$admin_age)))
  cat(sprintf("[%s] I=%d, J=%d, V=%d, age=[%g, %g]\n",
              label, sd_b$I, J_ITEMS, sd_b$V, age_range[1], age_range[2]))

  ## Per-kid mean log_r from video-level log_r_obs.
  kid_logr <- tibble(
    video      = seq_along(sd_b$log_r_obs),
    log_r_obs  = sd_b$log_r_obs,
    child_idx  = sd_b$video_to_child
  ) |>
    group_by(child_idx) |>
    summarise(log_r = mean(log_r_obs, na.rm = TRUE), .groups = "drop")
  cat(sprintf("[%s] kid log_r: median=%.2f, range=[%.2f, %.2f]\n",
              label, median(kid_logr$log_r), min(kid_logr$log_r), max(kid_logr$log_r)))

  ## Quartile assignment.
  q_breaks <- quantile(kid_logr$log_r, probs = c(0, 0.25, 0.5, 0.75, 1.0))
  q_breaks[1] <- q_breaks[1] - 0.01   # include the minimum
  q_breaks[length(q_breaks)] <- q_breaks[length(q_breaks)] + 0.01
  kid_logr$quartile <- cut(kid_logr$log_r, breaks = q_breaks,
                            labels = Q_LABELS, include.lowest = TRUE)
  # log_r summary per quartile (for the conditional sim)
  q_logr <- kid_logr |>
    group_by(quartile) |>
    summarise(log_r_q = median(log_r), n_kids = n(), .groups = "drop")
  cat(sprintf("[%s] quartile log_r medians:\n", label))
  print(q_logr)

  ## Empirical: bundle's x-sec, joined to kid's quartile.
  emp_df <- build_xsec_empirical(
    bundle_df  = bundle$df,
    age_range  = age_range,
    cache_path = cache_path
  )
  # Map empirical child_id (string) to bundle's child_idx (integer).
  # Bundle's admin_info has ii column matching child_idx; admin_key
  # contains child_id (e.g. "01_6_WG"). Easier: bundle$df has ii column.
  child_idx_map <- bundle$df |>
    distinct(across(intersect(c("subject_id", "child_id", "lab_subject_id"),
                              names(bundle$df))), ii) |>
    rename_with(~ "child_id", -ii)
  emp_df <- emp_df |>
    mutate(child_id = as.character(child_id)) |>
    inner_join(child_idx_map |> mutate(child_id = as.character(child_id)),
               by = "child_id") |>
    inner_join(kid_logr |> select(ii = child_idx, quartile, log_r),
               by = "ii")
  cat(sprintf("[%s] empirical: %d kids w/ quartile assigned\n",
              label, nrow(emp_df)))

  ## Draws.
  draws <- readRDS(file.path(SUMMARIES, paste0(tag, ".draws.rds")))
  if ("draws_df" %in% class(draws)) draws <- as.data.frame(draws)
  if (!"sigma_s" %in% names(draws)) draws$sigma_s <- 0
  if (!"s" %in% names(draws)) draws$s <- 0
  delta_j_df <- read_csv(file.path(SUMMARIES, paste0(tag, "_psi.csv")),
                          show_col_types = FALSE)
  med_col <- if ("delta_j_median" %in% names(delta_j_df)) "delta_j_median" else "psi_median"
  delta_j <- delta_j_df[[med_col]]

  ## Conditional simulation per quartile.
  set.seed(SEED)
  draw_idx <- sort(sample(seq_len(nrow(draws)), N_DRAWS_USE))
  AGE_GRID <- seq(age_range[1], age_range[2], by = 0.5)
  log_H <- sd_b$log_H; a0 <- sd_b$a0

  rows <- vector("list", length(draw_idx) * length(AGE_GRID) * nrow(q_logr))
  idx <- 0
  for (k in draw_idx) {
    sa <- as.numeric(draws$sigma_alpha[k])
    sz <- as.numeric(draws$sigma_zeta[k])
    dd <- as.numeric(draws$delta[k])
    sd_d <- as.numeric(draws$s[k])
    for (qi in seq_len(nrow(q_logr))) {
      log_r_q <- q_logr$log_r_q[qi]
      qlabel  <- as.character(q_logr$quartile[qi])
      for (a in AGE_GRID) {
        idx <- idx + 1
        set.seed(SEED + k * 100000 + qi * 1000 + as.integer(a * 10))
        # Fix log_r at quartile median; vary alpha and zeta from
        # their fitted posteriors. xi = log_r + log_alpha.
        log_alpha <- rnorm(N_KIDS_PER_AGE, 0, sa)
        zeta      <- rnorm(N_KIDS_PER_AGE, 0, sz)
        xi        <- log_r_q + log_alpha
        kappa     <- 1 + dd + zeta
        log_age   <- log(pmax(a - sd_d, 0.01) / a0)
        theta     <- xi + kappa * log_age
        eta       <- outer(theta + log_H, delta_j, "-")
        vocab     <- rowSums(plogis(eta))
        rows[[idx]] <- data.frame(quartile = qlabel,
                                  age = a,
                                  vocab = vocab)
      }
    }
  }
  sim_df <- bind_rows(rows)
  sim_med <- sim_df |>
    group_by(quartile, age) |>
    summarise(med = median(vocab, na.rm = TRUE),
              q25 = quantile(vocab, 0.25, na.rm = TRUE),
              q75 = quantile(vocab, 0.75, na.rm = TRUE),
              .groups = "drop")
  sim_med$quartile <- factor(sim_med$quartile, levels = Q_LABELS)

  ## Plot.
  ggplot() +
    geom_ribbon(data = sim_med,
                aes(age, ymin = q25, ymax = q75, fill = quartile),
                alpha = 0.15, show.legend = FALSE) +
    geom_line(data = sim_med,
              aes(age, med, colour = quartile),
              linewidth = 1.1) +
    geom_jitter(data = emp_df,
                aes(age, vocab, colour = quartile),
                width = 0.2, alpha = 0.70, size = 1.2) +
    scale_colour_manual(values = INPUT_PALETTE, name = "Input quartile") +
    scale_fill_manual(values = INPUT_PALETTE, name = "Input quartile") +
    coord_cartesian(xlim = age_range + c(-0.4, 0.4),
                    ylim = c(0, J_ITEMS)) +
    labs(x = "Age (months)",
         y = sprintf("Productive vocabulary (of %d items)", J_ITEMS),
         title = title,
         subtitle = sprintf(
           "Lines: median trajectory per input quartile (ribbon = IQR over alpha+zeta); dots: kids coloured by their measured log r quartile. I=%d.",
           sd_b$I)) +
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

combined <- (p_bv | p_sd) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Vocabulary trajectory by input quartile",
    subtitle = "IO model (M_best, no s, no s_i). Per-kid log r (mean over recordings) split into quartiles; model trajectory conditional on each quartile's median log r."
  ) &
  theme(legend.position = "bottom")

out_png <- file.path(OUT_DIR, "m_best_input_quartile_io.png")
ggsave(out_png, combined, width = 12, height = 5.5, dpi = 200)
cat(sprintf("\nWrote: %s\n", out_png))
