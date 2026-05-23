## Processing-speed-conditional vocab quantile plot for the Stanford-
## linked / Peekbank fit. Analogous to model/scripts/quantile_demo_io_input.R
## but stratifies kids by their LWL reaction-time intercept rather than
## by their observed input rate.
##
## Per-kid processing speed score: mean(log_rt_t - mu_rtslope * log_age_t)
## across their LWL trials -- "what's this kid's RT once you take out
## the population age trend?" Lower = faster.
##
## Model conditional trajectory: invert the proc-Stan link
##   log_rt_intercept_i = mu_rt - gamma_rt * log_alpha_i
## per posterior draw to recover log_alpha_q from each quartile's median
## intercept, then propagate forward with xi ~ N(mu_r + log_alpha_q,
## sigma_r), zeta ~ N(0, sigma_zeta). Bundle's mu_r / sigma_r are pinned
## (no IO channel here), so log_r variation enters as the sigma_r
## width of xi around the conditional mean.
##
## Output: outputs/figs/longitudinal/m_best_rt_quartile_proc.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/empirical_xsec_helper.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
})

OUT_DIR    <- file.path(PATHS$figs_dir, "longitudinal")
SUMMARIES  <- file.path(PATHS$fits_dir, "summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TAG_PROC       <- "long_proc_no_freq_slopes"
SEED           <- 20260523
N_DRAWS_USE    <- 200
N_KIDS_PER_AGE <- 200
Q_LABELS       <- c("Q1 (fastest)", "Q2", "Q3", "Q4 (slowest)")

# Fast = blue, slow = red. (Note ordering reversed vs. input panels;
# in input, "high" is red because more input = more vocab; in RT,
# "low log_rt" = faster processing = more vocab.)
RT_PALETTE <- c("Q1 (fastest)" = "#2c7bb6",
                "Q2"           = "#abd9e9",
                "Q3"           = "#fdae61",
                "Q4 (slowest)" = "#d7191c")

bundle <- readRDS(file.path(PATHS$fits_dir, "stanford_linked_subset_data.rds"))
sd_b   <- bundle$stan_data
J_ITEMS  <- sd_b$J
age_range <- c(floor(min(sd_b$admin_age)), ceiling(max(sd_b$admin_age)))
cat(sprintf("Proc bundle: I=%d, J=%d, N_lwl=%d, age=[%g, %g]\n",
            sd_b$I, J_ITEMS, sd_b$N_lwl, age_range[1], age_range[2]))

## ---- Per-kid RT intercept (age-adjusted) ------------------------
# We need a single number per kid summarizing their RT performance.
# Subtract the population age trend (mu_rtslope * log_age) from each
# trial, then average per kid -- that approximates the per-kid RT
# intercept the model is identifying.
draws <- readRDS(file.path(SUMMARIES, paste0(TAG_PROC, ".draws_full.rds")))
draws <- as.data.frame(draws)

mu_rtslope_med <- median(draws$mu_rtslope)
mu_rt_med      <- median(draws$mu_rt)
gamma_rt_med   <- median(draws$gamma_rt)
cat(sprintf("Posterior medians: mu_rt=%.3f, mu_rtslope=%.3f, gamma_rt=%.3f\n",
            mu_rt_med, mu_rtslope_med, gamma_rt_med))

trials <- tibble(
  child_idx = sd_b$lwl_to_child,
  log_age   = sd_b$lwl_log_age,
  log_rt    = sd_b$lwl_log_rt
)
kid_rt <- trials |>
  mutate(rt_intercept = log_rt - mu_rtslope_med * log_age) |>
  group_by(child_idx) |>
  summarise(rt_intercept_kid = mean(rt_intercept),
            n_trials = n(),
            .groups = "drop")
cat(sprintf("Per-kid RT intercept: median=%.2f, range=[%.2f, %.2f] (across %d kids)\n",
            median(kid_rt$rt_intercept_kid), min(kid_rt$rt_intercept_kid),
            max(kid_rt$rt_intercept_kid), nrow(kid_rt)))

q_breaks <- quantile(kid_rt$rt_intercept_kid, probs = c(0, 0.25, 0.5, 0.75, 1.0))
q_breaks[1] <- q_breaks[1] - 0.01
q_breaks[length(q_breaks)] <- q_breaks[length(q_breaks)] + 0.01
kid_rt$quartile <- cut(kid_rt$rt_intercept_kid, breaks = q_breaks,
                        labels = Q_LABELS, include.lowest = TRUE)
q_rt <- kid_rt |>
  group_by(quartile) |>
  summarise(rt_q   = median(rt_intercept_kid),
            n_kids = n(),
            .groups = "drop")
cat("Quartile RT intercept medians:\n"); print(q_rt)

## ---- Empirical: bundle x-sec joined to RT quartile --------------
cache_path <- file.path(PATHS$fits_dir, "demo_xsec_empirical_stanford.rds")
emp_df <- build_xsec_empirical(
  bundle_df  = bundle$df,
  age_range  = age_range,
  cache_path = cache_path
)
child_idx_map <- bundle$df |>
  distinct(across(intersect(c("lab_subject_id","subject_id","child_id"),
                            names(bundle$df))), ii) |>
  rename_with(~ "child_id", -ii)
emp_df <- emp_df |>
  mutate(child_id = as.character(child_id)) |>
  inner_join(child_idx_map |> mutate(child_id = as.character(child_id)),
             by = "child_id") |>
  inner_join(kid_rt |> select(ii = child_idx, quartile, rt_intercept_kid),
             by = "ii")
cat(sprintf("Empirical kids w/ RT quartile assigned: %d\n", nrow(emp_df)))

## ---- Load delta_j and other scalars ----------------------------
delta_j_df <- read_csv(file.path(SUMMARIES, paste0(TAG_PROC, "_psi.csv")),
                        show_col_types = FALSE)
med_col <- if ("delta_j_median" %in% names(delta_j_df)) "delta_j_median" else "psi_median"
delta_j <- delta_j_df[[med_col]]
if (!"sigma_s" %in% names(draws)) draws$sigma_s <- 0
if (!"s" %in% names(draws)) draws$s <- 0

## ---- Conditional simulation per RT quartile --------------------
set.seed(SEED)
draw_idx <- sort(sample(seq_len(nrow(draws)), N_DRAWS_USE))
AGE_GRID <- seq(age_range[1], age_range[2], by = 0.5)
log_H   <- sd_b$log_H
a0      <- sd_b$a0
mu_r_fx <- sd_b$mu_r
sigma_r_fx <- sd_b$sigma_r

rows <- vector("list", length(draw_idx) * length(AGE_GRID) * nrow(q_rt))
idx <- 0
for (k in draw_idx) {
  sz <- as.numeric(draws$sigma_zeta[k])
  dd <- as.numeric(draws$delta[k])
  sd_d <- as.numeric(draws$s[k])
  gam_k <- as.numeric(draws$gamma_rt[k])
  mu_rt_k <- as.numeric(draws$mu_rt[k])
  for (qi in seq_len(nrow(q_rt))) {
    rt_q   <- q_rt$rt_q[qi]
    qlabel <- as.character(q_rt$quartile[qi])
    # Invert Stan link to recover this quartile's implied log_alpha
    # under THIS posterior draw. log_rt_intercept = mu_rt - gamma_rt*log_alpha
    # => log_alpha = (mu_rt - rt_q) / gamma_rt.
    log_alpha_q <- (mu_rt_k - rt_q) / gam_k
    for (a in AGE_GRID) {
      idx <- idx + 1
      set.seed(SEED + k * 100000 + qi * 1000 + as.integer(a * 10))
      # xi for simulated kids in this quartile: log_r_dev varies
      # (sigma_r), log_alpha is fixed at the quartile's implied value.
      xi <- rnorm(N_KIDS_PER_AGE, mean = mu_r_fx + log_alpha_q,
                   sd = sigma_r_fx)
      zeta <- rnorm(N_KIDS_PER_AGE, 0, sz)
      kappa   <- 1 + dd + zeta
      log_age <- log(pmax(a - sd_d, 0.01) / a0)
      theta   <- xi + kappa * log_age
      eta     <- outer(theta + log_H, delta_j, "-")
      vocab   <- rowSums(plogis(eta))
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

## ---- Plot ------------------------------------------------------
p <- ggplot() +
  geom_ribbon(data = sim_med,
              aes(age, ymin = q25, ymax = q75, fill = quartile),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(data = sim_med,
            aes(age, med, colour = quartile),
            linewidth = 1.1) +
  geom_jitter(data = emp_df,
              aes(age, vocab, colour = quartile),
              width = 0.2, alpha = 0.75, size = 1.4) +
  scale_colour_manual(values = RT_PALETTE, name = "RT quartile") +
  scale_fill_manual(values = RT_PALETTE, name = "RT quartile") +
  coord_cartesian(xlim = age_range + c(-0.4, 0.4),
                  ylim = c(0, J_ITEMS)) +
  labs(x = "Age (months)",
       y = sprintf("Productive vocabulary (of %d items)", J_ITEMS),
       title = "Vocabulary trajectory by processing-speed quartile",
       subtitle = sprintf(
         "Stanford-linked (Peekbank LWL). Lines: median per-quartile trajectory (ribbon = IQR over log_r_dev + zeta). Dots: kids by RT quartile. I=%d, N_lwl=%d.",
         sd_b$I, sd_b$N_lwl)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        legend.position = "bottom")

out_png <- file.path(OUT_DIR, "m_best_rt_quartile_proc.png")
ggsave(out_png, p, width = 7.5, height = 5.5, dpi = 200)
cat(sprintf("\nWrote: %s\n", out_png))
