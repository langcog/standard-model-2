## Quality-multiplier structural analysis.
##
## Reviewer reflex: "kids' tokens are higher quality than LLM tokens —
## maybe what looks like acceleration (delta) and per-child growth-rate
## variance (sigma_zeta) is actually variation in TOKEN QUALITY across
## kids and across ages."
##
## This script lays out the structural equivalence:
##
## Pure accumulator with per-child, age-varying quality multiplier
## q_i(t) on input:
##
##   eta_ij = (xi_i + log q_{i0}) + (1 + gamma_i) * log_age - psi_j
##
## i.e., adding log q_i(t) = log q_{i0} + gamma_i * log_age to the
## linear predictor is OBSERVATIONALLY EQUIVALENT to the M_best
## parameterization with population delta = mean(gamma_i) and
## per-child zeta_i = gamma_i - mean(gamma_i).
##
## So under "everything is quality" the fit is mathematically identical
## to M_best. Only the labels change. The quality exponents gamma_i
## inherit M_best's (delta, sigma_zeta) posterior:
##
##   gamma_i = delta + zeta_i ~ N(delta, sigma_zeta^2)
##
## We visualize what this implies: how steep does per-child quality
## have to grow with age to absorb what M_best attributes to delta+zeta?
##
## Output:
##   outputs/figs/longitudinal/quality_multiplier_implied.png
##   outputs/figs/longitudinal/quality_multiplier_summary.csv

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Pull M_best posterior --------------------------------------
draws <- readRDS(file.path(PATHS$fits_dir, "summaries",
                           "long_no_freq_slopes.draws.rds"))
sd <- readRDS(file.path(PATHS$fits_dir, "long_subset_data.rds"))$stan_data

as_num <- function(x) as.numeric(unlist(x))
delta_med    <- median(as_num(draws$delta))
sigma_zeta_med <- median(as_num(draws$sigma_zeta))
mu_r <- sd$mu_r
sigma_r <- sd$sigma_r
a0   <- sd$a0
s_med <- median(as_num(draws$s))

cat(sprintf("M_best posterior: delta=%.2f, sigma_zeta=%.2f\n",
            delta_med, sigma_zeta_med))

# ---- Sample per-kid gamma exponents from M_best posterior ----------
N_KIDS_SHOW <- 120
set.seed(2026)
zeta_kid <- rnorm(N_KIDS_SHOW, mean = 0, sd = sigma_zeta_med)
gamma_kid <- delta_med + zeta_kid     # gamma_i for the quality interpretation

# Three reference exponents to highlight
ref_quantiles <- c(0.05, 0.50, 0.95)
ref_gamma <- quantile(gamma_kid, ref_quantiles)

# ---- Build q_i(t) trajectories on the age grid ---------------------
AGE_GRID <- seq(12, 30, by = 0.25)
log_age_rel <- log(pmax(AGE_GRID - s_med, 0.01) / a0)

# q_i(t) = (t/a0)^gamma_i, on the log-quality scale. We center on a0
# so q_i(a0) = 1 for all kids, isolating the AGE-VARYING portion.
log_q_mat <- outer(log_age_rel, gamma_kid)  # length(age) x N_KIDS_SHOW

# Reference trajectories: 5/50/95th percentile gamma
log_q_ref <- outer(log_age_rel, as.numeric(ref_gamma))

df_kid <- expand.grid(age = AGE_GRID, kid = seq_len(N_KIDS_SHOW)) |>
  mutate(log_q = as.vector(log_q_mat))

df_ref <- expand.grid(age = AGE_GRID, pct = names(ref_gamma)) |>
  mutate(label = sprintf("%s percentile (gamma=%.1f)",
                         rep(names(ref_gamma), each = length(AGE_GRID)),
                         rep(as.numeric(ref_gamma), each = length(AGE_GRID))),
         log_q = as.vector(log_q_ref))

# Quantitative summary at age endpoints
endpoints <- AGE_GRID[c(1, length(AGE_GRID))]
endpoint_qs <- sapply(as.numeric(ref_gamma), function(g) {
  log_q_endpoint <- log(pmax(endpoints - s_med, 0.01) / a0) * g
  exp(diff(log_q_endpoint))    # ratio q(30) / q(12)
})

summary_tbl <- tibble::tibble(
  percentile = names(ref_gamma),
  gamma_i    = as.numeric(ref_gamma),
  q_30_over_q_12 = endpoint_qs
)
write.csv(summary_tbl,
          file.path(OUT_DIR, "quality_multiplier_summary.csv"),
          row.names = FALSE)
cat("\n=== Implied quality growth from age 12 to 30 mo ===\n")
print(summary_tbl, digits = 4, row.names = FALSE)

# ---- Plot ----------------------------------------------------------
# Two panels: (a) log_q over age — linear in log_age, fanning;
# (b) implied q ratio q(t)/q(a0) on linear scale.
p_logq <- ggplot() +
  geom_line(data = df_kid, aes(age, log_q, group = kid),
            alpha = 0.10, colour = "grey25", linewidth = 0.4) +
  geom_line(data = df_ref, aes(age, log_q,
                               group = pct, colour = pct),
            linewidth = 1.1) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3,
             linetype = "dashed") +
  scale_colour_manual(
    values = c(`5%` = "#1f78b4", `50%` = "firebrick", `95%` = "#33a02c"),
    name = "kid percentile"
  ) +
  labs(x = "Age (months)",
       y = expression(log * " " * q[i](t) * "  (relative to "*a[0]*")"),
       title = "Implied per-child quality multipliers q_i(t)",
       subtitle = sprintf("Under the structural equivalence: gamma_i = delta + zeta_i ~ N(%.2f, %.2f^2). Quality must grow as (t/a_0)^gamma_i.",
                          delta_med, sigma_zeta_med)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 9, colour = "grey25"))

p_ratio <- ggplot() +
  geom_line(data = df_kid, aes(age, exp(log_q), group = kid),
            alpha = 0.10, colour = "grey25", linewidth = 0.4) +
  geom_line(data = df_ref, aes(age, exp(log_q),
                               group = pct, colour = pct),
            linewidth = 1.1) +
  geom_hline(yintercept = 1, colour = "grey75", linewidth = 0.3,
             linetype = "dashed") +
  scale_y_log10(labels = scales::label_log(),
                breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000)) +
  scale_colour_manual(
    values = c(`5%` = "#1f78b4", `50%` = "firebrick", `95%` = "#33a02c"),
    name = "kid percentile"
  ) +
  annotate("text", x = 30, y = max(exp(df_ref$log_q)),
           hjust = 1, vjust = 0,
           label = sprintf("Median kid: q(30)/q(12) = %.0fx",
                           summary_tbl$q_30_over_q_12[2]),
           size = 3.0, colour = "firebrick") +
  labs(x = "Age (months)",
       y = expression(q[i](t) * " / " * q[i](a[0]) * "  (log scale)"),
       title = "Same trajectories on multiplicative scale",
       subtitle = "If acceleration is 'really' age-varying token quality, this is how much quality must grow per kid.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 9, colour = "grey25"))

composite <- p_logq + p_ratio + plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Quality-multiplier interpretation: same model, different label",
    caption = sprintf(
      "Posterior of M_best gives gamma_i = delta + zeta_i ~ N(%.2f, %.2f^2). For 'token quality grows with age' to absorb M_best's acceleration, median kid needs q(30)/q(12) ~ %.0fx; +/- 1 SD spread covers %.0fx to %.0fx.",
      delta_med, sigma_zeta_med,
      summary_tbl$q_30_over_q_12[2],
      summary_tbl$q_30_over_q_12[1],
      summary_tbl$q_30_over_q_12[3])
  ) & theme(plot.title = element_text(face = "bold"),
            plot.caption = element_text(size = 8, colour = "grey45"))

out_png <- file.path(OUT_DIR, "quality_multiplier_implied.png")
ggsave(out_png, composite, width = 12, height = 5.8, dpi = 150)
cat(sprintf("\nWrote: %s\n", out_png))

cat("\nDone.\n")
