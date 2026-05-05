## Three figures for the "M_best explainer" slide insert.
##
## Outputs (under outputs/figs/longitudinal/):
##   m_best_spaghetti.png        — per-child trajectory fan, ability + vocab panels
##   m_best_variance_decomp.png  — between-child variance composition by age
##
## A coefficient table goes directly in the LaTeX slide; no figure for it.
##
## Pulls posterior from a chosen fit (default `long_no_freq_slopes`,
## English M_best). Pass a different fit_name to retarget.
##
## Run:
##   Rscript model/scripts/m_best_explainer_plots.R                  # English M_best
##   Rscript model/scripts/m_best_explainer_plots.R long_no_freq_slopes_norwegian
##
## Per-item difficulties (ψ_j) are pulled from the M1 ψ_j extract on
## the same item set (`outputs/figs/longitudinal/psi_freq_regression_per_word.csv`)
## as a stand-in for the spaghetti's vocab-count axis. M1 and M_best
## differ in absolute ψ_j scale (M_best fits to a different time term)
## but the within-fit ordering and spread are similar; for an
## explanatory plot this is honest. Caption notes the substitution.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
FIT_NAME <- if (length(args) >= 1) args[1] else "long_no_freq_slopes"

# Stan-data path map (mirrors scaling_disanalogy_plot.R)
stan_data_paths <- list(
  long_no_freq_slopes           = file.path(PATHS$fits_dir, "long_subset_data.rds"),
  long_slopes                   = file.path(PATHS$fits_dir, "long_subset_data.rds"),
  long_no_freq_slopes_norwegian = file.path(PATHS$fits_dir, "long_subset_data_nor.rds"),
  long_slopes_norwegian         = file.path(PATHS$fits_dir, "long_subset_data_nor.rds")
)

# ---- Pull posterior + bundle constants -----------------------------
draws <- readRDS(file.path(PATHS$fits_dir, "summaries",
                           paste0(FIT_NAME, ".draws.rds")))
sd_path <- stan_data_paths[[FIT_NAME]]
if (is.null(sd_path) || !file.exists(sd_path)) {
  sd_path <- file.path(PATHS$fits_dir, "long_subset_data.rds")
}
sd <- readRDS(sd_path)$stan_data

as_num <- function(x) as.numeric(unlist(x))
P <- list(
  fit_name = FIT_NAME,
  delta       = as_num(draws$delta),
  sigma_alpha = as_num(draws$sigma_alpha),
  sigma_xi    = as_num(draws$sigma_xi),
  sigma_zeta  = as_num(draws$sigma_zeta),
  rho_xi_zeta = as_num(draws$rho_xi_zeta),
  s           = as_num(draws$s),
  pi_alpha    = as_num(draws$pi_alpha),
  mu_r        = sd$mu_r,
  sigma_r     = sd$sigma_r,
  log_H       = sd$log_H,
  a0          = sd$a0
)
P$delta_med   <- median(P$delta)
P$sxi_med     <- median(P$sigma_xi)
P$sa_med      <- median(P$sigma_alpha)
P$sz_med      <- median(P$sigma_zeta)
P$rho_med     <- median(P$rho_xi_zeta)
P$s_med       <- median(P$s)
P$pia_med     <- median(P$pi_alpha)

cat(sprintf("Fit '%s': δ=%.2f σ_α=%.2f σ_ξ=%.2f σ_ζ=%.2f ρ=%.2f π_α=%.2f, a0=%.0f\n",
            FIT_NAME, P$delta_med, P$sa_med, P$sxi_med, P$sz_med,
            P$rho_med, P$pia_med, P$a0))

# Per-item ψ_j (M1 stand-in for vocab calibration)
psi_csv <- file.path(OUT_DIR, "psi_freq_regression_per_word.csv")
if (!file.exists(psi_csv) && grepl("norwegian", FIT_NAME)) {
  warning("Norwegian ψ_j extract not on disk; using English ψ_j as proxy for vocab calibration.")
}
psi_df <- read.csv(psi_csv, stringsAsFactors = FALSE)
psi_j  <- psi_df$psi_median
class_j <- psi_df$lexical_class
N_ITEMS <- length(psi_j)

# ---- Figure A: spaghetti trajectories on ability + vocab axes -------
##
## Sample N kids from the population MVN (ξ, ζ); compute θ(t) per kid;
## convert θ → expected vocab by summing sigmoid(θ - ψ_j) across CDI items.
N_KIDS_SHOW <- 120
# Restrict age range to where the longitudinal CDI sample is dense.
# Earlier ages have severe ζ-driven fan-out that's a model extrapolation
# artifact, not data; for an explainer slide we use 12-30.
AGE_GRID    <- seq(12, 30, by = 0.25)
set.seed(2026)

Sigma <- matrix(c(P$sxi_med^2,                    P$rho_med * P$sxi_med * P$sz_med,
                  P$rho_med * P$sxi_med * P$sz_med, P$sz_med^2), 2, 2)
L <- chol(Sigma)
z <- matrix(rnorm(2 * N_KIDS_SHOW), 2, N_KIDS_SHOW)
effs <- t(L) %*% z
xi_kid   <- P$mu_r + effs[1, ]
zeta_kid <- effs[2, ]

age_eff <- pmax(AGE_GRID - P$s_med, 0.01)
log_age_rel <- log(age_eff / P$a0)

theta_mat <- outer(rep(1, length(AGE_GRID)), xi_kid) +
             outer(log_age_rel, 1 + P$delta_med + zeta_kid)
# Subtract mu_r so y-axis is "ability above population log-rate intercept"
theta_centered <- theta_mat - P$mu_r

# Vocab: for each (age, kid), sum sigmoid(theta - psi_j) across items.
# Vectorize: at each age, vocab[i] = sum_j sigma(theta_mat[a, i] - psi_j)
vocab_mat <- matrix(NA_real_, nrow = length(AGE_GRID), ncol = N_KIDS_SHOW)
for (a in seq_along(AGE_GRID)) {
  # theta_mat[a, ] is length N_KIDS_SHOW; we want sum over items of sigmoid
  # broadcast: outer(theta_mat[a, ], psi_j, FUN = "-") gives N_KIDS x N_ITEMS
  diff <- outer(theta_mat[a, ], psi_j, FUN = "-")
  vocab_mat[a, ] <- rowSums(plogis(diff))
}

# Population mean and percentiles per age
pct_band <- function(M, q) apply(M, 1, quantile, probs = q, na.rm = TRUE)
theta_med <- pct_band(theta_centered, 0.5)
vocab_med <- pct_band(vocab_mat, 0.5)
vocab_p25 <- pct_band(vocab_mat, 0.25); vocab_p75 <- pct_band(vocab_mat, 0.75)
vocab_p05 <- pct_band(vocab_mat, 0.05); vocab_p95 <- pct_band(vocab_mat, 0.95)

# Long-form for ggplot
df_kid <- expand.grid(age = AGE_GRID, kid = seq_len(N_KIDS_SHOW)) |>
  mutate(theta = as.vector(theta_centered),
         vocab = as.vector(vocab_mat))

df_pop <- tibble::tibble(age = AGE_GRID,
                         theta_med = theta_med,
                         vocab_med = vocab_med,
                         vocab_p25 = vocab_p25, vocab_p75 = vocab_p75,
                         vocab_p05 = vocab_p05, vocab_p95 = vocab_p95)

p_theta <- ggplot() +
  geom_line(data = df_kid, aes(age, theta, group = kid),
            alpha = 0.10, colour = "grey25", linewidth = 0.4) +
  geom_line(data = df_pop, aes(age, theta_med),
            colour = "firebrick", linewidth = 1.3) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  labs(x = "Age (months)",
       y = expression("Ability  " * theta[i*","*t] - mu[r] * "  (logits)"),
       title = "Per-child ability trajectories",
       subtitle = sprintf("120 kids sampled from posterior population MVN; red = population median.\nM_best (%s): delta=%.2f, sigma_xi=%.2f, sigma_zeta=%.2f, rho(xi,zeta)=%.2f.",
                          FIT_NAME, P$delta_med, P$sxi_med, P$sz_med, P$rho_med)) +
  coord_cartesian(xlim = c(12, 30), ylim = c(-10, 10)) +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(size = 9, colour = "grey25"))

p_vocab <- ggplot() +
  geom_ribbon(data = df_pop, aes(age, ymin = vocab_p05, ymax = vocab_p95),
              alpha = 0.18, fill = "steelblue") +
  geom_ribbon(data = df_pop, aes(age, ymin = vocab_p25, ymax = vocab_p75),
              alpha = 0.30, fill = "steelblue") +
  geom_line(data = df_kid, aes(age, vocab, group = kid),
            alpha = 0.10, colour = "grey25", linewidth = 0.4) +
  geom_line(data = df_pop, aes(age, vocab_med),
            colour = "navy", linewidth = 1.3) +
  labs(x = "Age (months)",
       y = sprintf("Predicted productive vocab (out of %d items)", N_ITEMS),
       title = "Same kids, vocabulary scale",
       subtitle = sprintf("Vocab = sum_j sigmoid(theta_it - psi_j) across the %d CDI items in the longitudinal sample.\nBlue ribbons: 25/75 and 5/95 percentiles across kids.", N_ITEMS)) +
  coord_cartesian(xlim = c(12, 30), ylim = c(0, N_ITEMS)) +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(size = 9, colour = "grey25"))

composite <- p_theta + p_vocab + plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = sprintf("M_best (%s): per-child trajectory fan", FIT_NAME),
    caption = "Item difficulties psi_j for vocab calibration: posterior medians from the M1 (time-only) fit, used here as a stand-in (M_best psi_j not separately extracted; absolute scale would shift slightly, ranking and spread match)."
  ) & theme(plot.title = element_text(face = "bold"),
            plot.caption = element_text(size = 7.5, colour = "grey45"))

out_png <- file.path(OUT_DIR, "m_best_spaghetti.png")
ggsave(out_png, composite, width = 12, height = 5.5, dpi = 150)
cat("Wrote:", out_png, "\n")

# ---- Figure B: between-child variance composition by age ------------
##
## Var_i(θ_{it}) = σ_ξ² + log_age² · σ_ζ² + 2ρ σ_ξ σ_ζ log_age,
## with σ_ξ² = σ_α² + σ_r².
## At a₀ (log_age=0), Var = σ_ξ² and the input/efficiency share is exact.
## Off-pivot, slope variance contributes σ_ζ²·log_age² and the cross
## term 2ρσ_ξσ_ζ·log_age (can be negative). We show the four
## contributions stacked (with the cross term split into ± parts).

age_panel <- c(12, 16, 19, 24, 30)  # show variance composition at these ages
sr <- P$sigma_r
sa <- P$sa_med
sxi <- P$sxi_med
sz <- P$sz_med
rho <- P$rho_med

decomp <- lapply(age_panel, function(t) {
  log_age <- log(max(t - P$s_med, 0.01) / P$a0)
  v_input  <- sr^2
  v_eff    <- sa^2
  v_slope  <- (sz * log_age)^2
  v_cross  <- 2 * rho * sxi * sz * log_age
  total    <- v_input + v_eff + v_slope + v_cross
  tibble::tibble(
    age = t, log_age = log_age,
    component = c("Input (sigma_r^2)", "Efficiency (sigma_alpha^2)",
                  "Slope (sigma_zeta^2 * log_age^2)",
                  "Cross (2 rho * sigma_xi * sigma_zeta * log_age)"),
    value = c(v_input, v_eff, v_slope, v_cross),
    total = total
  )
}) |> bind_rows()

decomp$component <- factor(decomp$component,
                           levels = c("Input (sigma_r^2)",
                                      "Efficiency (sigma_alpha^2)",
                                      "Slope (sigma_zeta^2 * log_age^2)",
                                      "Cross (2 rho * sigma_xi * sigma_zeta * log_age)"))

# For stacked bar, separate ± values: positives stack up, negatives stack down
decomp$sign <- ifelse(decomp$value >= 0, "pos", "neg")
totals_df <- decomp |>
  group_by(age) |>
  summarise(
    total_var   = sum(value),
    pos_top     = sum(pmax(value, 0)),
    .groups = "drop"
  )

p_decomp <- ggplot(decomp, aes(factor(age), value, fill = component)) +
  geom_col(position = "stack", colour = "white", linewidth = 0.3) +
  geom_text(data = totals_df,
            aes(factor(age), pos_top,
                label = sprintf("Var(theta) = %.1f", total_var)),
            inherit.aes = FALSE,
            vjust = -0.6, size = 3.2, fontface = "bold", colour = "grey25") +
  scale_fill_manual(values = c(
    "Input (sigma_r^2)"                              = "#a6cee3",
    "Efficiency (sigma_alpha^2)"                     = "#1f78b4",
    "Slope (sigma_zeta^2 * log_age^2)"               = "#fb9a99",
    "Cross (2 rho * sigma_xi * sigma_zeta * log_age)" = "#fdbf6f"
  )) +
  labs(x = "Age (months)",
       y = expression("Variance contribution to  " * Var[i]( theta[it] )),
       fill = NULL,
       title = "Between-child variance composition by age (M_best)",
       subtitle = sprintf("sigma_r=%.2f (input), sigma_alpha=%.2f (efficiency), sigma_zeta=%.2f (growth-rate spread), rho(xi,zeta)=%.2f.\nAt a_0=%.0f mo, log_age=0 -> slope and cross terms vanish; pi_alpha = sigma_alpha^2 / (sigma_alpha^2 + sigma_r^2) = %.2f exact.",
                          sr, sa, sz, rho, P$a0, P$pia_med)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 9, colour = "grey25"))

out_png2 <- file.path(OUT_DIR, "m_best_variance_decomp.png")
ggsave(out_png2, p_decomp, width = 10, height = 5.5, dpi = 150)
cat("Wrote:", out_png2, "\n")

# ---- Save tabular data for the LaTeX coefficient table -------------
##
## Posterior medians + 95% CrI for the headline scalars, plus units and
## plain-English interpretation. The slide will format this as a table.
coef_tbl <- tibble::tibble(
  param = c("\\mu_r", "\\sigma_r", "\\sigma_\\alpha", "\\sigma_\\xi",
            "\\pi_\\alpha", "\\delta", "\\sigma_\\zeta", "\\rho(\\xi,\\zeta)",
            "s"),
  posterior = c(
    sprintf("%.2f (pinned)", P$mu_r),
    sprintf("%.3f (pinned)", P$sigma_r),
    sprintf("%.2f [%.2f, %.2f]", P$sa_med, quantile(P$sigma_alpha, 0.025), quantile(P$sigma_alpha, 0.975)),
    sprintf("%.2f [%.2f, %.2f]", P$sxi_med, quantile(P$sigma_xi, 0.025), quantile(P$sigma_xi, 0.975)),
    sprintf("%.2f [%.2f, %.2f]", P$pia_med, quantile(P$pi_alpha, 0.025), quantile(P$pi_alpha, 0.975)),
    sprintf("%.2f [%.2f, %.2f]", P$delta_med, quantile(P$delta, 0.025), quantile(P$delta, 0.975)),
    sprintf("%.2f [%.2f, %.2f]", P$sz_med, quantile(P$sigma_zeta, 0.025), quantile(P$sigma_zeta, 0.975)),
    sprintf("%+.2f [%+.2f, %+.2f]", P$rho_med, quantile(P$rho_xi_zeta, 0.025), quantile(P$rho_xi_zeta, 0.975)),
    sprintf("%.2f [%.2f, %.2f]", P$s_med, quantile(P$s, 0.025), quantile(P$s, 0.975))
  ),
  units = c("log tokens/hr", "log scale", "logits", "logits",
            "proportion", "dimensionless", "logits / log_age", "correlation",
            "months"),
  meaning = c(
    sprintf("Geometric mean input rate (~%.0f tok/hr).", exp(P$mu_r)),
    "Between-child input variation (Sperry/HR/WF).",
    "Between-child efficiency variation.",
    "Total intercept variance: σ_ξ² = σ_α² + σ_r².",
    "Share of intercept variance from efficiency, not input.",
    "Population age-acceleration: ability ∝ T^(1+δ).",
    "Between-child variation in growth rate around δ.",
    "Sign and strength of intercept-slope coupling.",
    "Global accumulation start age (pinned near 0)."
  )
)
write.csv(coef_tbl, file.path(OUT_DIR, "m_best_coef_table.csv"), row.names = FALSE)
cat("Wrote:", file.path(OUT_DIR, "m_best_coef_table.csv"), "\n")

cat("\nDone.\n")
