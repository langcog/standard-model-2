## Four-panel theta-space architecture demo at I=200. Companion to
## model/scripts/quantile_demo_4panel_I200.R: same architecture
## build-up but plotted in latent-ability theta(t) space rather than
## vocab-count space.
##
## Posterior values: read from the I=200 fits:
##   panel 1 (pure):  no posterior needed, all sigmas = 0
##   panel 2 (+ alpha): long_demo_alpha (sigma_alpha free)
##   panel 3 (+ kappa_pop + sigma_zeta): long_demo_kappa
##   panel 4 (M_best): long_no_freq_slopes_english_I200
##
## All panels simulate N_KIDS trajectories from the variant's posterior
## medians; theta(t) = xi + log_H + kappa * log((t - s)/a0) per kid.
## Identical math to the quantile-fan plot, just stops at theta.
##
## Output: outputs/figs/longitudinal/theta_spaghetti_4panel_I200.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BUNDLE  <- readRDS(file.path(PATHS$fits_dir, "long_subset_data_I200.rds"))
sd_b    <- BUNDLE$stan_data
mu_r    <- sd_b$mu_r
sigma_r <- sd_b$sigma_r
log_H   <- sd_b$log_H
a0      <- sd_b$a0
cat(sprintf("Bundle: mu_r=%.2f sigma_r=%.2f log_H=%.2f a0=%g\n",
            mu_r, sigma_r, log_H, a0))

AGE_GRID <- seq(12, 32, by = 0.5)   # I=200 era spans roughly this range
N_KIDS   <- 100
SEED     <- 20260523

# Posterior medians per variant. For variants where a parameter wasn't
# free in that fit, the fitted value is ~0 (or near it).
read_meds <- function(tag) {
  path <- file.path(PATHS$fits_dir, "summaries", paste0(tag, ".summary.rds"))
  if (!file.exists(path)) stop("Missing summary: ", path)
  s <- readRDS(path)
  setNames(s$median, s$variable)
}

m_alpha <- read_meds("long_demo_alpha")
m_kappa <- read_meds("long_demo_kappa")
m_best  <- read_meds("long_no_freq_slopes_english_I200")

cat("Posterior medians:\n")
cat(sprintf("  alpha:  sigma_alpha=%.2f\n", m_alpha[["sigma_alpha"]]))
cat(sprintf("  kappa:  sigma_alpha=%.2f delta=%.2f\n",
            m_kappa[["sigma_alpha"]], m_kappa[["delta"]]))
cat(sprintf("  m_best: sigma_alpha=%.2f sigma_zeta=%.2f delta=%.2f\n",
            m_best[["sigma_alpha"]], m_best[["sigma_zeta"]], m_best[["delta"]]))

MODELS <- list(
  list(label = "1. Pure accumulator",
       sigma_alpha = 0, sigma_zeta = 0, delta = 0),
  list(label = "2. + alpha (efficiency)",
       sigma_alpha = m_alpha[["sigma_alpha"]],
       sigma_zeta = 0, delta = 0),
  list(label = "3. + kappa (age dynamics)",
       sigma_alpha = m_kappa[["sigma_alpha"]],
       sigma_zeta = m_kappa[["sigma_zeta"]],
       delta = m_kappa[["delta"]]),
  list(label = "4. M_best (alpha + zeta + delta)",
       sigma_alpha = m_best[["sigma_alpha"]],
       sigma_zeta = m_best[["sigma_zeta"]],
       delta = m_best[["delta"]])
)

## ---- Simulate theta_i(t) per panel ----------------------------------
simulate_panel <- function(m) {
  set.seed(SEED)
  sigma_xi <- sqrt(sigma_r^2 + m$sigma_alpha^2)
  xi   <- rnorm(N_KIDS, mu_r, sigma_xi)
  zeta <- rnorm(N_KIDS, 0, m$sigma_zeta)
  kappa <- 1 + m$delta + zeta

  out <- expand.grid(kid = seq_len(N_KIDS), age = AGE_GRID,
                     KEEP.OUT.ATTRS = FALSE)
  out$theta <- xi[out$kid] + log_H +
               kappa[out$kid] * log(pmax(out$age, 0.01) / a0)
  out$panel <- m$label
  out
}

panel_data <- bind_rows(lapply(MODELS, simulate_panel))
mean_trajs <- panel_data |>
  group_by(panel, age) |>
  summarise(theta_mean = mean(theta), .groups = "drop")
quantile_bands <- panel_data |>
  group_by(panel, age) |>
  summarise(q10 = quantile(theta, 0.10),
            q50 = quantile(theta, 0.50),
            q90 = quantile(theta, 0.90),
            .groups = "drop")

panel_levels <- sapply(MODELS, function(m) m$label)
panel_data$panel    <- factor(panel_data$panel, levels = panel_levels)
mean_trajs$panel    <- factor(mean_trajs$panel, levels = panel_levels)
quantile_bands$panel <- factor(quantile_bands$panel, levels = panel_levels)

THETA_HALF <- 15   # typical delta_j; theta = 15 -> 50% on a typical word.

p_main <- ggplot(panel_data, aes(age, theta, group = kid)) +
  geom_hline(yintercept = THETA_HALF, colour = "#666666", linewidth = 0.5,
             linetype = "dotted") +
  geom_line(alpha = 0.10, colour = "grey25", linewidth = 0.25) +
  geom_line(data = quantile_bands, aes(age, q10, group = NULL),
            colour = "#1f78b4", linewidth = 0.8, linetype = "dashed") +
  geom_line(data = quantile_bands, aes(age, q90, group = NULL),
            colour = "#e31a1c", linewidth = 0.8, linetype = "dashed") +
  geom_line(data = mean_trajs, aes(age, theta_mean, group = NULL),
            colour = "black", linewidth = 1.0) +
  annotate("text", x = 12.5, y = THETA_HALF + 0.6,
           label = "theta_50% (typical word)",
           hjust = 0, vjust = 0, colour = "#666666", size = 2.9, alpha = 0.85,
           data = data.frame(panel = factor(panel_levels[1],
                                            levels = panel_levels))) +
  facet_wrap(~ panel, ncol = 2) +
  scale_x_continuous(breaks = c(12, 16, 20, 24, 28, 32)) +
  coord_cartesian(xlim = c(12, 32), ylim = c(0, 25)) +
  labs(x = "Age (months)",
       y = expression("Latent ability " ~ theta[i](t) ~ " (logit units)"),
       title = "Building up M_best (I=200, theta-space view)",
       subtitle = paste0(
         "100 simulated kids/panel. Grey: per-kid trajectories. Black: pop mean. ",
         "Blue/red dashed: 10/90th percentiles. Dotted: theta where typical word reaches 50%."
       )) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        panel.spacing = unit(1.2, "lines"))

out_png <- file.path(OUT_DIR, "theta_spaghetti_4panel_I200.png")
ggsave(out_png, p_main, width = 11, height = 9, dpi = 200)
cat(sprintf("\nWrote: %s\n", out_png))
