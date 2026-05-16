## Theta-space spaghetti plot: per-kid latent ability trajectories
## across the 5-panel model build.
##
## Each panel shows N_KIDS simulated kid trajectories of theta_i(t)
## generated from the population distribution implied by that model's
## posterior medians. Theta is the kid-level latent ability in logit
## units; sigmoid(theta - delta_j) is the per-word acquisition
## probability.
##
## The five models progressively turn on:
##   1. pure              -- no variation, no acceleration
##   2. + alpha           -- per-kid efficiency variation (vertical spread)
##   3. + kappa_pop       -- population acceleration (steeper slope)
##   4. + sigma_zeta      -- per-kid acceleration (fan opens with age)
##   5. + sigma_s         -- per-kid trajectory phase (asymmetric fan)
##
## The point: show what the model SAYS about kid trajectories, not
## what vocab fan it predicts. Cleaner reading of structural diffs.
##
## Output: outputs/figs/longitudinal/theta_spaghetti.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

AGE_GRID <- seq(8, 36, by = 0.5)
N_KIDS <- 100
SEED <- 20260516
mu_r    <- 7.337734    # external prior; same as bundle
sigma_r <- 0.5338644
log_H   <- log(365)
a0      <- 19          # English median admin age
s_pop   <- 0.5         # current pinning; will become 6 after refit

## ---- Posterior values per model -----------------------------------
##
## For models 1-4 we use the current English M_best (long_no_freq_slopes)
## fit's posteriors where applicable, with the relevant SDs zeroed out
## as the variant specifies. For model 5 we use long_no_freq_slopes_si_signed.
## After the s=6 refit lands, these values will update -- the script
## reads them from summaries when available, else falls back to these
## defaults.
read_or_default <- function(tag, defaults) {
  path <- file.path("fits/summaries", paste0(tag, ".summary.rds"))
  if (!file.exists(path)) return(defaults)
  s <- readRDS(path)
  med <- setNames(s$median, s$variable)
  vals <- defaults
  for (k in names(vals)) {
    if (k %in% names(med)) vals[[k]] <- as.numeric(med[[k]])
  }
  vals
}

# Defaults reflect current best-known posteriors at s=0.5; will be
# automatically overridden once new fits are extracted.
v_best <- read_or_default("long_no_freq_slopes",
                          list(sigma_alpha = 1.71, sigma_zeta = 3.47,
                               delta = 9.39, sigma_s = 0))
v_signed <- read_or_default("long_no_freq_slopes_si_signed",
                            list(sigma_alpha = 1.56, sigma_zeta = 3.51,
                                 delta = 9.62, sigma_s = 1.40))

MODELS <- list(
  list(label = "1. Pure accumulator",
       sigma_alpha = 0, sigma_zeta = 0, sigma_s = 0,
       delta = 0),
  list(label = "2. + α (efficiency)",
       sigma_alpha = v_best$sigma_alpha,
       sigma_zeta = 0, sigma_s = 0,
       delta = 0),
  list(label = "3. + κ_pop (pop. acceleration)",
       sigma_alpha = v_best$sigma_alpha,
       sigma_zeta = 0, sigma_s = 0,
       delta = v_best$delta),
  list(label = "4. + σ_ζ (per-kid acceleration)",
       sigma_alpha = v_best$sigma_alpha,
       sigma_zeta = v_best$sigma_zeta,
       sigma_s = 0,
       delta = v_best$delta),
  list(label = "5. + σ_s (per-kid phase)",
       sigma_alpha = v_signed$sigma_alpha,
       sigma_zeta = v_signed$sigma_zeta,
       sigma_s = v_signed$sigma_s,
       delta = v_signed$delta)
)

## ---- Simulate theta_i(t) per panel --------------------------------
simulate_panel <- function(m) {
  set.seed(SEED)
  sigma_xi <- sqrt(sigma_r^2 + m$sigma_alpha^2)
  xi   <- rnorm(N_KIDS, mu_r, sigma_xi)
  zeta <- rnorm(N_KIDS, 0, m$sigma_zeta)
  s_i  <- rnorm(N_KIDS, 0, m$sigma_s)
  s_i  <- s_i - mean(s_i)  # sum-to-zero centered (matches the model)
  kappa <- 1 + m$delta + zeta

  # theta_i(t) = xi_i + log_H + kappa_i * log((t - s_pop - s_i) / a_0)
  out <- expand.grid(kid = seq_len(N_KIDS), age = AGE_GRID,
                     KEEP.OUT.ATTRS = FALSE)
  out$theta <- xi[out$kid] + log_H +
               kappa[out$kid] * log(pmax(out$age - s_pop - s_i[out$kid], 0.01) / a0)
  out$panel <- m$label
  out
}

cat("Simulating panels...\n")
panel_data <- bind_rows(lapply(MODELS, simulate_panel))

# Population mean trajectory per panel (bold line)
mean_trajs <- panel_data |>
  group_by(panel, age) |>
  summarise(theta_mean = mean(theta), .groups = "drop")

# Quantile bands per panel
quantile_bands <- panel_data |>
  group_by(panel, age) |>
  summarise(q10 = quantile(theta, 0.10),
            q50 = quantile(theta, 0.50),
            q90 = quantile(theta, 0.90),
            .groups = "drop")

## ---- Plot ----------------------------------------------------------
panel_levels <- sapply(MODELS, function(m) m$label)
panel_data$panel    <- factor(panel_data$panel, levels = panel_levels)
mean_trajs$panel    <- factor(mean_trajs$panel, levels = panel_levels)
quantile_bands$panel <- factor(quantile_bands$panel, levels = panel_levels)

# Y-axis reference line at theta = median(delta_j). From English data,
# median per-word difficulty is around 15 logits. At theta = 15, a typical
# word has 50% acquisition probability. Add as annotation.
THETA_HALF <- 15

p_main <- ggplot(panel_data, aes(age, theta, group = kid)) +
  # Reference line: theta where typical-difficulty word has 50% prob.
  geom_hline(yintercept = THETA_HALF, colour = "#666666", linewidth = 0.5,
             linetype = "dotted") +
  geom_line(alpha = 0.10, colour = "grey25", linewidth = 0.25) +
  geom_line(data = quantile_bands, aes(age, q10, group = NULL),
            colour = "#1f78b4", linewidth = 0.8, linetype = "dashed") +
  geom_line(data = quantile_bands, aes(age, q90, group = NULL),
            colour = "#e31a1c", linewidth = 0.8, linetype = "dashed") +
  geom_line(data = mean_trajs, aes(age, theta_mean, group = NULL),
            colour = "black", linewidth = 1.0) +
  # Annotate the reference line on panel 1 only.
  annotate("text", x = 9, y = THETA_HALF + 0.6,
           label = "θ_50% (typical word)",
           hjust = 0, vjust = 0, colour = "#666666", size = 2.9, alpha = 0.85,
           data = data.frame(panel = factor(panel_levels[1],
                                            levels = panel_levels))) +
  facet_wrap(~ panel, ncol = 5) +
  scale_x_continuous(breaks = c(8, 16, 24, 32)) +
  coord_cartesian(xlim = c(8, 36), ylim = c(5, 25)) +
  labs(x = "Age (months)",
       y = expression("Latent ability " ~ theta[i](t) ~ " (logit units)"),
       title = "Predicted latent ability across the 5-model build",
       subtitle = paste0(
         "Each panel: 100 simulated kids (grey), population mean (black), ",
         "10th/90th percentiles (blue/red dashed). Dotted: θ where a ",
         "typical-difficulty word has 50% acquisition probability. ",
         "Posterior values from current English fits (s=0.5); will refresh after s=6 refit."
       )) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey25"),
        panel.spacing.x = unit(1.2, "lines"))

out_png <- file.path(OUT_DIR, "theta_spaghetti.png")
ggsave(out_png, p_main, width = 16, height = 4.5, dpi = 200)
cat(sprintf("\nWrote: %s\n", out_png))

# Print the posterior values used
cat("\n=== Posterior values used per panel ===\n")
for (m in MODELS) {
  cat(sprintf("%s:\n", m$label))
  cat(sprintf("  sigma_alpha = %.2f, sigma_zeta = %.2f, sigma_s = %.2f, delta = %.2f\n",
              m$sigma_alpha, m$sigma_zeta, m$sigma_s, m$delta))
}
