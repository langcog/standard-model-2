## How much quality variation (sigma_q) would we need to push pi_alpha
## below various thresholds?
##
## The setup. The model has xi_i = log r_i + log alpha_i with
##   log r_i  ~ N(mu_r, sigma_r)        # input rate (measured)
##   log alpha_i ~ N(0, sigma_alpha)    # efficiency (latent)
## sigma_xi^2 = sigma_r^2 + sigma_alpha^2  (data-identified)
## pi_alpha = sigma_alpha^2 / sigma_xi^2
##
## Now add a constant-across-development per-child quality
## multiplier q_i:
##   xi_i = log r_i + log q_i + log alpha'_i
##   log q_i ~ N(0, sigma_q)      # quality (latent; unmeasured)
##   log alpha'_i ~ N(0, sigma_alpha')
## sigma_xi^2 unchanged = sigma_r^2 + sigma_q^2 + sigma_alpha'^2
## pi_alpha_TRUE = sigma_alpha'^2 / sigma_xi^2
##
## The model currently absorbs sigma_q into its sigma_alpha estimate
## (no way to separate q from alpha without external quality data).
## So:
##   sigma_alpha^2 (what we report) = sigma_q^2 + sigma_alpha'^2
##   pi_alpha (what we report) >= pi_alpha_TRUE (we OVER-attribute to
##     efficiency when quality variation exists).
##
## This script computes: for the EN M_best posterior (sigma_alpha = 1.81,
## sigma_r = 0.534), how big does sigma_q need to be to push the TRUE
## pi_alpha below a given threshold?
##
## Note: this analysis applies ONLY to constant-per-child quality.
## Age-varying quality q_i(t) = (t/a_0)^gamma_i is a different beast --
## it competes with kappa_i, not with sigma_alpha. See the appendix
## slide for that analysis.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(knitr); library(ggplot2)
})

OUT_DIR <- PATHS$outputs_dir

# Headline EN M_best posteriors (from long_no_freq_slopes.summary.rds).
SIGMA_ALPHA <- 1.81        # what the model reports
SIGMA_R     <- 0.534       # external pin (Sperry within-pool)
SIGMA_XI    <- sqrt(SIGMA_R^2 + SIGMA_ALPHA^2)  # ~1.887
PI_ALPHA_NOW <- SIGMA_ALPHA^2 / SIGMA_XI^2      # ~0.92

cat(sprintf(
"Current EN M_best:
  sigma_alpha = %.3f, sigma_r = %.3f
  sigma_xi    = sqrt(sigma_r^2 + sigma_alpha^2) = %.3f
  pi_alpha    = sigma_alpha^2 / sigma_xi^2     = %.3f\n\n",
  SIGMA_ALPHA, SIGMA_R, SIGMA_XI, PI_ALPHA_NOW))

## ---- Required sigma_q for each target pi_alpha_TRUE -----------
## For target T = pi_alpha_TRUE:
##   sigma_alpha'^2 = T * sigma_xi^2
##   sigma_q^2      = sigma_alpha^2 - sigma_alpha'^2
##                  = sigma_alpha^2 - T * sigma_xi^2
##   sigma_q        = sqrt(...) -- valid when nonnegative.
required_sigma_q <- function(target_pi) {
  s2_a_prime <- target_pi * SIGMA_XI^2
  s2_q       <- SIGMA_ALPHA^2 - s2_a_prime
  ifelse(s2_q < 0, NA, sqrt(s2_q))
}

targets <- c(0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20)
tab <- tibble(
  pi_alpha_target = targets,
  sigma_q_needed  = required_sigma_q(targets)
) |>
  mutate(
    multiplicative_factor_per_SD = exp(sigma_q_needed),
    spread_16_to_84              = exp(2 * sigma_q_needed),
    sigma_q_vs_sigma_r           = sigma_q_needed / SIGMA_R
  )

cat("=== Quality SD required to push pi_alpha below each threshold ===\n\n")
print(knitr::kable(tab, format = "pipe", digits = 2,
                    align = c("r","r","r","r","r")))
cat(sprintf(
"\nNotes:
  * sigma_q is on the log scale (same units as sigma_r, sigma_alpha).
  * 'multiplicative factor per SD' = exp(sigma_q): how much more
    'effective input' a kid at +1 SD of quality experiences vs the
    median kid, per token of measured input.
  * 'spread 16-to-84' = exp(2*sigma_q): the full ±1 SD range as a
    multiplicative factor.
  * 'sigma_q vs sigma_r' is the ratio to the within-pool measured
    input-rate variation (0.534). >1 means quality varies more than
    quantity does.
  * The current pi_alpha (~0.92) is what the model reports if all of
    the kid-level intercept variance beyond input rate is called
    'efficiency'. The 'true' pi_alpha is lower whenever any of that
    variance is actually constant per-child quality.\n"))

write_csv(tab, file.path(OUT_DIR, "quality_variation_table.csv"))
writeLines(capture.output(print(knitr::kable(tab, format = "pipe",
                                              digits = 2,
                                              align = c("r","r","r","r","r")))),
           file.path(OUT_DIR, "quality_variation_table.md"))
cat(sprintf("\nWrote %s\n", file.path(OUT_DIR, "quality_variation_table.csv")))

## ---- Plot: pi_alpha as a function of sigma_q ------------------
sweep <- tibble(sigma_q = seq(0, 1.5, by = 0.01)) |>
  mutate(
    sigma_alpha_prime = sqrt(pmax(SIGMA_ALPHA^2 - sigma_q^2, 0)),
    pi_alpha_true     = sigma_alpha_prime^2 / SIGMA_XI^2,
    pi_alpha_reported = PI_ALPHA_NOW
  )

label_pts <- tibble(
  sigma_q = required_sigma_q(c(0.90, 0.70, 0.50)),
  pi_alpha = c(0.90, 0.70, 0.50)
) |> filter(!is.na(sigma_q))

p <- ggplot(sweep, aes(sigma_q, pi_alpha_true)) +
  geom_hline(yintercept = PI_ALPHA_NOW, linetype = "dotted",
             colour = "grey40") +
  geom_vline(xintercept = SIGMA_R, linetype = "dotted",
             colour = "#1f78b4") +
  geom_line(linewidth = 1.0, colour = "#d7191c") +
  geom_point(data = label_pts,
             aes(x = sigma_q, y = pi_alpha),
             inherit.aes = FALSE,
             colour = "black", size = 2.5) +
  geom_text(data = label_pts,
             aes(x = sigma_q, y = pi_alpha,
                 label = sprintf("sigma_q=%.2f -> pi_alpha=%.2f",
                                 sigma_q, pi_alpha)),
             inherit.aes = FALSE,
             hjust = -0.05, vjust = -0.4, size = 3.2) +
  annotate("text", x = SIGMA_R + 0.02, y = 0.05,
           label = sprintf("sigma_r=%.3f\n(measured)", SIGMA_R),
           hjust = 0, size = 3.0, colour = "#1f78b4") +
  annotate("text", x = 0.05, y = PI_ALPHA_NOW + 0.025,
           label = sprintf("current pi_alpha = %.2f\n(if sigma_q = 0)",
                            PI_ALPHA_NOW),
           hjust = 0, size = 3.0, colour = "grey40") +
  scale_x_continuous(name = "Hypothetical quality SD on log scale (sigma_q)",
                     limits = c(0, 1.5), breaks = seq(0, 1.5, 0.25)) +
  scale_y_continuous(name = "True pi_alpha (= efficiency share of intercept variance)",
                     limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "How much constant-per-child quality variation would erode pi_alpha?",
       subtitle = sprintf(
         "EN M_best: sigma_alpha=%.2f, sigma_r=%.3f, sigma_xi=%.3f. As sigma_q grows, variance previously assigned to efficiency reattributes to quality.",
         SIGMA_ALPHA, SIGMA_R, SIGMA_XI)) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"))

out_png <- file.path(PATHS$figs_dir, "longitudinal",
                      "quality_variation_pi_alpha.png")
ggsave(out_png, p, width = 8, height = 5, dpi = 200)
cat(sprintf("Wrote %s\n", out_png))
