## Reduced-form test: does input RATE predict the per-child SLOPE?
##
## The accumulator puts input rate (log r_i) in the intercept only: more
## input is a level shift, not a change in the growth exponent. The
## "bootstrapping / Matthew effect" alternative says more input also
## makes you accelerate faster — i.e., log r_i should predict the
## per-child slope deviation ζ_i (γ > 0).
##
## This is only identifiable where input is MEASURED. The IO fits
## (BabyView, SEEDLingS) anchor per-child log_r_true to observed token
## counts, and also estimate per-child ζ. So we can ask, directly:
##   does log_r_true covary with ζ across kids?
##
## We do it across posterior draws (not just medians) so the answer
## carries honest uncertainty:
##   for each draw d:  γ_d = slope of lm(ζ[d,] ~ z(log_r_true[d,]))
##                      r_d = cor(ζ[d,], log_r_true[d,])
## then summarize γ_d, r_d over draws.
##
## γ in standardized units: SD of ζ per +1 SD of log_r_true.
##
## Output:
##   outputs/input_slope_check.csv
##   outputs/figs/longitudinal/input_slope_check.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(rstan); library(dplyr); library(tidyr); library(ggplot2)
  library(readr)
})

OUT_CSV <- file.path(PATHS$outputs_dir, "input_slope_check.csv")
OUT_PNG <- file.path(PATHS$figs_dir, "longitudinal", "input_slope_check.png")

FITS <- list(
  BabyView  = "fits/io_no_freq_slopes.rds",
  SEEDLingS = "fits/io_no_freq_slopes_seedlings.rds"
)

analyze_one <- function(label, path) {
  cat(sprintf("\n=== %s (%s) ===\n", label, path))
  fit <- readRDS(path)
  lrt  <- rstan::extract(fit, "log_r_true")$log_r_true   # [draws, I]
  zeta <- rstan::extract(fit, "zeta")$zeta               # [draws, I]
  xi   <- rstan::extract(fit, "xi")$xi                   # [draws, I]
  I <- ncol(lrt); nd <- nrow(lrt)
  cat(sprintf("  I=%d kids, %d draws\n", I, nd))

  # per-draw reduced-form slope (standardized predictor) + correlation
  gamma_d <- numeric(nd); cor_d <- numeric(nd)
  for (d in seq_len(nd)) {
    x <- as.numeric(scale(lrt[d, ]))      # z-scored input rate
    y <- zeta[d, ]
    gamma_d[d] <- coef(lm(y ~ x))[2]      # ζ units per +1 SD log_r_true
    cor_d[d]   <- cor(lrt[d, ], zeta[d, ])
  }

  qfun <- function(v) quantile(v, c(0.5, 0.025, 0.975))
  g <- qfun(gamma_d); r <- qfun(cor_d)
  # also: correlation of input rate with the INTERCEPT xi, as a sanity
  # check (should be strongly positive — input lives in the intercept)
  cor_xi_d <- sapply(seq_len(nd), function(d) cor(lrt[d, ], xi[d, ]))
  cxi <- qfun(cor_xi_d)

  cat(sprintf("  gamma (ζ SD-change per +1 SD log_r_true): %+.3f [%+.3f, %+.3f]\n",
              g[1], g[2], g[3]))
  cat(sprintf("  cor(log_r_true, zeta):  %+.3f [%+.3f, %+.3f]\n",
              r[1], r[2], r[3]))
  cat(sprintf("  cor(log_r_true, xi)  :  %+.3f [%+.3f, %+.3f]  (sanity: should be high)\n",
              cxi[1], cxi[2], cxi[3]))
  p_pos <- mean(gamma_d > 0)
  cat(sprintf("  P(gamma > 0) = %.2f\n", p_pos))

  # posterior-median per-child values for the scatter
  scatter <- tibble(
    dataset      = label,
    child        = seq_len(I),
    log_r_true   = apply(lrt,  2, median),
    zeta         = apply(zeta, 2, median),
    xi           = apply(xi,   2, median)
  )

  summary_row <- tibble(
    dataset = label, I = I,
    gamma_med = g[1], gamma_lo = g[2], gamma_hi = g[3],
    cor_zeta_med = r[1], cor_zeta_lo = r[2], cor_zeta_hi = r[3],
    cor_xi_med = cxi[1], cor_xi_lo = cxi[2], cor_xi_hi = cxi[3],
    p_gamma_pos = p_pos
  )
  list(summary = summary_row, scatter = scatter)
}

res <- lapply(names(FITS), function(nm) analyze_one(nm, FITS[[nm]]))
summ    <- bind_rows(lapply(res, `[[`, "summary"))
scatter <- bind_rows(lapply(res, `[[`, "scatter"))

write_csv(summ, OUT_CSV)
cat(sprintf("\nWrote %s\n", OUT_CSV))
cat("\n=== summary ===\n"); print(summ)

## ---- Plot: per-child zeta vs log_r_true, per dataset ----------------
lab_df <- summ |>
  mutate(label = sprintf("%s (N=%d)\nγ = %+.2f [%+.2f, %+.2f]",
                          dataset, I, gamma_med, gamma_lo, gamma_hi))
scatter <- scatter |> left_join(lab_df |> select(dataset, label), by = "dataset")

p <- ggplot(scatter, aes(log_r_true, zeta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(alpha = 0.6, size = 1.6, colour = "#1f78b4") +
  geom_smooth(method = "lm", se = TRUE, colour = "#d7191c",
              linewidth = 0.8, fill = "#d7191c", alpha = 0.15) +
  facet_wrap(~ label, scales = "free_x") +
  labs(x = "Per-child input rate  log r_i  (posterior median)",
       y = "Per-child slope deviation  ζ_i  (posterior median)",
       title = "Does input rate predict the per-child growth slope?",
       subtitle = "Reduced-form test of input-on-slope (γ). Accumulator predicts γ = 0 (input is a level effect). γ summarized over posterior draws; band = lm SE on medians.") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        strip.text = element_text(face = "bold"))

ggsave(OUT_PNG, p, width = 9, height = 4.5, dpi = 200)
cat(sprintf("Wrote %s\n", OUT_PNG))
