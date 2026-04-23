## Analysis of Wordbank v2 fit.
## Produces RQ-aligned plots and posterior summaries.

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

here <- function(...) file.path("/Users/mcfrank/Projects/standard_model_2", ...)
TAG <- "subset"   # change to "full" after scaling

fit <- readRDS(here(sprintf("model/wordbank_fit_%s.rds", TAG)))
dat <- readRDS(here(sprintf("model/wordbank_data_%s.rds", TAG)))
draws <- as_draws_df(fit)

word_info <- dat$word_info
class_levels <- dat$class_levels

## ---- RQ4: pi_alpha and variance decomposition -------------------------
cat("\n=== RQ4: Input- vs efficiency-driven variance ===\n")
print(summary(fit, pars = c("sigma_alpha", "sigma_xi", "pi_alpha"))$summary
        [, c("mean", "2.5%", "50%", "97.5%")], digits = 3)

## ---- RQ2 & RQ3: s and delta --------------------------------------------
cat("\n=== RQ2 (start time) & RQ3 (age rate-change) ===\n")
print(summary(fit, pars = c("s", "delta"))$summary
        [, c("mean", "2.5%", "50%", "97.5%")], digits = 3)

## ---- Class-level thresholds -------------------------------------------
cat("\n=== Class thresholds ===\n")
mu_tbl <- summary(fit, pars = paste0("mu_c[", seq_along(class_levels), "]"))$summary
tau_tbl <- summary(fit, pars = paste0("tau_c[", seq_along(class_levels), "]"))$summary
class_tbl <- tibble(
  class = class_levels,
  mu_median = mu_tbl[, "50%"],
  mu_lo = mu_tbl[, "2.5%"], mu_hi = mu_tbl[, "97.5%"],
  tau_median = tau_tbl[, "50%"],
  tau_lo = tau_tbl[, "2.5%"], tau_hi = tau_tbl[, "97.5%"]
)
print(class_tbl, n = Inf)

## ---- RQ1: per-word psi vs log_p ---------------------------------------
cat("\n=== RQ1: how well does frequency reconstruct word difficulty? ===\n")
psi_cols <- grep("^psi\\[", names(draws), value = TRUE)
psi_med  <- sapply(psi_cols, function(p) median(draws[[p]]))
psi_lo   <- sapply(psi_cols, function(p) quantile(draws[[p]], .025))
psi_hi   <- sapply(psi_cols, function(p) quantile(draws[[p]], .975))

psi_df <- word_info %>%
  mutate(psi_median = psi_med,
         psi_lo = psi_lo, psi_hi = psi_hi,
         class = class_levels[cc],
         log_p = log(prob))

## Correlation of psi with log_p (within-sample median)
r_global <- cor(psi_df$psi_median, psi_df$log_p)
r2_global <- r_global^2
cat(sprintf("Global: r(psi, log_p) = %.3f, R^2 = %.3f\n", r_global, r2_global))

## By class
cat("\nBy lexical class:\n")
psi_df %>% group_by(class) %>%
  summarise(n = n(),
            r = cor(psi_median, log_p),
            R2 = r^2) %>%
  print()

## Posterior distribution of correlation (integrate over uncertainty in psi)
n_samples <- min(500, nrow(draws))
psi_draws_mat <- as.matrix(draws[sample(nrow(draws), n_samples), psi_cols])
r_by_sample <- apply(psi_draws_mat, 1, function(v) cor(v, psi_df$log_p))
cat(sprintf("Posterior r(psi, log_p): median=%.3f [%.3f, %.3f]\n",
            median(r_by_sample),
            quantile(r_by_sample, 0.025), quantile(r_by_sample, 0.975)))
cat(sprintf("Posterior R^2: median=%.3f [%.3f, %.3f]\n",
            median(r_by_sample^2),
            quantile(r_by_sample^2, 0.025), quantile(r_by_sample^2, 0.975)))

## ---- Plots -------------------------------------------------------------
p_psi_logp <- ggplot(psi_df, aes(log_p, psi_median, colour = class)) +
  geom_errorbar(aes(ymin = psi_lo, ymax = psi_hi), width = 0, alpha = 0.3) +
  geom_point(size = 1.6, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey40",
              aes(group = 1), linewidth = 0.5, linetype = "dashed") +
  labs(title = "RQ1: word threshold psi vs. log frequency",
       subtitle = sprintf("Global r = %.2f, R^2 = %.2f. Dashed = global OLS.",
                          r_global, r2_global),
       x = "log p_j  (CHILDES)", y = "Posterior median psi_j") +
  theme_minimal(base_size = 12)
ggsave(here("model/wordbank_psi_vs_logp.png"), p_psi_logp,
       width = 7, height = 5, dpi = 150)

## Class mean bars
p_class <- ggplot(class_tbl, aes(reorder(class, mu_median), mu_median)) +
  geom_pointrange(aes(ymin = mu_lo, ymax = mu_hi)) +
  coord_flip() +
  labs(title = "Class-level thresholds (mu_c)",
       x = NULL, y = "Log-threshold (log-tokens)") +
  theme_minimal(base_size = 12)
ggsave(here("model/wordbank_class_means.png"), p_class,
       width = 6, height = 3.5, dpi = 150)

## Posterior density on pi_alpha
p_pi <- ggplot(tibble(pi = draws$pi_alpha), aes(pi)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  geom_vline(xintercept = median(draws$pi_alpha), linetype = "dashed") +
  labs(title = "RQ4: pi_alpha posterior",
       subtitle = "Proportion of child-level outcome variance from learning efficiency",
       x = "pi_alpha", y = "density") +
  xlim(0, 1) +
  theme_minimal(base_size = 12)
ggsave(here("model/wordbank_pi_alpha.png"), p_pi,
       width = 6, height = 4, dpi = 150)

## Posterior density on s and delta
p_s <- ggplot(tibble(s = draws$s), aes(s)) +
  geom_density(fill = "darkorange", alpha = 0.4) +
  labs(title = "RQ2: start-time posterior", x = "s (months)", y = "density") +
  theme_minimal(base_size = 12)
p_d <- ggplot(tibble(d = draws$delta), aes(d)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "RQ3: delta posterior", x = "delta (age rate-change exponent)", y = "density") +
  theme_minimal(base_size = 12)
ggsave(here("model/wordbank_s.png"), p_s, width = 6, height = 4, dpi = 150)
ggsave(here("model/wordbank_delta.png"), p_d, width = 6, height = 4, dpi = 150)

cat("\nFigures saved to model/wordbank_*.png\n")
