## Parameter recovery test: simulate data with known parameters, fit, and
## verify the posteriors cover truth.
##
## Usage:   Rscript model/scripts/01_recovery.R
## Outputs: fits/recovery.rds, outputs/figs/recovery_*.png
##          logs scalar recovery table to stdout.

source("model/R/config.R")
source("model/R/helpers.R")

FIT_CFG <- modifyList(DEFAULT_FIT_CONFIG,
                      list(iter = 1500, warmup = 750,
                           adapt_delta = 0.9))

# Simulate with known true values
sim <- simulate_data(
  I = 250, J = 150, C = 3,
  sigma_alpha_true = 0.5,
  mu_c_true  = c(6.5, 8.0, 9.5),
  tau_c_true = c(0.5, 0.7, 0.7),
  s_true = 4.5, delta_true = 0.1,
  seed = 42
)
stan_data <- build_stan_data_from_sim(sim)
cat(sprintf("Simulated: I=%d, J=%d, N=%d (mean y=%.3f)\n",
            sim$constants$I, sim$constants$J,
            nrow(sim$obs), mean(sim$obs$y)))

fit <- fit_variant(stan_data, "recovery", cfg = FIT_CFG)

cat("\n--- Sampler diagnostics ---\n")
print(check_hmc_diagnostics(fit))

cat("\n--- Scalar recovery ---\n")
scalar_tbl <- summarize_fit(fit)
true_vals <- c(sigma_alpha = sim$true$sigma_alpha,
               s           = sim$true$s,
               delta       = sim$true$delta,
               pi_alpha    = sim$true$sigma_alpha^2 /
                             (sim$true$sigma_alpha^2 + sim$true$sigma_r^2),
               sigma_xi    = sqrt(sim$true$sigma_alpha^2 + sim$true$sigma_r^2))
scalar_tbl$truth  <- true_vals[scalar_tbl$param]
scalar_tbl$in_ci  <- scalar_tbl$truth >= scalar_tbl$lo95 &
                    scalar_tbl$truth <= scalar_tbl$hi95
print(scalar_tbl[, c("param", "truth", "median", "lo95", "hi95", "in_ci",
                     "n_eff", "Rhat")], digits = 3)

# Per-word psi recovery
draws <- as_draws_df(fit)
psi_cols <- grep("^psi\\[", names(draws), value = TRUE)
psi_med <- sapply(psi_cols, function(p) median(draws[[p]]))
psi_lo  <- sapply(psi_cols, function(p) quantile(draws[[p]], .025))
psi_hi  <- sapply(psi_cols, function(p) quantile(draws[[p]], .975))
psi_df <- tibble(j = seq_along(sim$true$psi),
                 truth = sim$true$psi,
                 post_median = psi_med,
                 lo = psi_lo, hi = psi_hi,
                 class = factor(sim$true$cc))
psi_df$in_ci <- psi_df$truth >= psi_df$lo & psi_df$truth <= psi_df$hi
cat(sprintf("\npsi recovery: r(truth, posterior median) = %.3f\n",
            cor(psi_df$truth, psi_df$post_median)))
cat(sprintf("psi CI coverage: %.2f (target 0.95)\n", mean(psi_df$in_ci)))

p_psi <- ggplot(psi_df, aes(truth, post_median, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, alpha = .5) +
  geom_point(size = 1.6) +
  labs(title = "Recovery: word thresholds",
       x = "True psi_j", y = "Posterior median psi_j", colour = "Class") +
  theme_minimal(base_size = 12)
ggsave(file.path(PATHS$figs_dir, "recovery_psi.png"),
       p_psi, width = 5.5, height = 5, dpi = 150)

cat("\nDone. Fit: fits/recovery.rds. Figures: outputs/figs/recovery_*.png\n")
