## Parameter-recovery test for the reparameterized log-linear IRT
## accumulator (log_irt_v2.stan).
##
## Note on link: the explainer doc writes the likelihood as probit; this
## implementation uses the logit link (canonical "log-linear IRT"). They
## differ by a scale factor of ~1.7; interpretation of psi as a log-
## threshold is preserved up to that constant. We'll reconcile the doc
## when the fit is working.

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(ggplot2)
})
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

here <- function(...) file.path("/Users/mcfrank/Projects/standard_model_2/model", ...)
source(here("simulate.R"))

## ---- Config -------------------------------------------------------------
N_CHILDREN <- 250
N_WORDS    <- 150
N_CHAINS   <- 4
N_ITER     <- 1500
N_WARMUP   <- 750

## ---- 1. Simulate --------------------------------------------------------
sim <- simulate_data(
  I = N_CHILDREN, J = N_WORDS, C = 3,
  sigma_alpha_true = 0.5,
  mu_c_true  = c(6.5, 8.0, 9.5),
  tau_c_true = c(0.5, 0.7, 0.7),
  s_true = 4.5,
  delta_true = 0.1,
  seed = 42
)
cat(sprintf("Simulated: I=%d children, J=%d words, N=%d obs. Mean y=%.3f\n",
            sim$constants$I, sim$constants$J,
            nrow(sim$obs), mean(sim$obs$y)))

stan_data <- build_stan_data(sim)

## ---- 2. Fit v2 ----------------------------------------------------------
t0 <- Sys.time()
fit <- stan(
  file   = here("log_irt_v2.stan"),
  data   = stan_data,
  chains = N_CHAINS,
  iter   = N_ITER,
  warmup = N_WARMUP,
  seed   = 123,
  control = list(adapt_delta = 0.9, max_treedepth = 10)
)
cat(sprintf("Total sampling time: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## ---- 3. Diagnostics -----------------------------------------------------
cat("\n--- Sampler diagnostics ---\n")
print(check_hmc_diagnostics(fit))

scalar_pars <- c("sigma_alpha", "s", "delta", "pi_alpha", "sigma_xi")
class_pars  <- c(paste0("mu_c[", 1:3, "]"),
                 paste0("tau_c[", 1:3, "]"))
s1 <- summary(fit, pars = c(scalar_pars, class_pars))$summary
print(s1[, c("mean", "2.5%", "50%", "97.5%", "n_eff", "Rhat")])

## ---- 4. Recovery comparison --------------------------------------------
true_pi <- sim$true$sigma_alpha^2 /
           (sim$true$sigma_alpha^2 + sim$true$sigma_r^2)
true_sx <- sqrt(sim$true$sigma_alpha^2 + sim$true$sigma_r^2)

truth_tbl <- data.frame(
  param = c(scalar_pars, class_pars),
  truth = c(sim$true$sigma_alpha, sim$true$s, sim$true$delta,
            true_pi, true_sx,
            sim$true$mu_c, sim$true$tau_c),
  post_median = s1[, "50%"],
  post_lo     = s1[, "2.5%"],
  post_hi     = s1[, "97.5%"]
)
truth_tbl$in_ci <- truth_tbl$truth >= truth_tbl$post_lo &
                   truth_tbl$truth <= truth_tbl$post_hi

cat("\n--- Scalar recovery ---\n")
print(truth_tbl, row.names = FALSE, digits = 3)

## Word-level psi
draws <- as_draws_df(fit)
psi_cols <- grep("^psi\\[", names(draws), value = TRUE)
psi_med <- sapply(psi_cols, function(p) median(draws[[p]]))
psi_lo  <- sapply(psi_cols, function(p) quantile(draws[[p]], .025))
psi_hi  <- sapply(psi_cols, function(p) quantile(draws[[p]], .975))

psi_df <- data.frame(
  j = seq_along(sim$true$psi),
  truth = sim$true$psi,
  post_median = psi_med,
  lo = psi_lo, hi = psi_hi,
  class = factor(sim$true$cc)
)
psi_df$in_ci <- psi_df$truth >= psi_df$lo & psi_df$truth <= psi_df$hi

cat(sprintf("\npsi recovery: correlation(truth, posterior median) = %.3f\n",
            cor(psi_df$truth, psi_df$post_median)))
cat(sprintf("psi CI coverage: %.2f (target ~0.95)\n", mean(psi_df$in_ci)))

## Child-level log_alpha (via posterior-expected)
la_cols <- grep("^log_alpha_mean\\[", names(draws), value = TRUE)
la_med  <- sapply(la_cols, function(p) median(draws[[p]]))
la_df <- data.frame(i = seq_along(sim$true$log_alpha),
                    truth = sim$true$log_alpha,
                    post_median = la_med)
cat(sprintf("log_alpha_mean recovery: correlation(truth, posterior median) = %.3f\n",
            cor(la_df$truth, la_df$post_median)))

## Child-level xi recovery
xi_cols <- grep("^xi\\[", names(draws), value = TRUE)
xi_med  <- sapply(xi_cols, function(p) median(draws[[p]]))
true_xi <- with(sim$true, log_r + log_alpha)
cat(sprintf("xi recovery: correlation(truth, posterior median) = %.3f\n",
            cor(true_xi, xi_med)))

## ---- 5. Plots ----------------------------------------------------------
p_psi <- ggplot(psi_df, aes(truth, post_median, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, alpha = .5) +
  geom_point(size = 1.6) +
  labs(title = "Word-level psi recovery (v2)",
       x = "True psi_j", y = "Posterior median psi_j",
       colour = "Class") +
  theme_minimal(base_size = 12)
ggsave(here("recovery_psi.png"), p_psi, width = 5.5, height = 5, dpi = 150)

p_la <- ggplot(la_df, aes(truth, post_median)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 1.2, alpha = .7) +
  labs(title = "Child-level log_alpha recovery (v2)",
       x = "True log_alpha_i", y = "Posterior expected log_alpha_i") +
  theme_minimal(base_size = 12)
ggsave(here("recovery_log_alpha.png"), p_la, width = 5.5, height = 5, dpi = 150)

saveRDS(fit, here("recovery_fit.rds"))
saveRDS(sim, here("recovery_sim.rds"))

cat("\nDone. Figures: recovery_psi.png, recovery_log_alpha.png\n")
