## Analyze io_slopes / babyview fit. Three plots:
##  1. Per-child observed input distributions vs inferred log_r_true
##  2. Per-recording log_r_obs scatter (with reactivity inflation)
##  3. Per-child predicted vs observed vocab trajectory + scalar summaries
##
## Output: outputs/figs/seedlings_io_*.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(posterior); library(rstan)
})

fit_path <- file.path(PATHS$fits_dir, "io_slopes_seedlings.rds")
fit    <- readRDS(fit_path)
bundle <- load_dataset_bundle("seedlings")
sd_    <- bundle$stan_data
draws  <- as_draws_df(fit)

OUT_FIGS <- file.path(PATHS$figs_dir, "io")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("SEEDLingS io_slopes: I=%d, A=%d, V=%d, J=%d, N=%d\n",
            sd_$I, sd_$A, sd_$V, sd_$J, sd_$N))

# ---- 1. Scalar posteriors panel ---- #
scalars <- c("sigma_alpha", "pi_alpha", "sigma_xi" = NA,
             "sigma_zeta", "delta", "s",
             "beta_react", "sigma_within", "sigma_r")
sc_summary <- tibble(
  param  = c("sigma_alpha", "pi_alpha", "sigma_zeta", "delta",
             "beta_react", "sigma_within", "sigma_r"),
  median = NA_real_, lo95 = NA_real_, hi95 = NA_real_
)
for (i in seq_len(nrow(sc_summary))) {
  p <- sc_summary$param[i]
  if (p %in% names(draws)) {
    sc_summary$median[i] <- median(draws[[p]])
    sc_summary$lo95[i]   <- quantile(draws[[p]], 0.025)
    sc_summary$hi95[i]   <- quantile(draws[[p]], 0.975)
  }
}
sc_summary <- sc_summary %>%
  mutate(param = factor(param, levels = param))

p_scalars <- ggplot(sc_summary,
                    aes(y = param, x = median, xmin = lo95, xmax = hi95)) +
  geom_pointrange(color = "steelblue", size = 0.4) +
  labs(x = "posterior median, 95% CrI", y = NULL,
       title = "SEEDLingS io_slopes: scalar posteriors",
       subtitle = sprintf("I=%d kids, V=%d videos, J=%d items, N=%d obs",
                          sd_$I, sd_$V, sd_$J, sd_$N)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(family = "mono", size = 10))

# ---- 2. Per-recording observed log_r_obs vs per-child inferred ---- #
videos <- bundle$recordings %>% as_tibble()
log_r_true_med <- sapply(seq_len(sd_$I),
                          function(i) median(draws[[sprintf("log_r_true[%d]", i)]]))
log_r_true_lo  <- sapply(seq_len(sd_$I),
                          function(i) quantile(draws[[sprintf("log_r_true[%d]", i)]], 0.025))
log_r_true_hi  <- sapply(seq_len(sd_$I),
                          function(i) quantile(draws[[sprintf("log_r_true[%d]", i)]], 0.975))
beta_react_med <- median(draws$beta_react)

child_input <- videos %>%
  group_by(child_ii) %>%
  summarise(n_recordings = n(),
            mean_log_r_obs = mean(log_r_obs),
            sd_log_r_obs = sd(log_r_obs),
            .groups = "drop") %>%
  rename(ii = child_ii) %>%
  mutate(log_r_true_med = log_r_true_med[ii],
         log_r_true_lo  = log_r_true_lo[ii],
         log_r_true_hi  = log_r_true_hi[ii])

p_input <- ggplot() +
  # Per-video log_r_obs as small points
  geom_jitter(data = videos %>%
                rename(ii = child_ii) %>%
                mutate(log_r_true_med = log_r_true_med[ii]),
              aes(x = log_r_true_med, y = log_r_obs),
              alpha = 0.18, size = 0.4, width = 0.05) +
  # Per-child mean log_r_obs as red points
  geom_point(data = child_input,
             aes(x = log_r_true_med, y = mean_log_r_obs),
             color = "firebrick", size = 1.5) +
  # 95% CrI on log_r_true
  geom_errorbarh(data = child_input,
                 aes(y = mean_log_r_obs,
                     xmin = log_r_true_lo, xmax = log_r_true_hi),
                 color = "firebrick", height = 0, alpha = 0.5) +
  geom_abline(slope = 1, intercept = beta_react_med,
              linetype = "dashed", color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "gray50") +
  labs(x = expression(log~r[i]^{true} ~ "(posterior median, 95% CrI)"),
       y = expression(log~r[iv]^{obs} ~ "(per video)"),
       title = "Per-child observed input vs inferred true rate",
       subtitle = sprintf(paste0("Dashed green: y = x + beta_react (= %.2f); ",
                                  "dotted: y = x. ",
                                  "Reactivity inflation = %.0f%% (CrI: %.0f-%.0f%%)"),
                          beta_react_med,
                          (exp(beta_react_med) - 1) * 100,
                          (exp(quantile(draws$beta_react, 0.025)) - 1) * 100,
                          (exp(quantile(draws$beta_react, 0.975)) - 1) * 100)) +
  theme_minimal(base_size = 11)

# ---- 3. Population vocab trajectory (smooth grid) + observed bins ---- #
admin_info <- bundle$admin_info
df <- bundle$df %>% as_tibble()

# Population parameters. log_irt_io.stan keeps xi and zeta in
# independent priors (no MVN, so rho = 0 by construction).
# sigma_xi is not declared as a derived quantity, so reconstruct it.
sigma_alpha <- median(draws$sigma_alpha)
sigma_r     <- if ("sigma_r" %in% names(draws)) median(draws$sigma_r) else sd_$sigma_r
mu_r        <- if ("mu_r" %in% names(draws)) median(draws$mu_r) else sd_$mu_r
sigma_xi    <- sqrt(sigma_alpha^2 + sigma_r^2)
sigma_zeta  <- median(draws$sigma_zeta)
rho_xi_zeta <- 0    # independent in io model

psi_med  <- sapply(seq_len(sd_$J),
                   function(j) median(draws[[sprintf("psi[%d]", j)]]))
log_lambda_med <- if ("log_lambda[1]" %in% names(draws)) {
  sapply(seq_len(sd_$J),
         function(j) median(draws[[sprintf("log_lambda[%d]", j)]]))
} else {
  rep(0, sd_$J)
}
lambda_med <- exp(log_lambda_med)
log_p <- log(bundle$word_info$prob)
log_H <- sd_$log_H
a0    <- sd_$a0
s_med <- median(draws$s)
delta_med <- median(draws$delta)

# Sample population children
N_DRAWS  <- 500
AGE_GRID <- seq(8, 32, by = 0.25)
set.seed(20260429)
if (sigma_zeta > 0.01) {
  Sigma <- matrix(c(sigma_xi^2,
                    rho_xi_zeta * sigma_xi * sigma_zeta,
                    rho_xi_zeta * sigma_xi * sigma_zeta,
                    sigma_zeta^2), 2, 2)
  Z <- MASS::mvrnorm(N_DRAWS, mu = c(mu_r, 0), Sigma = Sigma)
  xi_draws   <- Z[, 1]; zeta_draws <- Z[, 2]
} else {
  xi_draws   <- rnorm(N_DRAWS, mu_r, sigma_xi)
  zeta_draws <- rep(0, N_DRAWS)
}

# Population trajectory under the *prior* MVN (fresh hypothetical kids):
# this is what the model says future kids look like.
prior_traj <- bind_rows(lapply(AGE_GRID, function(a) {
  ae <- max(a - s_med, 0.01); la <- log(ae / a0)
  base_n <- xi_draws + log_H + (1 + delta_med + zeta_draws) * la
  eta <- outer(base_n, log_p - psi_med, "+")
  eta <- sweep(eta, 2, lambda_med, "*")
  vocab_per_draw <- rowSums(plogis(eta))
  tibble(age = a,
         mean = mean(vocab_per_draw),
         lo10 = quantile(vocab_per_draw, 0.10),
         hi90 = quantile(vocab_per_draw, 0.90))
}))

# Empirical "fitted-kids" trajectory: each of the I=20 actual kids'
# fitted growth curve, averaged across kids at every age. Uses each
# kid's posterior-mean (xi_i, zeta_i). Different from the prior MVN
# when the model is partially under-identified between delta and
# mean(zeta_i).
xi_kids   <- sapply(seq_len(sd_$I),
                    function(i) median(draws[[sprintf("xi[%d]", i)]]))
zeta_kids <- sapply(seq_len(sd_$I),
                    function(i) median(draws[[sprintf("zeta[%d]", i)]]))
emp_traj <- bind_rows(lapply(AGE_GRID, function(a) {
  ae <- max(a - s_med, 0.01); la <- log(ae / a0)
  base_n <- xi_kids + log_H + (1 + delta_med + zeta_kids) * la
  eta <- outer(base_n, log_p - psi_med, "+")
  eta <- sweep(eta, 2, lambda_med, "*")
  vocab_per_kid <- rowSums(plogis(eta))
  tibble(age = a,
         mean = mean(vocab_per_kid),
         lo10 = quantile(vocab_per_kid, 0.10),
         hi90 = quantile(vocab_per_kid, 0.90))
}))

# Observed: bin admin totals by integer age
obs_bins <- admin_info %>% as_tibble() %>%
  mutate(obs_total = sapply(aa,
                            function(k) sum(df$produces[df$aa == k]))) %>%
  mutate(age_bin = round(age)) %>%
  group_by(age_bin) %>%
  summarise(n = n(), obs_mean = mean(obs_total), .groups = "drop") %>%
  filter(n >= 2)   # BabyView has fewer admins per age bin

# Also show every individual admin (translucent)
all_admins <- admin_info %>% as_tibble() %>%
  mutate(obs_total = sapply(aa,
                            function(k) sum(df$produces[df$aa == k])))

p_traj <- ggplot() +
  # Empirical fitted-kids ribbon
  geom_ribbon(data = emp_traj,
              aes(x = age, ymin = lo10, ymax = hi90),
              fill = "steelblue", alpha = 0.20) +
  geom_line(data = emp_traj, aes(x = age, y = mean,
                                  color = "fitted (mean across 20 kids)"),
            linewidth = 0.9) +
  geom_line(data = prior_traj, aes(x = age, y = mean,
                                    color = "prior MVN (mu_r, 0)"),
            linewidth = 0.9, linetype = "dashed") +
  # Observed admin totals
  geom_point(data = all_admins,
             aes(x = age, y = obs_total,
                 color = "observed admin"),
             alpha = 0.25, size = 0.9) +
  geom_point(data = obs_bins,
             aes(x = age_bin, y = obs_mean,
                 color = "observed bin mean"),
             size = 1.8) +
  scale_color_manual(values = c("fitted (mean across 20 kids)" = "steelblue",
                                "prior MVN (mu_r, 0)" = "darkorange",
                                "observed admin" = "firebrick",
                                "observed bin mean" = "firebrick"),
                     name = NULL) +
  labs(x = "age (months)", y = "mean vocab (out of J=200)",
       title = "SEEDLingS io_slopes: trajectories",
       subtitle = sprintf(paste0("Solid steel: average of 20 actual fitted kids' ",
                                  "growth curves.  ",
                                  "Dashed orange: hypothetical kid drawn from ",
                                  "prior MVN(mu_r, Sigma).  ",
                                  "Gap reveals delta vs mean(zeta) under-",
                                  "identification (N=%d kids)."),
                          sd_$I)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# ---- 4. Compose ---- #
out <- (p_scalars / p_input / p_traj) +
  plot_layout(heights = c(1, 1.2, 1.2))
ggsave(file.path(PATHS$figs_dir, "io", "seedlings_io_summary.png"),
       out, width = 9, height = 12, dpi = 150)
cat("Wrote outputs/figs/seedlings_io_summary.png\n")

# Save individual plots too for clean inclusion
ggsave(file.path(PATHS$figs_dir, "io", "seedlings_io_scalars.png"),
       p_scalars, width = 7, height = 4, dpi = 150)
ggsave(file.path(PATHS$figs_dir, "io", "seedlings_io_input.png"),
       p_input, width = 8, height = 5, dpi = 150)
ggsave(file.path(PATHS$figs_dir, "io", "seedlings_io_trajectory.png"),
       p_traj, width = 8, height = 5, dpi = 150)
cat("Wrote individual panel files too.\n")
