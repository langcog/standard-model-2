## Analyze io_slopes / babyview fit. Three plots:
##  1. Per-child observed input distributions vs inferred log_r_true
##  2. Per-recording log_r_obs scatter (with reactivity inflation)
##  3. Per-child predicted vs observed vocab trajectory + scalar summaries
##
## Output: model/figs/babyview_io_*.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(posterior); library(rstan)
})

fit_path <- file.path(PATHS$fits_dir, "io_slopes.rds")
fit    <- readRDS(fit_path)
bundle <- load_dataset_bundle("babyview")
sd_    <- bundle$stan_data
draws  <- as_draws_df(fit)

cat(sprintf("BabyView io_slopes: I=%d, A=%d, V=%d, J=%d, N=%d\n",
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
       title = "BabyView io_slopes: scalar posteriors",
       subtitle = sprintf("I=%d kids, V=%d videos, J=%d items, N=%d obs",
                          sd_$I, sd_$V, sd_$J, sd_$N)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(family = "mono", size = 10))

# ---- 2. Per-recording observed log_r_obs vs per-child inferred ---- #
videos <- bundle$videos %>% as_tibble()
log_r_true_med <- sapply(seq_len(sd_$I),
                          function(i) median(draws[[sprintf("log_r_true[%d]", i)]]))
log_r_true_lo  <- sapply(seq_len(sd_$I),
                          function(i) quantile(draws[[sprintf("log_r_true[%d]", i)]], 0.025))
log_r_true_hi  <- sapply(seq_len(sd_$I),
                          function(i) quantile(draws[[sprintf("log_r_true[%d]", i)]], 0.975))
beta_react_med <- median(draws$beta_react)

child_input <- videos %>%
  group_by(child_ii) %>%
  summarise(n_videos = n(),
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

# ---- 3. Per-child vocab trajectory: fitted vs observed ---- #
admin_info <- bundle$admin_info
df <- bundle$df %>% as_tibble()

# Compute predicted total per admin via posterior medians (cheap)
xi_med <- sapply(seq_len(sd_$I),
                  function(i) median(draws[[sprintf("xi[%d]", i)]]))
zeta_med <- sapply(seq_len(sd_$I),
                   function(i) median(draws[[sprintf("zeta[%d]", i)]]))
psi_med  <- sapply(seq_len(sd_$J),
                   function(j) median(draws[[sprintf("psi[%d]", j)]]))
log_p <- log(bundle$word_info$prob)
log_H <- sd_$log_H
a0    <- sd_$a0
s_med <- median(draws$s)
delta_med <- median(draws$delta)

ai <- pmax(admin_info$age - s_med, 0.01)
log_age <- log(ai / a0)
pred_total <- numeric(nrow(admin_info))
obs_total  <- numeric(nrow(admin_info))
for (k in seq_len(nrow(admin_info))) {
  i <- admin_info$ii[k]
  eta <- xi_med[i] + log_p + log_H +
         (1 + delta_med + zeta_med[i]) * log_age[k] - psi_med
  pred_total[k] <- sum(plogis(eta))
  obs_total[k]  <- sum(df$produces[df$aa == admin_info$aa[k]])
}

per_admin <- admin_info %>%
  as_tibble() %>%
  mutate(pred_total = pred_total, obs_total = obs_total)

# Per-child trajectory (line per kid)
p_traj <- ggplot(per_admin, aes(x = age, group = ii)) +
  geom_line(aes(y = pred_total), color = "steelblue", alpha = 0.4) +
  geom_point(aes(y = obs_total), color = "firebrick",
             alpha = 0.6, size = 1) +
  geom_line(aes(y = obs_total), color = "firebrick",
            alpha = 0.4, linetype = "dotted") +
  labs(x = "age (months)", y = "vocab (out of J=200)",
       title = "Per-child trajectories: fitted vs observed",
       subtitle = sprintf("Blue lines: fitted growth curves; red dots: observed admin totals (N=%d kids)",
                          sd_$I)) +
  theme_minimal(base_size = 11)

# ---- 4. Compose ---- #
out <- (p_scalars / p_input / p_traj) +
  plot_layout(heights = c(1, 1.2, 1.2))
ggsave(file.path(PATHS$figs_dir, "babyview_io_summary.png"),
       out, width = 9, height = 12, dpi = 150)
cat("Wrote model/figs/babyview_io_summary.png\n")

# Save individual plots too for clean inclusion
ggsave(file.path(PATHS$figs_dir, "babyview_io_scalars.png"),
       p_scalars, width = 7, height = 4, dpi = 150)
ggsave(file.path(PATHS$figs_dir, "babyview_io_input.png"),
       p_input, width = 8, height = 5, dpi = 150)
ggsave(file.path(PATHS$figs_dir, "babyview_io_trajectory.png"),
       p_traj, width = 8, height = 5, dpi = 150)
cat("Wrote individual panel files too.\n")
