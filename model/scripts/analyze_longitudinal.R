## Analyze the longitudinal 2PL+slopes fit.
##
## Produces four plots:
##   1. Posterior densities of headline scalars (sigma_alpha, sigma_zeta,
##      rho_xi_zeta, pi_alpha).
##   2. Scatter of posterior-mean (xi_i, zeta_i) with marginal densities —
##      the direct visual of the intercept-slope coupling.
##   3. Predicted vs. observed vocabulary trajectories (per-child sample).
##   4. Variance-of-ability vs. age (variance grows with age if sigma_zeta > 0).

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

fit    <- readRDS(file.path(PATHS$fits_dir, "long_2pl_slopes.rds"))
bundle <- readRDS(file.path(PATHS$fits_dir, "long_subset_data.rds"))

draws <- as_draws_df(fit)
sd_   <- bundle$stan_data
I <- sd_$I; A <- sd_$A; J <- sd_$J

## ---- 1. Headline scalar posteriors ----
pars <- c("sigma_alpha", "sigma_zeta", "rho_xi_zeta", "pi_alpha", "s", "delta")
scalars <- tibble::tibble(
  param = pars,
  median = sapply(pars, function(p) median(draws[[p]])),
  lo = sapply(pars, function(p) quantile(draws[[p]], .025)),
  hi = sapply(pars, function(p) quantile(draws[[p]], .975))
)
cat("\n--- Headline scalars ---\n"); print(scalars, digits = 3)

p1 <- draws %>%
  select(all_of(pars)) %>%
  tidyr::pivot_longer(everything(), names_to = "param", values_to = "val") %>%
  ggplot(aes(val)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  facet_wrap(~param, scales = "free", ncol = 3) +
  labs(title = "Longitudinal 2PL+slopes posteriors",
       x = NULL, y = "density") +
  theme_minimal(base_size = 11)
ggsave(file.path(PATHS$figs_dir, "long_2pl_slopes_1_scalars.png"),
       p1, width = 10, height = 5, dpi = 150)

## ---- 2. (xi, zeta) scatter with marginal densities ----
xi_cols   <- grep("^xi\\[",   names(draws), value = TRUE)
zeta_cols <- grep("^zeta\\[", names(draws), value = TRUE)
xi_mean   <- colMeans(draws[, xi_cols])
zeta_mean <- colMeans(draws[, zeta_cols])

xz <- tibble::tibble(xi = xi_mean, zeta = zeta_mean)
r_obs <- cor(xi_mean, zeta_mean)
p2 <- ggplot(xz, aes(xi, zeta)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = FALSE, colour = "firebrick",
              linewidth = 0.6) +
  labs(title = "Per-child (xi, zeta) posterior means",
       subtitle = sprintf("r(xi, zeta) = %.2f  |  posterior rho = %.2f",
                          r_obs, median(draws$rho_xi_zeta)),
       x = expression(xi[i]~"(intercept, log-token units)"),
       y = expression(zeta[i]~"(slope deviation)")) +
  theme_minimal(base_size = 11)
ggsave(file.path(PATHS$figs_dir, "long_2pl_slopes_2_xi_zeta_scatter.png"),
       p2, width = 5.5, height = 5, dpi = 150)

## ---- 3. Per-child vocab trajectories, sample ----
# Actual observed trajectories: sum y per (child, admin), colour by child_id
obs_admin <- tibble::tibble(aa = sd_$aa, y = sd_$y) %>%
  group_by(aa) %>%
  summarise(vocab = sum(y), .groups = "drop") %>%
  mutate(age = sd_$admin_age[aa],
         ii  = sd_$admin_to_child[aa])

# Select a random subset of children for plotting
set.seed(42)
sample_ii <- sample(unique(obs_admin$ii),
                    size = min(40, length(unique(obs_admin$ii))))
obs_sub <- obs_admin %>% filter(ii %in% sample_ii)

p3 <- ggplot(obs_sub, aes(age, vocab, group = ii)) +
  geom_line(alpha = 0.35, colour = "grey30") +
  geom_point(alpha = 0.7, size = 1.2) +
  labs(title = sprintf("Observed trajectories (random %d children)",
                       length(sample_ii)),
       subtitle = "Each line = one child's vocabulary as it grows",
       x = "Age (months)", y = "Words produced (of 200)") +
  theme_minimal(base_size = 11)
ggsave(file.path(PATHS$figs_dir, "long_2pl_slopes_3_trajectories.png"),
       p3, width = 6, height = 4.5, dpi = 150)

## ---- 4. Variance-of-ability vs. age ----
# Model says: Var(theta_i | age) = sigma_xi^2 + sigma_zeta^2 * L(a)^2 +
#             2 * rho * sigma_xi * sigma_zeta * L(a)
# (the covariance term matters when rho != 0)
s_d   <- draws$s
sig_x <- draws$sigma_xi
sig_z <- draws$sigma_zeta
rho   <- draws$rho_xi_zeta

ages <- seq(15, 35, by = 0.5)
pred_df <- purrr::map_dfr(ages, function(a) {
  L <- log(pmax(a - s_d, 0.01) / MODEL_CONSTANTS$a0)
  var_theta <- sig_x^2 + 2 * rho * sig_x * sig_z * L + sig_z^2 * L^2
  var_theta <- pmax(var_theta, 1e-6)
  tibble::tibble(age = a,
                 sd_med = median(sqrt(var_theta)),
                 sd_lo = quantile(sqrt(var_theta), .025),
                 sd_hi = quantile(sqrt(var_theta), .975))
})

# Observed: within-age bin SD of logit(vocab/J+1)
obs_by_age <- obs_admin %>%
  mutate(age_bin = round(age),
         lprop = qlogis(pmin(pmax((vocab + 0.5) / (J + 1), 1e-4), 1 - 1e-4))) %>%
  group_by(age_bin) %>%
  summarise(n = n(), sd_lprop = sd(lprop), .groups = "drop") %>%
  filter(n >= 5)

p4 <- ggplot(pred_df, aes(age)) +
  geom_ribbon(aes(ymin = sd_lo, ymax = sd_hi), alpha = 0.25,
              fill = "steelblue") +
  geom_line(aes(y = sd_med), colour = "steelblue", linewidth = 0.8) +
  geom_point(data = obs_by_age, aes(age_bin, sd_lprop, size = n),
             alpha = 0.7, show.legend = FALSE) +
  scale_size_area(max_size = 4) +
  labs(title = "Within-age ability SD vs. age",
       subtitle = "Dots = observed SD of logit(vocab prop); ribbon = model sqrt(Var(theta|age))",
       x = "Age (months)", y = "SD of ability (logit scale)") +
  theme_minimal(base_size = 11)
ggsave(file.path(PATHS$figs_dir, "long_2pl_slopes_4_variance_by_age.png"),
       p4, width = 7, height = 4.5, dpi = 150)

## Combined overview
combined <- (p1 / (p2 | p3) / p4) +
  plot_annotation(title = "Longitudinal 2PL+slopes: posterior summary",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          size = 14)))
ggsave(file.path(PATHS$figs_dir, "long_2pl_slopes_all.png"),
       combined, width = 11, height = 13, dpi = 120)

cat("\nFigures saved to model/figs/long_2pl_slopes_*.png\n")
