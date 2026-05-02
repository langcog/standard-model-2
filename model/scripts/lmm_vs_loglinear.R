## LMM (linear-in-age) vs. log-linear (1+delta)*log_age comparison.
##
## Predicted population vocabulary trajectory under both models,
## overlaid with observed admin-level bin means. If LMM and log-linear
## fit the data equally well, the (1+delta) coefficient is buying
## flexibility, not structural content.
##
## Output: model/figs/longitudinal/lmm_vs_loglinear.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(posterior); library(rstan)
})

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

bundle <- load_dataset_bundle("english")
sd_    <- bundle$stan_data
log_p  <- log(bundle$word_info$prob)
log_H  <- sd_$log_H
a0     <- sd_$a0

# ---- Pull both fits ---- #
fit_ll  <- readRDS(file.path(PATHS$fits_dir, "long_slopes.rds"))
fit_lmm <- readRDS(file.path(PATHS$fits_dir, "long_lmm_slopes.rds"))
d_ll    <- as_draws_df(fit_ll)
d_lmm   <- as_draws_df(fit_lmm)

# ---- Helper: predict population vocab as a function of age ---- #
# For each model, sample N_DRAWS hypothetical kids from the posterior
# MVN and compute predicted vocab at each age in a grid.
N_DRAWS  <- 500
AGE_GRID <- seq(8, 32, by = 0.25)
set.seed(20260502)

predict_vocab <- function(d, model_kind = c("loglinear", "lmm")) {
  model_kind <- match.arg(model_kind)
  med <- function(p) median(d[[p]])
  mu_r       <- sd_$mu_r
  sigma_xi   <- med("sigma_xi")
  sigma_zeta <- med("sigma_zeta")
  rho        <- med("rho_xi_zeta")
  Sigma <- matrix(c(sigma_xi^2,
                    rho * sigma_xi * sigma_zeta,
                    rho * sigma_xi * sigma_zeta,
                    sigma_zeta^2), 2, 2)
  Z <- MASS::mvrnorm(N_DRAWS, mu = c(mu_r, 0), Sigma = Sigma)
  xi   <- Z[, 1]
  zeta <- Z[, 2]
  psi  <- sapply(seq_len(sd_$J),
                 function(j) median(d[[sprintf("psi[%d]", j)]]))

  bind_rows(lapply(AGE_GRID, function(a) {
    if (model_kind == "loglinear") {
      s     <- med("s"); delta <- med("delta")
      log_age <- log(max(a - s, 0.01) / a0)
      base_n  <- xi + log_H + (1 + delta + zeta) * log_age
    } else {
      beta_age <- med("beta_age")
      dt <- a - a0
      base_n <- xi + log_H + (beta_age + zeta) * dt
    }
    eta <- outer(base_n, log_p - psi, "+")     # N_DRAWS x J
    vocab <- rowSums(plogis(eta))
    tibble(age = a, mean = mean(vocab),
           lo10 = quantile(vocab, 0.10),
           hi90 = quantile(vocab, 0.90))
  })) %>% mutate(model = model_kind)
}

cat("Computing population trajectories...\n")
pred_ll  <- predict_vocab(d_ll,  "loglinear")
pred_lmm <- predict_vocab(d_lmm, "lmm")
pred <- bind_rows(pred_ll, pred_lmm) %>%
  mutate(model = factor(model,
                        levels = c("loglinear", "lmm"),
                        labels = c("log-linear (1+delta)*log_age",
                                   "LMM (linear-in-age)")))

# ---- Observed bin means for reference ---- #
df <- bundle$df %>% as_tibble()
admin_info <- bundle$admin_info %>% as_tibble()
admin_totals <- df %>% group_by(aa) %>%
  summarise(produced = sum(produces), .groups = "drop") %>%
  left_join(admin_info %>% select(aa, age), by = "aa")

obs_bin <- admin_totals %>% mutate(age_bin = round(age)) %>%
  group_by(age_bin) %>%
  summarise(n = n(), obs_mean = mean(produced),
            obs_se = sd(produced) / sqrt(n),
            .groups = "drop") %>%
  filter(n >= 5)

# ---- Panel 1: population trajectory comparison ---- #
p_traj <- ggplot() +
  geom_ribbon(data = pred,
              aes(x = age, ymin = lo10, ymax = hi90, fill = model),
              alpha = 0.18) +
  geom_line(data = pred,
            aes(x = age, y = mean, color = model),
            linewidth = 1.0) +
  geom_pointrange(data = obs_bin,
                  aes(x = age_bin, y = obs_mean,
                      ymin = obs_mean - obs_se, ymax = obs_mean + obs_se),
                  color = "black", size = 0.25) +
  scale_color_manual(values = c("log-linear (1+delta)*log_age" = "#1f77b4",
                                 "LMM (linear-in-age)" = "#d62728")) +
  scale_fill_manual(values = c("log-linear (1+delta)*log_age" = "#1f77b4",
                                "LMM (linear-in-age)" = "#d62728")) +
  labs(x = "age (months)",
       y = "predicted vocab (out of J = 200)",
       title = "(1) Population mean vocab trajectory: log-linear vs. LMM",
       subtitle = paste0("Black points: observed admin-level mean ",
                          "+/- SE.  Ribbons: 80% population spread ",
                          "from prior MVN.")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

# ---- Panel 2: per-kid posterior alignment ---- #
# Compare per-kid xi_i and zeta_i (in their respective parameterizations)
xi_ll   <- sapply(seq_len(sd_$I),
                  function(i) median(d_ll[[sprintf("xi[%d]", i)]]))
zeta_ll <- sapply(seq_len(sd_$I),
                  function(i) median(d_ll[[sprintf("zeta[%d]", i)]]))
xi_lmm  <- sapply(seq_len(sd_$I),
                  function(i) median(d_lmm[[sprintf("xi[%d]", i)]]))
zeta_lmm <- sapply(seq_len(sd_$I),
                   function(i) median(d_lmm[[sprintf("zeta[%d]", i)]]))

# Convert log-linear zeta to "logits/month at a_0" for direct comparison:
# d theta / dt at t=a_0 in log-linear is (1+delta+zeta)/a_0.
# In LMM it's (beta_age + zeta).
# So compare the per-kid implied slope at a_0:
#   loglinear:  (1 + delta + zeta_ll) / a_0
#   lmm:        beta_age + zeta_lmm
delta_ll <- median(d_ll$delta); s_ll <- median(d_ll$s)
beta_age <- median(d_lmm$beta_age)
slope_ll  <- (1 + delta_ll + zeta_ll) / a0
slope_lmm <- beta_age + zeta_lmm

per_kid <- tibble(
  ii = seq_len(sd_$I),
  xi_ll = xi_ll, xi_lmm = xi_lmm,
  slope_ll = slope_ll, slope_lmm = slope_lmm
)

p_xi <- ggplot(per_kid, aes(x = xi_ll, y = xi_lmm)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(alpha = 0.6, size = 1.4, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick",
              linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("r = %+.3f", cor(xi_ll, xi_lmm)),
           hjust = -0.1, vjust = 1.5, color = "firebrick", size = 3.6) +
  labs(x = expression(xi[i] ~ "  (log-linear)"),
       y = expression(xi[i] ~ "  (LMM)"),
       title = "(2a) Per-child intercept agreement") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

p_slope <- ggplot(per_kid, aes(x = slope_ll, y = slope_lmm)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(alpha = 0.6, size = 1.4, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick",
              linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("r = %+.3f", cor(slope_ll, slope_lmm)),
           hjust = -0.1, vjust = 1.5, color = "firebrick", size = 3.6) +
  labs(x = "implied slope at a_0 (logits/mo, log-linear)",
       y = "implied slope at a_0 (logits/mo, LMM)",
       title = "(2b) Per-child slope agreement at a_0") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

out <- p_traj / (p_xi | p_slope) + plot_layout(heights = c(1.4, 1))
ggsave(file.path(OUT_FIGS, "lmm_vs_loglinear.png"),
       out, width = 11, height = 9, dpi = 200)
cat("Wrote lmm_vs_loglinear.png\n")

cat(sprintf("\nPer-kid agreement:\n  r(xi_ll, xi_lmm)         = %+.3f\n  r(slope_ll, slope_lmm)  = %+.3f\n",
            cor(xi_ll, xi_lmm), cor(slope_ll, slope_lmm)))
