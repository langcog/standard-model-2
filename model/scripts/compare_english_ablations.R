## Compare the English longitudinal ablation set in three plots:
##  1. Forest plot of scalar posteriors across variants
##  2. Population mean vocab trajectory: fitted vs observed by variant
##  3. Per-child posterior xi/zeta scatter, faceted by variant
##
## Inputs (5 fits):
##   model/fits/long_slopes.rds            (lean reference)
##   model/fits/long_baseline.rds          (drops slopes)
##   model/fits/long_fix_delta_slopes.rds  (pins delta=0)
##   model/fits/long_free_s_slopes.rds     (frees s)
##   model/fits/long_2pl_slopes.rds        (adds 2PL)
##
## Output:
##   model/figs/english_ablations_scalars.png
##   model/figs/english_ablations_trajectory.png
##   model/figs/english_ablations_xi_zeta.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(posterior); library(rstan)
})

VARIANTS <- c(
  "long_slopes"            = "Lean reference\n(slopes only)",
  "long_baseline"          = "Drop slopes\n(no zeta)",
  "long_fix_delta_slopes"  = "Pin delta = 0\n(no acceleration)",
  "long_free_s_slopes"     = "Free s\n(start time)",
  "long_2pl_slopes"        = "Add 2PL\n(item lambda)"
)

bundle <- load_dataset_bundle("english")
sd_    <- bundle$stan_data
admin_info <- bundle$admin_info %>% mutate(aa = aa)
df <- bundle$df %>% as_tibble()

# ---- 1. Extract posterior summaries from every fit ---- #
extract_scalars <- function(variant) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) {
    warning(sprintf("Missing fit: %s", path)); return(NULL)
  }
  fit <- readRDS(path)
  draws <- as_draws_df(fit)
  pars <- intersect(c("sigma_alpha","sigma_zeta","sigma_lambda",
                      "pi_alpha","s","delta","rho_xi_zeta","sigma_xi"),
                    names(draws))
  out <- tibble(
    variant = variant,
    param   = pars,
    mean    = sapply(pars, function(p) mean(draws[[p]])),
    median  = sapply(pars, function(p) median(draws[[p]])),
    lo95    = sapply(pars, function(p) quantile(draws[[p]], 0.025)),
    hi95    = sapply(pars, function(p) quantile(draws[[p]], 0.975))
  )
  out
}

cat("Extracting scalar posteriors...\n")
scalars <- bind_rows(lapply(names(VARIANTS), extract_scalars)) %>%
  mutate(variant_label = factor(VARIANTS[variant],
                                levels = unname(VARIANTS)),
         variant       = factor(variant, levels = names(VARIANTS)))
print(scalars, n = 50)

# ---- 2. Forest plot of scalars ---- #
PARAM_ORDER <- c("sigma_alpha", "sigma_xi", "pi_alpha",
                 "sigma_zeta", "rho_xi_zeta",
                 "delta", "s", "sigma_lambda")
scalars_p <- scalars %>%
  filter(param %in% PARAM_ORDER) %>%
  mutate(param = factor(param, levels = PARAM_ORDER))

p_scalars <- ggplot(scalars_p,
                    aes(y = variant_label, x = median,
                        xmin = lo95, xmax = hi95,
                        color = variant)) +
  geom_pointrange(size = 0.3) +
  facet_wrap(~param, scales = "free_x", nrow = 2) +
  scale_color_brewer(palette = "Set1", guide = "none") +
  labs(x = "posterior median, 95% CrI", y = NULL,
       title = "English longitudinal ablations: scalar posteriors",
       subtitle = "5 variants of the lean baseline") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))

ggsave(file.path(PATHS$figs_dir, "english_ablations_scalars.png"),
       p_scalars, width = 11, height = 6, dpi = 150)
cat("Wrote english_ablations_scalars.png\n")

# ---- 3. Per-admin predicted vocab vs observed ---- #
# Strategy: under each fit, compute posterior-median P(produces) for
# every (admin, item) cell, then sum within each admin to get
# predicted total. Compare to observed admin total.
get_admin_predictions <- function(variant) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) return(NULL)
  fit <- readRDS(path)
  # Posterior of eta is too large; instead use posterior medians of
  # each parameter and reconstruct eta. This is the "median fit" not
  # a proper PPC, but cheap and adequate for trajectory inspection.
  draws <- as_draws_df(fit)
  med <- function(p) median(draws[[p]])
  # Per-child / per-item / per-admin quantities:
  xi   <- sapply(seq_len(sd_$I), function(i) median(draws[[sprintf("xi[%d]", i)]]))
  zeta <- if ("zeta[1]" %in% names(draws))
            sapply(seq_len(sd_$I),
                   function(i) median(draws[[sprintf("zeta[%d]", i)]]))
          else rep(0, sd_$I)
  psi  <- sapply(seq_len(sd_$J),
                 function(j) median(draws[[sprintf("psi[%d]", j)]]))
  log_lambda <- if ("log_lambda[1]" %in% names(draws))
                  sapply(seq_len(sd_$J),
                         function(j) median(draws[[sprintf("log_lambda[%d]", j)]]))
                else rep(0, sd_$J)
  s     <- med("s")
  delta <- med("delta")
  log_p <- log(bundle$word_info$prob)
  log_H <- sd_$log_H
  a0    <- sd_$a0

  # Per-admin-item linear predictor
  # eta = lambda * (xi + log_p + log_H + (1 + delta + zeta) * log((age - s)/a0) - psi)
  ai <- sapply(admin_info$age, function(a) max(a - s, 0.01))
  log_age <- log(ai / a0)

  # Build a per-admin total predicted score by summing inv_logit(eta)
  ii_per_admin <- admin_info$ii
  pred_total <- numeric(nrow(admin_info))
  obs_total  <- numeric(nrow(admin_info))
  for (k in seq_len(nrow(admin_info))) {
    i <- ii_per_admin[k]
    eta <- exp(log_lambda) * (xi[i] + log_p + log_H +
                              (1 + delta + zeta[i]) * log_age[k] - psi)
    pred_total[k] <- sum(plogis(eta))
    # Observed: count of produces in this admin
    obs_total[k] <- sum(df$produces[df$aa == admin_info$aa[k]])
  }

  tibble(variant = variant,
         age = admin_info$age,
         pred_total = pred_total,
         obs_total  = obs_total,
         ii         = ii_per_admin)
}

cat("Computing admin-level predictions for each variant...\n")
preds <- bind_rows(lapply(names(VARIANTS), function(v) {
  cat("  ", v, "\n")
  get_admin_predictions(v)
})) %>%
  mutate(variant_label = factor(VARIANTS[variant],
                                levels = unname(VARIANTS)))

# Population mean vocab trajectory: bin by age, mean obs and pred
traj <- preds %>%
  mutate(age_bin = round(age)) %>%
  group_by(variant_label, age_bin) %>%
  summarise(n = n(),
            obs_mean  = mean(obs_total),
            pred_mean = mean(pred_total),
            obs_sd    = sd(obs_total),
            pred_sd   = sd(pred_total),
            .groups = "drop") %>%
  filter(n >= 5)

p_traj <- ggplot(traj, aes(x = age_bin)) +
  geom_ribbon(aes(ymin = pred_mean - pred_sd,
                  ymax = pred_mean + pred_sd),
              fill = "steelblue", alpha = 0.18) +
  geom_line(aes(y = pred_mean, color = "fitted"), linewidth = 0.8) +
  geom_point(aes(y = obs_mean, color = "observed"),
             size = 1.2, alpha = 0.85) +
  facet_wrap(~variant_label, nrow = 1) +
  scale_color_manual(values = c(fitted = "steelblue",
                                observed = "firebrick"),
                     name = NULL) +
  labs(x = "age (months)", y = "mean vocab (out of J=200)",
       title = "Population mean vocab trajectory: fitted vs observed",
       subtitle = sprintf("English longitudinal, 5 variants. Bin: rounded age, n>=5 admins")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave(file.path(PATHS$figs_dir, "english_ablations_trajectory.png"),
       p_traj, width = 14, height = 4, dpi = 150)
cat("Wrote english_ablations_trajectory.png\n")

# ---- 4. Per-child xi-zeta scatter (only variants with slopes) ---- #
get_xi_zeta <- function(variant) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) return(NULL)
  fit <- readRDS(path)
  draws <- as_draws_df(fit)
  if (!any(grepl("^zeta\\[", names(draws)))) return(NULL)
  xi   <- sapply(seq_len(sd_$I),
                 function(i) median(draws[[sprintf("xi[%d]", i)]]))
  zeta <- sapply(seq_len(sd_$I),
                 function(i) median(draws[[sprintf("zeta[%d]", i)]]))
  tibble(variant = variant, ii = seq_len(sd_$I), xi = xi, zeta = zeta)
}

xz <- bind_rows(lapply(names(VARIANTS), get_xi_zeta)) %>%
  mutate(variant_label = factor(VARIANTS[variant],
                                levels = unname(VARIANTS)))

if (nrow(xz) > 0) {
  rho_per_variant <- xz %>% group_by(variant_label) %>%
    summarise(rho = cor(xi, zeta), .groups = "drop")
  p_xz <- ggplot(xz, aes(xi, zeta)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "firebrick",
                linewidth = 0.6) +
    geom_text(data = rho_per_variant,
              aes(label = sprintf("r = %.2f", rho)),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
              inherit.aes = FALSE, size = 3.2) +
    facet_wrap(~variant_label, nrow = 1) +
    labs(x = expression(xi[i]), y = expression(zeta[i]),
         title = "Per-child (xi, zeta) posterior medians",
         subtitle = "Variants without per-child slopes are omitted") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(PATHS$figs_dir, "english_ablations_xi_zeta.png"),
         p_xz, width = 12, height = 3.5, dpi = 150)
  cat("Wrote english_ablations_xi_zeta.png\n")
}

cat("\nDone.\n")
