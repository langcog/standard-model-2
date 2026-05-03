## Compare the English longitudinal ablation set in three plots:
##  1. Forest plot of scalar posteriors across variants
##  2. Population mean vocab trajectory: fitted vs observed by variant
##  3. Per-child posterior xi/zeta scatter, faceted by variant
##
## Inputs (5 fits):
##   fits/long_slopes.rds            (lean reference)
##   fits/long_baseline.rds          (drops slopes)
##   fits/long_fix_delta_slopes.rds  (pins delta=0)
##   fits/long_free_s_slopes.rds     (frees s)
##   fits/long_2pl_slopes.rds        (adds 2PL)
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

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

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

ggsave(file.path(PATHS$figs_dir, "longitudinal", "english_ablations_scalars.png"),
       p_scalars, width = 11, height = 6, dpi = 150)
cat("Wrote english_ablations_scalars.png\n")

# ---- 3. Population-level vocab trajectory (smooth grid) ---- #
# For each variant, predict expected vocab at each age in a regular
# grid by integrating over the population distribution of (xi, zeta).
# Uses posterior medians of all parameters; samples N_DRAWS children
# from MVN(mu_r, Sigma) where Sigma comes from sigma_xi, sigma_zeta,
# rho_xi_zeta. This is the "fixed effects + population variance"
# prediction, not the per-admin (which would condition on each
# individual child's posterior).
N_DRAWS <- 500
AGE_GRID <- seq(11, 32, by = 0.25)

get_population_trajectory <- function(variant) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) return(NULL)
  fit <- readRDS(path)
  draws <- as_draws_df(fit)
  med <- function(p) median(draws[[p]])

  mu_r        <- sd_$mu_r       # data input, not a parameter
  sigma_xi    <- med("sigma_xi")
  has_zeta    <- "sigma_zeta" %in% names(draws) &&
                 median(draws$sigma_zeta) > 0.01
  sigma_zeta  <- if (has_zeta) med("sigma_zeta") else 0
  rho_xi_zeta <- if (has_zeta && "rho_xi_zeta" %in% names(draws))
                   med("rho_xi_zeta") else 0

  psi <- sapply(seq_len(sd_$J),
                function(j) median(draws[[sprintf("psi[%d]", j)]]))
  log_lambda <- if ("log_lambda[1]" %in% names(draws))
                  sapply(seq_len(sd_$J),
                         function(j) median(draws[[sprintf("log_lambda[%d]", j)]]))
                else rep(0, sd_$J)
  lambda <- exp(log_lambda)

  s     <- med("s")
  delta <- med("delta")
  log_p <- log(bundle$word_info$prob)
  log_H <- sd_$log_H
  a0    <- sd_$a0

  # Sample N_DRAWS children from population distribution. If
  # sigma_zeta is effectively pinned at 0 (no slopes), draw only xi.
  if (sigma_zeta > 0.01) {
    Sigma <- matrix(c(sigma_xi^2,
                      rho_xi_zeta * sigma_xi * sigma_zeta,
                      rho_xi_zeta * sigma_xi * sigma_zeta,
                      sigma_zeta^2), 2, 2)
    Z <- MASS::mvrnorm(N_DRAWS, mu = c(mu_r, 0), Sigma = Sigma)
    xi_draws   <- Z[, 1]
    zeta_draws <- Z[, 2]
  } else {
    xi_draws   <- rnorm(N_DRAWS, mu_r, sigma_xi)
    zeta_draws <- rep(0, N_DRAWS)
  }

  # For each age + each draw, compute predicted vocab total
  # Vectorize: for each age, build matrix [N_DRAWS x J] of eta values
  out_rows <- vector("list", length(AGE_GRID))
  for (k in seq_along(AGE_GRID)) {
    a   <- AGE_GRID[k]
    ae  <- max(a - s, 0.01)
    la  <- log(ae / a0)
    # eta_{nj} = lambda_j * (xi_n + log_p_j + log_H +
    #                        (1 + delta + zeta_n) * la - psi_j)
    base_n <- xi_draws + log_H + (1 + delta + zeta_draws) * la
    # outer term: base_n  +  log_p[j] - psi[j]
    eta <- outer(base_n, log_p - psi, "+")    # N_DRAWS x J
    eta <- sweep(eta, 2, lambda, "*")
    p   <- plogis(eta)
    vocab_per_draw <- rowSums(p)
    out_rows[[k]] <- tibble(age = a,
                            mean = mean(vocab_per_draw),
                            lo10 = quantile(vocab_per_draw, 0.10),
                            hi90 = quantile(vocab_per_draw, 0.90))
  }
  bind_rows(out_rows) %>% mutate(variant = variant)
}

cat("Computing population-level trajectories (prior MVN)...\n")
set.seed(20260429)
pop_traj <- bind_rows(lapply(names(VARIANTS), function(v) {
  cat("  ", v, "\n")
  get_population_trajectory(v)
})) %>%
  mutate(variant_label = factor(VARIANTS[variant],
                                levels = unname(VARIANTS)),
         source = "prior MVN (mu_r, 0)")

# Empirical fitted-kids trajectory: each actual kid's posterior-mean
# growth curve, averaged at every grid age. Reveals delta vs mean(zeta)
# under-identification when this differs from the prior-MVN curve above.
get_empirical_trajectory <- function(variant) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) return(NULL)
  fit <- readRDS(path)
  draws <- as_draws_df(fit)
  med <- function(p) median(draws[[p]])

  xi_kids   <- sapply(seq_len(sd_$I),
                       function(i) median(draws[[sprintf("xi[%d]", i)]]))
  zeta_kids <- if ("zeta[1]" %in% names(draws))
                 sapply(seq_len(sd_$I),
                        function(i) median(draws[[sprintf("zeta[%d]", i)]]))
               else rep(0, sd_$I)
  psi <- sapply(seq_len(sd_$J),
                function(j) median(draws[[sprintf("psi[%d]", j)]]))
  log_lambda <- if ("log_lambda[1]" %in% names(draws)) {
    sapply(seq_len(sd_$J),
           function(j) median(draws[[sprintf("log_lambda[%d]", j)]]))
  } else { rep(0, sd_$J) }
  lambda <- exp(log_lambda)
  s     <- med("s"); delta <- med("delta")
  log_p <- log(bundle$word_info$prob)
  log_H <- sd_$log_H; a0 <- sd_$a0

  bind_rows(lapply(AGE_GRID, function(a) {
    ae <- max(a - s, 0.01); la <- log(ae / a0)
    base_n <- xi_kids + log_H + (1 + delta + zeta_kids) * la
    eta <- outer(base_n, log_p - psi, "+")
    eta <- sweep(eta, 2, lambda, "*")
    vocab_per_kid <- rowSums(plogis(eta))
    tibble(age = a, mean = mean(vocab_per_kid),
           lo10 = quantile(vocab_per_kid, 0.10),
           hi90 = quantile(vocab_per_kid, 0.90))
  })) %>% mutate(variant = variant)
}

cat("Computing empirical fitted-kids trajectories...\n")
emp_traj <- bind_rows(lapply(names(VARIANTS), function(v) {
  cat("  ", v, "\n")
  get_empirical_trajectory(v)
})) %>%
  mutate(variant_label = factor(VARIANTS[variant],
                                levels = unname(VARIANTS)),
         source = "fitted (mean of N kids)")

# Observed: bin actual admin totals by integer age
obs_bins <- df %>%
  group_by(aa) %>%
  summarise(age = first(admin_info$age[admin_info$aa == first(aa)]),
            obs_total = sum(produces), .groups = "drop") %>%
  mutate(age_bin = round(age)) %>%
  group_by(age_bin) %>%
  summarise(n = n(), obs_mean = mean(obs_total), .groups = "drop") %>%
  filter(n >= 5)

p_traj <- ggplot() +
  # Empirical fitted-kids ribbon (steel)
  geom_ribbon(data = emp_traj,
              aes(x = age, ymin = lo10, ymax = hi90),
              fill = "steelblue", alpha = 0.18) +
  geom_line(data = emp_traj,
            aes(x = age, y = mean, color = "fitted (mean of kids)"),
            linewidth = 0.8) +
  # Prior MVN (orange dashed)
  geom_line(data = pop_traj,
            aes(x = age, y = mean, color = "prior MVN (mu_r, 0)"),
            linewidth = 0.7, linetype = "dashed") +
  geom_point(data = obs_bins,
             aes(x = age_bin, y = obs_mean,
                 color = "observed (admin mean)"),
             size = 1.2, alpha = 0.85) +
  facet_wrap(~variant_label, ncol = 3) +
  scale_color_manual(values = c("fitted (mean of kids)" = "steelblue",
                                "prior MVN (mu_r, 0)" = "darkorange",
                                "observed (admin mean)" = "firebrick"),
                     name = NULL) +
  labs(x = "age (months)", y = "mean vocab (out of J=200)",
       title = "Population vocab trajectory: fitted kids vs prior MVN vs observed",
       subtitle = paste0("Solid steel: fitted-kids mean.  ",
                         "Dashed orange: hypothetical kid from prior MVN.  ",
                         "Steel ribbon: 80% population spread.")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave(file.path(PATHS$figs_dir, "longitudinal", "english_ablations_trajectory.png"),
       p_traj, width = 10, height = 6.5, dpi = 150)
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
  ggsave(file.path(PATHS$figs_dir, "longitudinal", "english_ablations_xi_zeta.png"),
         p_xz, width = 12, height = 3.5, dpi = 150)
  cat("Wrote english_ablations_xi_zeta.png\n")
}

cat("\nDone.\n")
