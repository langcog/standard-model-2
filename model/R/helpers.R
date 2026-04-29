## Core helpers for the log-linear IRT accumulator model.
## Loaded by every driver script via source("model/R/helpers.R").
##
## Provides:
##   load_wordbank_data()
##   load_input_rate_prior()
##   subsample_wordbank()
##   build_stan_data()
##   simulate_data()
##   fit_variant()
##   summarize_fit()
##   variant_hyperpriors()
##   plot_psi_vs_logp(), plot_class_means(), plot_posterior_density()

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(posterior)
  library(ggplot2)
  library(tibble)
})
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ===========================================================================
# Data loading
# ===========================================================================

## Load English cross-sectional Wordbank data, drawn from the
## WG-and-WS-combined long_items.rds. One admin per (child, form):
## the most recent admin per child in each form. Children may appear
## twice if they have data on both WG and WS — they're treated as
## separate "persons" for the cross-sectional fit (admin_key as id).
##
## This replaces the older WS-only engWS_preprocessed.Rdata.
load_wordbank_data <- function(
  long_items_path = file.path(PATHS$fits_dir, "long_items.rds"),
  language = "English (American)"
) {
  long <- readRDS(long_items_path)
  long %>%
    filter(language == !!language,
           !is.na(prob), prob > 0,
           !is.na(produces)) %>%
    rename(lexical_class = lexical_category) %>%
    # Treat each (child, form, age) admin as its own "person" for
    # cross-sectional purposes; collapses into a single row per child
    # if they only have one admin.
    mutate(person = paste(child_id, form, age, sep = "_")) %>%
    select(person, age, item, lexical_class, prob, produces, child_id, form)
}

load_input_rate_prior <- function(path = PATHS$input_rate,
                                  column = "adult_child_tokens_hr") {
  hr <- read.csv(path)
  x <- hr[[column]]
  x <- x[!is.na(x)]
  list(mu_r = mean(log(x)), sigma_r = sd(log(x)), n = length(x))
}

# ===========================================================================
# Subsampling + Stan data construction
# ===========================================================================

subsample_wordbank <- function(df, n_children, n_items,
                               seed = 20250420, age_bins = 5,
                               class_col = "lexical_class") {
  set.seed(seed)

  persons <- df %>% distinct(person, age) %>%
    mutate(age_bin = cut(age, breaks = age_bins)) %>%
    group_by(age_bin) %>%
    slice_sample(n = max(1, floor(n_children / age_bins))) %>%
    ungroup() %>%
    pull(person)
  if (length(persons) > n_children) persons <- sample(persons, n_children)

  items <- df %>% distinct(item, !!sym(class_col)) %>%
    group_by(!!sym(class_col)) %>%
    slice_sample(n = max(1, floor(n_items / 5))) %>%
    ungroup() %>%
    pull(item)
  if (length(items) > n_items) items <- sample(items, n_items)

  df %>% filter(person %in% persons, item %in% items)
}

build_stan_data <- function(df, prior_r,
                            priors = DEFAULT_PRIORS,
                            constants = MODEL_CONSTANTS) {
  df <- df %>%
    mutate(ii = as.integer(factor(person)),
           jj = as.integer(factor(item)),
           cc = as.integer(factor(lexical_class)))

  class_levels <- levels(factor(df$lexical_class))
  child_age <- df %>% distinct(ii, age) %>% arrange(ii) %>% pull(age)
  word_info <- df %>% distinct(jj, item, prob, cc) %>% arrange(jj)

  stan_data <- c(
    list(
      N = nrow(df),
      I = max(df$ii), J = max(df$jj), C = max(df$cc),
      ii = df$ii, jj = df$jj, cc = word_info$cc,
      y  = df$produces,
      age = child_age,
      log_p = log(word_info$prob),
      log_H = constants$log_H,
      a0    = constants$a0,
      mu_r = prior_r$mu_r,
      sigma_r = prior_r$sigma_r
    ),
    priors
  )

  list(stan_data = stan_data,
       word_info = word_info,
       class_levels = class_levels,
       child_age = child_age,
       df = df)
}

# Variants opt in to model components ABOVE the lean baseline defined
# in DEFAULT_PRIORS (Rasch + frequency + per-class psi + free delta,
# with s pinned at 0 and no per-child slopes).
#
# Naming convention: variants are named for what they ADD.
#   baseline      - lean defaults; the cross-sectional default
#   slopes        - + per-child slopes zeta (default for longitudinal)
#   2pl           - + 2PL item discrimination
#   2pl_slopes    - + both
#   free_s        - frees the start time s (RQ2 robustness)
#   fix_delta     - pins delta = 0 (RQ3 ablation; "no acceleration")
#   no_freq       - drops frequency by setting beta=0 in prep (handled
#                   in data prep, not here; left as a sentinel)
#
# Long-form variants (long_*) prefix is preserved for clarity in the
# longitudinal pipeline but resolves to the same overrides.
variant_hyperpriors <- function(name) {
  # Strip pipeline prefixes ("long_proc_", "long_", "io_") so all
  # pipelines share the same variant grammar. Order matters: longer
  # prefixes first so they match before their substrings.
  base <- sub("^(long_proc_|long_|io_)", "", name)

  switch(base,
    baseline      = list(),
    slopes        = list(sigma_zeta_prior_sd = 1),
    `2pl`         = list(sigma_lambda_prior_sd = 1),
    `2pl_slopes`  = list(sigma_lambda_prior_sd = 1,
                         sigma_zeta_prior_sd = 1),
    free_s        = list(s_prior_mean = 4.5, s_prior_sd = 2),
    free_s_slopes = list(s_prior_mean = 4.5, s_prior_sd = 2,
                         sigma_zeta_prior_sd = 1),
    fix_delta     = list(delta_prior_mean = 0, delta_prior_sd = 0.001),
    fix_delta_slopes = list(delta_prior_mean = 0, delta_prior_sd = 0.001,
                            sigma_zeta_prior_sd = 1),
    # Legacy variants for re-loading old fits / explicit comparison
    fix_s         = list(s_prior_mean = 2, s_prior_sd = 0.001),
    both_fixed    = list(delta_prior_mean = 0, delta_prior_sd = 0.001,
                         s_prior_mean = 2,    s_prior_sd = 0.001),
    `2pl_fix_delta` = list(sigma_lambda_prior_sd = 1,
                           delta_prior_mean = 0, delta_prior_sd = 0.001),
    stop(sprintf("Unknown variant: %s", name))
  )
}

# ===========================================================================
# Simulation (for recovery test)
# ===========================================================================

simulate_data <- function(I = 250, J = 150, C = 3,
                          mu_r = log(1198), sigma_r = 0.4,
                          sigma_alpha_true = 0.5,
                          mu_c_true  = c(6.5, 8.0, 9.5),
                          tau_c_true = c(0.5, 0.7, 0.7),
                          s_true = 4.5, delta_true = 0.1,
                          log_H = log(365), a0 = 20,
                          age_range = c(12, 30),
                          log_p_range = c(log(1e-5), log(1e-3)),
                          seed = 42) {
  set.seed(seed)
  stopifnot(length(mu_c_true) == C, length(tau_c_true) == C)

  cc <- sort(rep(seq_len(C), length.out = J))
  cc <- sample(cc)
  log_p <- runif(J, log_p_range[1], log_p_range[2])
  psi   <- rnorm(J, mu_c_true[cc], tau_c_true[cc])

  log_r     <- rnorm(I, mu_r, sigma_r)
  log_alpha <- rnorm(I, 0, sigma_alpha_true)
  ages      <- runif(I, age_range[1], age_range[2])

  idx <- expand.grid(ii = seq_len(I), jj = seq_len(J),
                     KEEP.OUT.ATTRS = FALSE)
  ae <- pmax(ages[idx$ii] - s_true, 0.01)
  eta <- log_r[idx$ii] + log_alpha[idx$ii] +
         log_p[idx$jj] + log_H +
         (1 + delta_true) * log(ae / a0) - psi[idx$jj]
  y <- rbinom(length(eta), 1, plogis(eta))

  list(
    obs = data.frame(ii = idx$ii, jj = idx$jj,
                     age = ages[idx$ii],
                     log_p = log_p[idx$jj], y = y),
    true = list(log_r = log_r, log_alpha = log_alpha, ages = ages,
                psi = psi, cc = cc, log_p = log_p,
                sigma_alpha = sigma_alpha_true,
                mu_c = mu_c_true, tau_c = tau_c_true,
                s = s_true, delta = delta_true,
                mu_r = mu_r, sigma_r = sigma_r),
    constants = list(log_H = log_H, a0 = a0, I = I, J = J, C = C)
  )
}

build_stan_data_from_sim <- function(sim, priors = DEFAULT_PRIORS) {
  obs <- sim$obs
  c(
    list(N = nrow(obs), I = sim$constants$I, J = sim$constants$J,
         C = sim$constants$C,
         ii = obs$ii, jj = obs$jj, cc = sim$true$cc,
         y = obs$y, age = sim$true$ages, log_p = sim$true$log_p,
         log_H = sim$constants$log_H, a0 = sim$constants$a0,
         mu_r = sim$true$mu_r, sigma_r = sim$true$sigma_r),
    priors
  )
}

# ===========================================================================
# Fitting
# ===========================================================================

fit_variant <- function(stan_data, tag,
                        cfg = DEFAULT_FIT_CONFIG,
                        model_path = PATHS$stan_model,
                        fits_dir = PATHS$fits_dir,
                        force = FALSE) {
  fit_file <- file.path(fits_dir, sprintf("%s.rds", tag))
  if (!force && file.exists(fit_file)) {
    message(sprintf("[%s] already fit at %s, loading.", tag, fit_file))
    return(readRDS(fit_file))
  }

  message(sprintf("[%s] fitting (chains=%d, iter=%d, warmup=%d)...",
                  tag, cfg$chains, cfg$iter, cfg$warmup))
  t0 <- Sys.time()
  fit <- stan(
    file    = model_path,
    data    = stan_data,
    chains  = cfg$chains,
    iter    = cfg$iter,
    warmup  = cfg$warmup,
    seed    = cfg$seed,
    control = list(adapt_delta = cfg$adapt_delta,
                   max_treedepth = cfg$max_treedepth)
  )
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("[%s] sampling time: %.1f min", tag, dt))
  saveRDS(fit, fit_file)
  fit
}

# ===========================================================================
# Summaries
# ===========================================================================

summarize_fit <- function(fit, pars = c("sigma_alpha", "s", "delta",
                                        "pi_alpha", "sigma_xi")) {
  s <- summary(fit, pars = pars)$summary
  tibble(
    param = rownames(s),
    mean = s[, "mean"],
    lo95 = s[, "2.5%"], median = s[, "50%"], hi95 = s[, "97.5%"],
    n_eff = s[, "n_eff"], Rhat = s[, "Rhat"]
  )
}

class_threshold_table <- function(fit, class_levels) {
  C <- length(class_levels)
  mu  <- summary(fit, pars = paste0("mu_c[",  seq_len(C), "]"))$summary
  tau <- summary(fit, pars = paste0("tau_c[", seq_len(C), "]"))$summary
  tibble(
    class = class_levels,
    mu_median = mu[, "50%"], mu_lo = mu[, "2.5%"], mu_hi = mu[, "97.5%"],
    tau_median = tau[, "50%"], tau_lo = tau[, "2.5%"], tau_hi = tau[, "97.5%"]
  )
}

extract_psi_df <- function(fit, word_info, class_levels) {
  draws <- as_draws_df(fit)
  psi_cols <- grep("^psi\\[", names(draws), value = TRUE)
  word_info %>%
    mutate(psi_median = sapply(psi_cols, function(p) median(draws[[p]])),
           psi_lo = sapply(psi_cols, function(p) quantile(draws[[p]], .025)),
           psi_hi = sapply(psi_cols, function(p) quantile(draws[[p]], .975)),
           class = class_levels[cc],
           log_p = log(prob))
}

# ===========================================================================
# Plots
# ===========================================================================

plot_psi_vs_logp <- function(psi_df, save_path = NULL, tag = "") {
  r <- cor(psi_df$psi_median, psi_df$log_p)
  p <- ggplot(psi_df, aes(log_p, psi_median, colour = class)) +
    geom_errorbar(aes(ymin = psi_lo, ymax = psi_hi), width = 0, alpha = 0.3) +
    geom_point(size = 1.6, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey40",
                aes(group = 1), linewidth = 0.5, linetype = "dashed") +
    labs(title = sprintf("RQ1: psi vs log p%s",
                         if (nzchar(tag)) sprintf(" (%s)", tag) else ""),
         subtitle = sprintf("r = %.2f, R^2 = %.2f", r, r^2),
         x = "log p_j", y = "Posterior median psi_j") +
    theme_minimal(base_size = 12)
  if (!is.null(save_path)) ggsave(save_path, p, width = 7, height = 5, dpi = 150)
  invisible(p)
}

plot_class_means <- function(class_tbl, save_path = NULL, tag = "") {
  p <- ggplot(class_tbl, aes(reorder(class, mu_median), mu_median)) +
    geom_pointrange(aes(ymin = mu_lo, ymax = mu_hi)) +
    coord_flip() +
    labs(title = sprintf("Class thresholds%s",
                         if (nzchar(tag)) sprintf(" (%s)", tag) else ""),
         x = NULL, y = "Log-threshold (log-tokens)") +
    theme_minimal(base_size = 12)
  if (!is.null(save_path)) ggsave(save_path, p, width = 6, height = 3.5, dpi = 150)
  invisible(p)
}

plot_posterior_density <- function(fit, param, save_path = NULL,
                                   xlab = param, tag = "",
                                   vline = NULL, xlim = NULL) {
  draws <- as_draws_df(fit)
  p <- ggplot(tibble(x = draws[[param]]), aes(x)) +
    geom_density(fill = "steelblue", alpha = 0.4) +
    labs(title = sprintf("%s%s", param,
                         if (nzchar(tag)) sprintf(" (%s)", tag) else ""),
         x = xlab, y = "density") +
    theme_minimal(base_size = 12)
  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linetype = "dashed")
  if (!is.null(xlim)) p <- p + xlim(xlim)
  if (!is.null(save_path)) ggsave(save_path, p, width = 6, height = 4, dpi = 150)
  invisible(p)
}
