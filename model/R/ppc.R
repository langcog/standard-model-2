## Posterior-predictive-check helpers.
##
## Exposes four functions, each taking (fit, bundle, ...) and returning a
## ggplot:
##   ppc_word_growth(fit, bundle, n_words = 12, n_draws = 100)
##   ppc_vocab_trajectory(fit, bundle, n_draws = 100)
##   ppc_vocab_distribution(fit, bundle, age_bins = c(18, 24, 30))
##   ppc_calibration(fit, bundle, n_bins = 20)
##
## Plus ppc_suite(fit, bundle, ...) which combines all four via patchwork
## and saves to a single PNG.
##
## All predictions assume the Stan model log_irt.stan (xi_i is a transformed
## parameter; psi_j, s, delta, sigma_xi, log_H, a0 accessible either from
## the fit or from MODEL_CONSTANTS + bundle).

suppressPackageStartupMessages({
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(posterior)
  library(tibble)
})

# ----------------------------------------------------------------------------
# Core prediction helpers
# ----------------------------------------------------------------------------

## Thin posterior draws to a manageable number for plotting.
thin_posterior <- function(fit, n_draws = 100) {
  d <- as_draws_df(fit)
  idx <- round(seq(1, nrow(d), length.out = n_draws))
  d[idx, , drop = FALSE]
}

## Compute eta_ij for specific (child ability xi, word psi, age, log_p)
## given a single draw of (s, delta) and optional lambda (2PL discrimination)
## and optional zeta (per-child slope deviation).
eta_from_parts <- function(xi, log_p, psi, age, s, delta,
                           log_H, a0, lambda = 1, zeta = 0) {
  ae <- pmax(age - s, 0.01)
  base <- xi + log_p + log_H + (1 + delta + zeta) * log(ae / a0) - psi
  lambda * base
}

## Extract lambda_j per draw. Returns a matrix [n_draws x J].
## Falls back to all-1 if the fit didn't include lambda (older fits).
extract_lambda_mat <- function(draws, J) {
  lam_cols <- grep("^lambda\\[", names(draws), value = TRUE)
  if (length(lam_cols) == 0) {
    matrix(1, nrow = nrow(draws), ncol = J)
  } else {
    as.matrix(draws[, lam_cols, drop = FALSE])
  }
}

## Extract zeta_i per draw. Returns a matrix [n_draws x I].
## Falls back to all-0 if the fit didn't include zeta (older fits).
extract_zeta_mat <- function(draws, I) {
  z_cols <- grep("^zeta\\[", names(draws), value = TRUE)
  if (length(z_cols) == 0) {
    matrix(0, nrow = nrow(draws), ncol = I)
  } else {
    as.matrix(draws[, z_cols, drop = FALSE])
  }
}

## Population-average P(y=1) at (age, word j) by MC integration over the
## posterior-inferred xi distribution (N(mu_r, sigma_xi^2)).
## Returns a matrix [n_draws x n_ages] of marginal P.
predict_word_curve_marginal <- function(draws, word_j, log_p_j, ages,
                                        log_H, a0, mu_r,
                                        n_xi_samples = 300) {
  psi_col <- sprintf("psi[%d]", word_j)
  psi_draws      <- draws[[psi_col]]
  s_draws        <- draws$s
  delta_draws    <- draws$delta
  sigma_xi_draws <- draws$sigma_xi

  lam_col <- sprintf("lambda[%d]", word_j)
  lambda_draws <- if (lam_col %in% names(draws)) draws[[lam_col]] else rep(1, length(psi_draws))

  # sigma_zeta may not exist in older fits; default to 0
  sig_z_draws <- if ("sigma_zeta" %in% names(draws)) draws$sigma_zeta else rep(0, length(psi_draws))

  out <- matrix(NA_real_, nrow = length(psi_draws), ncol = length(ages))
  for (d in seq_along(psi_draws)) {
    xi_samp   <- rnorm(n_xi_samples, mu_r, sigma_xi_draws[d])
    zeta_samp <- rnorm(n_xi_samples, 0, sig_z_draws[d])
    for (a_idx in seq_along(ages)) {
      eta <- eta_from_parts(xi_samp, log_p_j, psi_draws[d],
                            ages[a_idx], s_draws[d], delta_draws[d],
                            log_H, a0,
                            lambda = lambda_draws[d],
                            zeta = zeta_samp)
      out[d, a_idx] <- mean(plogis(eta))
    }
  }
  out
}

## For each posterior draw, compute predicted P(y=1) for every observed
## (child, word) cell conditional on that child's xi_i draw.
## Returns a matrix [n_draws x N]. Uses lambda_j if present (2PL fit).
predict_cellwise_cond <- function(draws, bundle, log_H, a0) {
  sd_ <- bundle$stan_data
  xi_cols  <- grep("^xi\\[",  names(draws), value = TRUE)
  psi_cols <- grep("^psi\\[", names(draws), value = TRUE)

  xi_mat   <- as.matrix(draws[, xi_cols,  drop = FALSE])   # n_draws x I
  psi_mat  <- as.matrix(draws[, psi_cols, drop = FALSE])   # n_draws x J
  lam_mat  <- extract_lambda_mat(draws, J = ncol(psi_mat)) # n_draws x J
  zeta_mat <- extract_zeta_mat(draws,  I = ncol(xi_mat))   # n_draws x I
  s_d      <- draws$s
  delta_d  <- draws$delta

  ii <- sd_$ii; jj <- sd_$jj
  log_p_obs <- sd_$log_p[jj]
  age_obs   <- sd_$age[ii]

  n_draws <- nrow(draws); N <- length(ii)
  out <- matrix(NA_real_, nrow = n_draws, ncol = N)
  for (d in seq_len(n_draws)) {
    xi_per_obs   <- xi_mat[d,   ii]
    psi_per_obs  <- psi_mat[d,  jj]
    lam_per_obs  <- lam_mat[d,  jj]
    zeta_per_obs <- zeta_mat[d, ii]
    ae <- pmax(age_obs - s_d[d], 0.01)
    base <- xi_per_obs + log_p_obs + log_H +
            (1 + delta_d[d] + zeta_per_obs) * log(ae / a0) - psi_per_obs
    out[d, ] <- plogis(lam_per_obs * base)
  }
  out
}

# ----------------------------------------------------------------------------
# Plot 1 - per-word growth curves
# ----------------------------------------------------------------------------

## Pick a representative set of words: n/class_count per class, chosen at
## evenly spaced psi quantiles within class.
select_words_for_plot <- function(bundle, fit, n_words = 12) {
  draws    <- as_draws_df(fit)
  psi_cols <- grep("^psi\\[", names(draws), value = TRUE)
  psi_med  <- sapply(psi_cols, function(p) median(draws[[p]]))

  info <- bundle$word_info
  info$psi_med <- psi_med
  info$class   <- bundle$class_levels[info$cc]

  per_class <- max(1, floor(n_words / length(unique(info$class))))
  qs <- seq(0.1, 0.9, length.out = per_class)

  info %>% group_by(class) %>%
    arrange(psi_med) %>%
    slice(round(qs * (n() - 1)) + 1) %>%
    ungroup() %>%
    arrange(class, psi_med)
}

ppc_word_growth <- function(fit, bundle,
                            n_words = 12,
                            n_draws = 100,
                            constants = MODEL_CONSTANTS) {
  draws <- thin_posterior(fit, n_draws)

  picks <- select_words_for_plot(bundle, fit, n_words)
  ages <- seq(16, 30, by = 0.5)

  sd_ <- bundle$stan_data

  # Observed proportions by age for each picked word
  obs_tbl <- tibble(
    person_age = sd_$age[sd_$ii],
    jj = sd_$jj,
    y  = sd_$y
  ) %>%
    filter(jj %in% picks$jj) %>%
    mutate(age_bin = round(person_age)) %>%
    group_by(jj, age_bin) %>%
    summarise(n = n(), p = mean(y), .groups = "drop") %>%
    left_join(picks %>% select(jj, item, class), by = "jj") %>%
    mutate(label = sprintf("%s (%s)", item, class))

  # Predicted curves
  pred_tbl <- purrr::map_dfr(seq_len(nrow(picks)), function(k) {
    row <- picks[k, ]
    M <- predict_word_curve_marginal(
      draws, word_j = row$jj, log_p_j = log(row$prob), ages = ages,
      log_H = constants$log_H, a0 = constants$a0,
      mu_r = sd_$mu_r, n_xi_samples = 200)
    tibble(jj = row$jj, item = row$item, class = row$class,
           age = ages,
           med = apply(M, 2, median),
           lo  = apply(M, 2, quantile, .025),
           hi  = apply(M, 2, quantile, .975))
  }) %>%
    mutate(label = sprintf("%s (%s)", item, class))

  # Order facets by posterior-median psi
  facet_order <- picks %>% mutate(label = sprintf("%s (%s)", item,
    bundle$class_levels[cc])) %>% pull(label)

  ggplot(pred_tbl, aes(age)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = class), alpha = 0.2) +
    geom_line(aes(y = med, colour = class), linewidth = 0.6) +
    geom_point(data = obs_tbl, aes(age_bin, p, size = n),
               alpha = 0.6, show.legend = FALSE) +
    facet_wrap(~factor(label, levels = facet_order), ncol = 4) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_size_area(max_size = 2.4) +
    labs(title = "PPC: per-word growth curves",
         subtitle = "Dots = observed proportion; line/ribbon = posterior predictive (pop avg)",
         x = "Age (months)", y = "P(produces)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

# ----------------------------------------------------------------------------
# Plot 2 - aggregate vocabulary trajectory
# ----------------------------------------------------------------------------

ppc_vocab_trajectory <- function(fit, bundle,
                                  n_draws = 100,
                                  constants = MODEL_CONSTANTS) {
  draws <- thin_posterior(fit, n_draws)
  sd_   <- bundle$stan_data

  # Observed: total vocab per child
  obs_df <- tibble(ii = sd_$ii, y = sd_$y) %>%
    group_by(ii) %>%
    summarise(vocab = sum(y), .groups = "drop") %>%
    mutate(age = sd_$age[ii])

  # Predicted: for each draw, sum P(y=1) across items per child.
  # Each draw yields a predicted vocab per child.
  P <- predict_cellwise_cond(draws, bundle,
                             log_H = constants$log_H, a0 = constants$a0)
  # Aggregate per child
  I <- sd_$I
  pred_vocab <- matrix(NA_real_, nrow = nrow(draws), ncol = I)
  for (i in seq_len(I)) {
    cols <- which(sd_$ii == i)
    pred_vocab[, i] <- rowSums(P[, cols, drop = FALSE])
  }
  pred_df <- tibble(
    ii = seq_len(I),
    age = sd_$age,
    med = apply(pred_vocab, 2, median),
    lo  = apply(pred_vocab, 2, quantile, .025),
    hi  = apply(pred_vocab, 2, quantile, .975)
  )

  ggplot(pred_df, aes(age)) +
    geom_point(data = obs_df, aes(age, vocab), alpha = 0.25, size = 1.2) +
    geom_linerange(aes(ymin = lo, ymax = hi), alpha = 0.25, colour = "firebrick") +
    geom_point(aes(y = med), colour = "firebrick", size = 1.2) +
    labs(title = "PPC: child vocabulary vs. age",
         subtitle = "Grey = observed total vocab; red = posterior predicted (per child, 95% CI)",
         x = "Age (months)", y = sprintf("Words produced (of %d)", sd_$J)) +
    theme_minimal(base_size = 11)
}

# ----------------------------------------------------------------------------
# Plot 3 - vocabulary distribution at age slices
# ----------------------------------------------------------------------------

## Vocab distribution at age slices.
## Computes the MARGINAL posterior-predictive: sample hypothetical children
## from xi ~ N(mu_r, sigma_xi^2) under each posterior draw, predict their
## vocab, compare distribution to observed at each age slice. This is the
## true generative test of the population model, unlike a conditional PPC
## using each observed child's own posterior xi_i.
ppc_vocab_distribution <- function(fit, bundle,
                                    age_bins = c(18, 24, 30),
                                    half_window = 2,
                                    n_draws = 100,
                                    n_synthetic = 400,
                                    constants = MODEL_CONSTANTS) {
  draws <- thin_posterior(fit, n_draws)
  sd_   <- bundle$stan_data
  psi_cols <- grep("^psi\\[", names(draws), value = TRUE)
  psi_mat  <- as.matrix(draws[, psi_cols, drop = FALSE])  # n_draws x J
  lam_mat  <- extract_lambda_mat(draws, J = ncol(psi_mat))
  s_d      <- draws$s
  delta_d  <- draws$delta
  sig_xi_d <- draws$sigma_xi
  sig_z_d  <- if ("sigma_zeta" %in% names(draws)) draws$sigma_zeta else rep(0, nrow(draws))
  mu_r     <- sd_$mu_r
  log_p    <- sd_$log_p

  # Observed vocab by child
  obs_df <- tibble(ii = sd_$ii, y = sd_$y) %>%
    group_by(ii) %>%
    summarise(vocab = sum(y), .groups = "drop") %>%
    mutate(age = sd_$age[ii])

  # For each age slice, simulate n_synthetic hypothetical children per draw
  # and record their predicted total vocab.
  panels <- purrr::map_dfr(age_bins, function(a_c) {
    in_band <- which(abs(sd_$age - a_c) <= half_window)
    if (length(in_band) == 0) return(NULL)
    obs_here <- obs_df$vocab[in_band]

    pred <- numeric(0)
    for (d in seq_len(nrow(draws))) {
      xi_sim   <- rnorm(n_synthetic, mu_r, sig_xi_d[d])
      zeta_sim <- rnorm(n_synthetic, 0, sig_z_d[d])
      ae <- max(a_c - s_d[d], 0.01)
      logL <- log(ae / constants$a0)
      # per-child slope contribution (varies across synthetic kids)
      slope_contrib <- (1 + delta_d[d] + zeta_sim) * logL       # length n_synthetic
      base_word <- log_p + constants$log_H - psi_mat[d, ]        # length J
      # eta_ij = lambda_j * (xi_i + slope_contrib_i + base_word_j)
      child_part <- xi_sim + slope_contrib                       # length n_synthetic
      eta <- outer(child_part, rep(1, length(base_word))) +
             outer(rep(1, n_synthetic), base_word)
      eta <- eta * rep(lam_mat[d, ], each = n_synthetic)
      dim(eta) <- c(n_synthetic, length(base_word))
      probs <- plogis(eta)
      # Sample binary outcomes; sum rows to get per-child vocab
      y_sim <- rbinom(length(probs), 1, probs)
      dim(y_sim) <- dim(probs)
      pred <- c(pred, rowSums(y_sim))
    }
    bind_rows(
      tibble(age_center = a_c, kind = "observed",  vocab = obs_here),
      tibble(age_center = a_c, kind = "predicted (marginal)", vocab = pred)
    )
  })

  ggplot(panels, aes(vocab, fill = kind, colour = kind)) +
    geom_density(alpha = 0.35, adjust = 1.2) +
    facet_wrap(~sprintf("%d ± %d mo", age_center, half_window), ncol = 3,
               scales = "free_y") +
    labs(title = "PPC: marginal vocab distribution at age slices",
         subtitle = paste0("Observed vs. generated vocab for hypothetical children ",
                           "(xi ~ N(mu_r, sigma_xi^2) from posterior)"),
         x = sprintf("Words produced (of %d)", sd_$J), y = "density",
         fill = NULL, colour = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

# ----------------------------------------------------------------------------
# Plot 4 - calibration
# ----------------------------------------------------------------------------

ppc_calibration <- function(fit, bundle, n_bins = 20,
                             n_draws = 100,
                             constants = MODEL_CONSTANTS) {
  draws <- thin_posterior(fit, n_draws)
  sd_   <- bundle$stan_data

  P <- predict_cellwise_cond(draws, bundle,
                             log_H = constants$log_H, a0 = constants$a0)
  pmean <- colMeans(P)

  tibble(p = pmean, y = sd_$y) %>%
    mutate(bin = cut(p, breaks = seq(0, 1, length.out = n_bins + 1),
                     include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(p_mean = mean(p), y_mean = mean(y), n = n(),
              .groups = "drop") %>%
    ggplot(aes(p_mean, y_mean)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey40") +
    geom_point(aes(size = n), alpha = 0.8) +
    scale_size_area(max_size = 6) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "PPC: calibration (binned)",
         subtitle = "Each dot = bin of predicted P(y=1); size = # cells",
         x = "Mean predicted P(y=1)", y = "Observed proportion y=1") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
}

# ----------------------------------------------------------------------------
# Plot 5 - 2PL diagnostic: lambda_j by class and by log_p
# ----------------------------------------------------------------------------

## Plot 6 — Variance-of-ability by age.
## Directly tests the per-child-slopes claim: if sigma_zeta > 0, the
## variance of effective ability theta_i = xi_i + (1+delta+zeta_i) * L(a_i)
## across children grows with |L(a)|.
## Compares observed within-age vocab variance (on logit scale) to model
## predicted variance, for both the per-child conditional posterior and
## the model-implied marginal distribution.
ppc_variance_by_age <- function(fit, bundle,
                                 n_draws = 100,
                                 age_bin_width = 1,
                                 constants = MODEL_CONSTANTS) {
  draws <- thin_posterior(fit, n_draws)
  sd_   <- bundle$stan_data

  # Observed: within each age bin, compute sample SD of logit(vocab/J)
  obs_df <- tibble(ii = sd_$ii, y = sd_$y) %>%
    group_by(ii) %>%
    summarise(vocab = sum(y), .groups = "drop") %>%
    mutate(age = sd_$age[ii],
           age_bin = round(age / age_bin_width) * age_bin_width,
           logit_v = qlogis(pmin(pmax((vocab + 0.5) / (sd_$J + 1), 1e-4), 1 - 1e-4))) %>%
    group_by(age_bin) %>%
    summarise(n = n(), sd_logit_v = sd(logit_v), .groups = "drop") %>%
    filter(n >= 5)

  # Model-implied: theta_i = xi_i + (1+delta+zeta_i) * L(a)
  # Var(theta | a) = sigma_xi^2 + sigma_zeta^2 * L(a)^2
  s_d      <- draws$s
  sig_xi_d <- draws$sigma_xi
  sig_z_d  <- if ("sigma_zeta" %in% names(draws)) draws$sigma_zeta else rep(0, nrow(draws))

  ages_seq <- seq(min(sd_$age), max(sd_$age), by = 0.25)
  pred_df <- purrr::map_dfr(seq_along(ages_seq), function(k) {
    a <- ages_seq[k]
    logL <- log(pmax(a - s_d, 0.01) / constants$a0)
    var_theta <- sig_xi_d^2 + sig_z_d^2 * logL^2
    tibble(age = a,
           sd_med = median(sqrt(var_theta)),
           sd_lo  = quantile(sqrt(var_theta), 0.025),
           sd_hi  = quantile(sqrt(var_theta), 0.975))
  })

  ggplot(pred_df, aes(age)) +
    geom_ribbon(aes(ymin = sd_lo, ymax = sd_hi), alpha = 0.25,
                fill = "steelblue") +
    geom_line(aes(y = sd_med), colour = "steelblue", linewidth = 0.8) +
    geom_point(data = obs_df, aes(age_bin, sd_logit_v, size = n),
               alpha = 0.7, show.legend = FALSE) +
    scale_size_area(max_size = 4) +
    labs(title = "PPC: within-age ability SD vs. age",
         subtitle = "Dots = observed SD of logit vocab proportion within 1-mo bins; ribbon = model sqrt(Var(theta|age))",
         x = "Age (months)", y = "SD of ability (logit scale)") +
    theme_minimal(base_size = 11)
}

ppc_lambda_diagnostics <- function(fit, bundle) {
  draws <- as_draws_df(fit)
  lam_cols <- grep("^lambda\\[", names(draws), value = TRUE)
  if (length(lam_cols) == 0) {
    # Old fit without 2PL; return NULL so callers can skip.
    return(NULL)
  }
  lam_med <- sapply(lam_cols, function(p) median(draws[[p]]))
  lam_lo  <- sapply(lam_cols, function(p) quantile(draws[[p]], .025))
  lam_hi  <- sapply(lam_cols, function(p) quantile(draws[[p]], .975))

  info <- bundle$word_info %>%
    mutate(lambda_med = lam_med,
           lambda_lo = lam_lo, lambda_hi = lam_hi,
           class = bundle$class_levels[cc],
           log_p = log(prob))

  # Panel A: lambda_j by class (box + points)
  pA <- ggplot(info, aes(class, lambda_med, colour = class)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_jitter(width = 0.25, alpha = 0.6, size = 1.3) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA) +
    scale_y_log10() +
    labs(title = "lambda_j by lexical class",
         subtitle = "Dashed line = 1 (no discrimination effect)",
         x = NULL, y = "lambda_j (log scale)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")

  # Panel B: lambda_j vs log_p
  r_lp <- cor(info$log_p, log(info$lambda_med))
  pB <- ggplot(info, aes(log_p, lambda_med, colour = class)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_point(alpha = 0.7, size = 1.3) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey40",
                linewidth = 0.5, aes(group = 1)) +
    scale_y_log10() +
    labs(title = "lambda_j vs. log frequency",
         subtitle = sprintf("r(log lambda, log p) = %.2f", r_lp),
         x = "log p_j (CHILDES)", y = "lambda_j (log scale)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  # Panel C: posterior density of sigma_lambda
  sig_draws <- draws$sigma_lambda
  pC <- ggplot(tibble(x = sig_draws), aes(x)) +
    geom_density(fill = "darkorchid", alpha = 0.4) +
    labs(title = "sigma_lambda posterior",
         subtitle = "SD of log lambda_j across words; ~0 = effectively Rasch",
         x = "sigma_lambda", y = "density") +
    theme_minimal(base_size = 11)

  list(by_class = pA, by_logp = pB, sigma = pC,
       info = info)
}

# ----------------------------------------------------------------------------
# Driver - ppc_suite
# ----------------------------------------------------------------------------

ppc_suite <- function(fit, bundle,
                      save_dir = NULL,
                      variant_tag = "",
                      n_draws = 100,
                      age_bins = c(18, 24, 30),
                      constants = MODEL_CONSTANTS) {
  message("[ppc] word growth curves...")
  p1 <- ppc_word_growth(fit, bundle, n_draws = n_draws,
                        constants = constants)
  message("[ppc] vocabulary trajectory...")
  p2 <- ppc_vocab_trajectory(fit, bundle, n_draws = n_draws,
                              constants = constants)
  message("[ppc] vocabulary distribution...")
  p3 <- ppc_vocab_distribution(fit, bundle, age_bins = age_bins,
                                n_draws = n_draws,
                                constants = constants)
  message("[ppc] calibration...")
  p4 <- ppc_calibration(fit, bundle, n_draws = n_draws,
                         constants = constants)

  message("[ppc] variance-by-age...")
  pV <- ppc_variance_by_age(fit, bundle, n_draws = n_draws,
                             constants = constants)

  # If the fit has 2PL (lambda not pinned), also compute lambda diagnostics.
  lam_plots <- ppc_lambda_diagnostics(fit, bundle)

  if (!is.null(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    tag <- if (nzchar(variant_tag)) sprintf("_%s", variant_tag) else ""
    ggsave(file.path(save_dir, sprintf("ppc%s_1_word_growth.png", tag)),
           p1, width = 11, height = 8, dpi = 150)
    ggsave(file.path(save_dir, sprintf("ppc%s_2_vocab_trajectory.png", tag)),
           p2, width = 7, height = 5, dpi = 150)
    ggsave(file.path(save_dir, sprintf("ppc%s_3_vocab_distribution.png", tag)),
           p3, width = 10, height = 4, dpi = 150)
    ggsave(file.path(save_dir, sprintf("ppc%s_4_calibration.png", tag)),
           p4, width = 5, height = 5, dpi = 150)
    ggsave(file.path(save_dir, sprintf("ppc%s_8_variance_by_age.png", tag)),
           pV, width = 7, height = 5, dpi = 150)
    if (!is.null(lam_plots)) {
      ggsave(file.path(save_dir, sprintf("ppc%s_5_lambda_by_class.png", tag)),
             lam_plots$by_class, width = 6, height = 4, dpi = 150)
      ggsave(file.path(save_dir, sprintf("ppc%s_6_lambda_vs_logp.png", tag)),
             lam_plots$by_logp, width = 7, height = 5, dpi = 150)
      ggsave(file.path(save_dir, sprintf("ppc%s_7_sigma_lambda.png", tag)),
             lam_plots$sigma, width = 6, height = 4, dpi = 150)
    }
    # Combined overview
    if (is.null(lam_plots)) {
      combined <- (p1) / (p2 | p3) / (p4 | pV) +
        plot_annotation(
          title = sprintf("Posterior-predictive checks%s",
                          if (nzchar(variant_tag)) sprintf(" (%s)", variant_tag) else ""),
          theme = theme(plot.title = element_text(face = "bold", size = 14))
        )
      ggsave(file.path(save_dir, sprintf("ppc%s_all.png", tag)),
             combined, width = 14, height = 18, dpi = 120)
    } else {
      combined <- (p1) / (p2 | p3) / (p4 | pV) /
                  (lam_plots$sigma | lam_plots$by_class | lam_plots$by_logp) +
        plot_annotation(
          title = sprintf("Posterior-predictive checks + 2PL diagnostics%s",
                          if (nzchar(variant_tag)) sprintf(" (%s)", variant_tag) else ""),
          theme = theme(plot.title = element_text(face = "bold", size = 14))
        )
      ggsave(file.path(save_dir, sprintf("ppc%s_all.png", tag)),
             combined, width = 14, height = 22, dpi = 120)
    }
  }
  invisible(list(
    panels = list(word = p1, trajectory = p2,
                  distribution = p3, calibration = p4,
                  variance_by_age = pV,
                  lambda_by_class = if (!is.null(lam_plots)) lam_plots$by_class else NULL,
                  lambda_by_logp  = if (!is.null(lam_plots)) lam_plots$by_logp else NULL,
                  sigma_lambda    = if (!is.null(lam_plots)) lam_plots$sigma   else NULL)))
}
