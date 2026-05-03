## Schematic plots for explaining the model and its properties.
##
## NOT FITTED — these are deterministic generations from chosen "true"
## parameters under the full model and four ablations. Useful for
## explaining the model in talks / paper figures, and for showing how
## removing each ingredient changes predicted patterns.
##
## Outputs (under outputs/figs/schematic/):
##   A_growth_curves.png   - vocab(age) families, 5 model variants
##   B_per_word_growth.png - per-word growth curves, full vs 1PL
##   C_psi_vs_freq.png     - item difficulty vs frequency, freq-driven vs random
##
## Run:  Rscript model/scripts/schematic_plots.R

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT_DIR <- file.path(PATHS$figs_dir, "schematic")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Choose schematic "true" parameters ----------------------------
set.seed(2026)
N_KIDS  <- 600          # for stable percentile bands
N_ITEMS <- 200
CLASSES <- c("nouns", "verbs", "adjectives", "function_words")
N_CLASS <- length(CLASSES)
AGE_GRID <- seq(8, 36, by = 1)
PERCENTILES <- c(0.05, 0.25, 0.50, 0.75, 0.95)

# Parameters chosen so the FULL MODEL panel roughly matches the
# Wordbank median trajectory (vocab ~150 by age 30). Ablations then
# show as dramatic departures from this baseline.
PARAMS <- list(
  mu_r        = 7.34,
  sigma_r     = 0.5,
  sigma_alpha = 1.5,             # boosted from 1.0 → wider population spread
  sigma_zeta  = 0.6,
  rho_xi_zeta = 0.5,
  s           = 4,
  delta       = 2.5,             # boosted from 1.5 → steeper population growth
  log_H       = log(365),
  a0          = 20,
  # Class-level mean ψ shifted down by 0.5 → words a touch easier
  mu_psi_by_class = c(nouns = 5.0, verbs = 6.0, adjectives = 6.5, function_words = 8.5),
  tau_psi   = 0.7,
  sigma_lambda = 0.35,
  beta_p    = 1.0
)

# ---- 2. Helper: simulate a population given a parameter set -----------
simulate_world <- function(params,
                           ablation = c("full", "rasch_1pl", "no_accel",
                                        "random_psi", "no_start", "no_slopes"),
                           n_kids = N_KIDS, n_items = N_ITEMS,
                           seed = 42) {
  ablation <- match.arg(ablation)
  set.seed(seed)

  # Item draws --------------------------------------------------
  cc <- sample(rep(seq_len(N_CLASS), length.out = n_items))
  log_p <- rnorm(n_items, mean = 0, sd = 1) - 8           # plausible Zipfian-ish
  log_p <- pmin(log_p, -3)                                # cap at -3
  # Class-specific log_p shift: function words higher freq, nouns moderate
  shift <- c(nouns = 0, verbs = -0.5, adjectives = -1, function_words = 1.5)
  log_p <- log_p + shift[cc]

  if (ablation == "random_psi") {
    psi <- rnorm(n_items, mean = mean(params$mu_psi_by_class), sd = params$tau_psi)
  } else {
    psi <- params$mu_psi_by_class[cc] + rnorm(n_items, 0, params$tau_psi)
  }

  # 2PL discrimination
  lambda <- if (ablation == "rasch_1pl") {
    rep(1, n_items)
  } else {
    exp(rnorm(n_items, 0, params$sigma_lambda))
  }

  # Children ----------------------------------------------------
  # bivariate (xi_centered, zeta) with rho; zero zeta if no_slopes ablation
  sz <- if (ablation == "no_slopes") 0 else params$sigma_zeta
  Sigma <- matrix(c(params$sigma_alpha^2 + params$sigma_r^2,
                    params$rho_xi_zeta * sqrt(params$sigma_alpha^2 + params$sigma_r^2) * sz,
                    params$rho_xi_zeta * sqrt(params$sigma_alpha^2 + params$sigma_r^2) * sz,
                    sz^2 + 1e-12), 2, 2)
  L <- chol(Sigma)
  z <- matrix(rnorm(2 * n_kids), 2, n_kids)
  effs <- t(L) %*% z
  xi <- params$mu_r + effs[1, ]
  zeta <- effs[2, ]

  # Effective parameters under each ablation
  delta_eff <- if (ablation == "no_accel") 0   else params$delta
  s_eff     <- if (ablation == "no_start") 0   else params$s
  beta_p    <- params$beta_p

  list(cc = cc, log_p = log_p, psi = psi, lambda = lambda,
       xi = xi, zeta = zeta,
       delta = delta_eff, s = s_eff, beta_p = beta_p,
       ablation = ablation,
       params = params, n_kids = n_kids, n_items = n_items)
}

# Compute predicted P(y=1) for one child × one item × one age vector
predict_p <- function(world, child_i, item_j, ages) {
  ae <- pmax(ages - world$s, 0.01)
  log_age <- log(ae / world$params$a0)
  base <- world$xi[child_i] + world$beta_p * world$log_p[item_j] +
          world$params$log_H +
          (1 + world$delta + world$zeta[child_i]) * log_age -
          world$psi[item_j]
  plogis(world$lambda[item_j] * base)
}

# Predicted vocabulary (sum of P over items) for one child across ages
predict_vocab <- function(world, child_i, ages) {
  ae <- pmax(ages - world$s, 0.01)
  log_age <- log(ae / world$params$a0)
  log_age_mat <- matrix(log_age, nrow = world$n_items, ncol = length(ages),
                        byrow = TRUE)
  zeta_term <- (1 + world$delta + world$zeta[child_i]) * log_age_mat
  base <- world$xi[child_i] + world$beta_p * world$log_p +
          world$params$log_H + zeta_term - world$psi
  eta <- world$lambda * base
  probs <- plogis(eta)
  colSums(probs)
}

# ---- 3. Figure A: percentile bands of vocab growth ---------------------
##
## For each variant, simulate N_KIDS kids' trajectories and at each age
## compute population percentiles. Show as ribbons + median lines.
## Add a "Wordbank (real data)" panel with empirical percentiles for
## direct visual comparison.
fig_A <- function() {
  variants <- c("full", "rasch_1pl", "no_accel", "no_slopes",
                "random_psi", "no_start")
  variant_labels <- c(
    full        = "Full model",
    rasch_1pl   = "1PL  (uniform discrimination)",
    no_accel    = "No acceleration  (δ = 0)",
    no_slopes   = "No per-child slopes  (σ_ζ = 0)",
    random_psi  = "Difficulty independent of frequency",
    no_start    = "No start time  (s = 0)"
  )

  # Schematic variants
  rows <- list()
  for (v in variants) {
    world <- simulate_world(PARAMS, ablation = v)
    vocab_mat <- sapply(seq_len(world$n_kids),
                        function(i) predict_vocab(world, i, AGE_GRID))
    # vocab_mat is age x kids; compute percentiles per age
    pct <- apply(vocab_mat, 1, quantile, probs = PERCENTILES)
    rows[[length(rows) + 1]] <- data.frame(
      panel = variant_labels[v],
      age = AGE_GRID,
      p05 = pct["5%", ],  p25 = pct["25%", ], p50 = pct["50%", ],
      p75 = pct["75%", ], p95 = pct["95%", ]
    )
  }
  df_sim <- do.call(rbind, rows)
  df_sim$panel <- factor(df_sim$panel,
                         levels = c(variant_labels[variants],
                                    "Wordbank (observed)"))

  # Real data: pool committed Wordbank cross-sectional + longitudinal
  wb <- tryCatch({
    sd_  <- readRDS(file.path(PATHS$fits_dir, "subset_data.rds"))$stan_data
    long <- readRDS(file.path(PATHS$fits_dir, "long_subset_data.rds"))$stan_data
    csec <- data.frame(age = sd_$age[sd_$ii],
                       jj  = sd_$jj, y = sd_$y) |>
      dplyr::group_by(child = sd_$ii) |>
      dplyr::summarise(age = first(age), vocab = sum(y), .groups = "drop")
    long_v <- data.frame(age = long$admin_age[long$aa],
                         jj  = long$jj, y = long$y) |>
      dplyr::group_by(adm = long$aa) |>
      dplyr::summarise(age = first(age), vocab = sum(y), .groups = "drop")
    bind_rows(csec, long_v)
  }, error = function(e) NULL)

  if (!is.null(wb) && nrow(wb) > 0) {
    wb <- wb |> mutate(age_bin = round(age))
    wb_pct <- wb |>
      dplyr::group_by(age_bin) |>
      dplyr::filter(dplyr::n() >= 5) |>
      dplyr::summarise(
        p05 = quantile(vocab, 0.05),
        p25 = quantile(vocab, 0.25),
        p50 = quantile(vocab, 0.50),
        p75 = quantile(vocab, 0.75),
        p95 = quantile(vocab, 0.95),
        .groups = "drop"
      ) |>
      mutate(age = age_bin, panel = "Wordbank (observed)") |>
      dplyr::select(panel, age, p05, p25, p50, p75, p95)
    df_sim <- bind_rows(df_sim, wb_pct)
  }

  ggplot(df_sim, aes(age)) +
    geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.20, fill = "steelblue") +
    geom_ribbon(aes(ymin = p25, ymax = p75), alpha = 0.30, fill = "steelblue") +
    geom_line(aes(y = p50), colour = "navy", linewidth = 0.7) +
    facet_wrap(~panel, ncol = 4) +
    coord_cartesian(ylim = c(0, N_ITEMS)) +
    labs(title = "Vocabulary growth: model variants and observed data",
         subtitle = sprintf("Bands = 5/25/75/95th percentiles, line = median. Schematic kids: N = %d per panel.",
                            N_KIDS),
         x = "Age (months)", y = "Vocabulary (count of items)") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
}

# ---- 4. Figure B: per-word growth, SINGLE CHILD, 2PL vs 1PL ------------
##
## Single child (median ξ, ζ = 0): for each of 8 items, plot P(produces)
## across age under full (2PL) and 1PL ablation. With one fixed child,
## the per-word slope difference is much sharper than under population
## averaging.
fig_B <- function() {
  world_full <- simulate_world(PARAMS, ablation = "full",      seed = 99)
  world_1pl  <- simulate_world(PARAMS, ablation = "rasch_1pl", seed = 99)

  # Replace child 1 with a "median" hypothetical child for both worlds
  median_xi <- median(world_full$xi)
  world_full$xi[1] <- median_xi; world_full$zeta[1] <- 0
  world_1pl$xi[1]  <- median_xi; world_1pl$zeta[1]  <- 0

  # Pick 8 items: easiest + hardest of each class, by ψ in the full world
  set.seed(2)
  picks <- tibble::tibble(
    j = seq_len(world_full$n_items),
    class = CLASSES[world_full$cc],
    log_p = world_full$log_p,
    psi   = world_full$psi,
    lambda = world_full$lambda
  ) |>
    dplyr::group_by(class) |>
    dplyr::arrange(psi) |>
    dplyr::slice(c(1, dplyr::n())) |>     # easiest + hardest in each class
    dplyr::mutate(extreme = c("easy", "hard")) |>
    dplyr::ungroup() |>
    dplyr::mutate(label = sprintf("%s — %s\n(λ = %.2f)", class, extreme, lambda))

  rows <- list()
  for (k in seq_len(nrow(picks))) {
    j <- picks$j[k]
    p_full <- sapply(AGE_GRID, function(a) predict_p(world_full, 1, j, a))
    p_1pl  <- sapply(AGE_GRID, function(a) predict_p(world_1pl,  1, j, a))
    rows[[length(rows) + 1]] <- data.frame(
      label = picks$label[k], age = AGE_GRID,
      p_full = p_full, p_1pl = p_1pl
    )
  }
  df <- do.call(rbind, rows) |>
    tidyr::pivot_longer(c(p_full, p_1pl), names_to = "model", values_to = "p") |>
    dplyr::mutate(model = ifelse(model == "p_full", "Full (2PL)", "1PL"))

  ggplot(df, aes(age, p, colour = model, linetype = model)) +
    geom_line(linewidth = 0.85) +
    facet_wrap(~label, ncol = 4) +
    scale_colour_manual(values = c("Full (2PL)" = "firebrick", "1PL" = "steelblue")) +
    scale_linetype_manual(values = c("Full (2PL)" = "solid", "1PL" = "dashed")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Per-word growth, single average child, 2PL vs 1PL",
         subtitle = "Single hypothetical child (median ξ, ζ = 0). 2PL gives word-specific slopes via λⱼ; 1PL forces all words to share the population slope.",
         x = "Age (months)", y = "P(produces)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
}

# ---- 5. Figure C: ψ_j vs log p_j, frequency-driven vs random -----------
fig_C <- function() {
  world_full <- simulate_world(PARAMS, ablation = "full",       seed = 7)
  world_rand <- simulate_world(PARAMS, ablation = "random_psi", seed = 7)
  # Same log_p across both (set.seed in simulate_world ensures this)

  df <- bind_rows(
    tibble::tibble(world = "Frequency-driven ψ (full model)",
                   class = CLASSES[world_full$cc],
                   log_p = world_full$log_p, psi = world_full$psi),
    tibble::tibble(world = "Random ψ (ablation)",
                   class = CLASSES[world_rand$cc],
                   log_p = world_rand$log_p, psi = world_rand$psi)
  )

  ggplot(df, aes(log_p, psi, colour = class)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                colour = "grey30", linewidth = 0.5,
                linetype = "dashed") +
    facet_wrap(~world, ncol = 2) +
    labs(title = "Item difficulty vs frequency",
         subtitle = "Left: ψ structured by class with class-typical means (full model). Right: ψ drawn independently of frequency (ablation).",
         x = "log frequency  log p_j",
         y = "log threshold  ψ_j",
         colour = "class") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

# ---- 6. Render & save -------------------------------------------------
pA <- fig_A()
ggsave(file.path(OUT_DIR, "A_growth_curves.png"), pA,
       width = 11, height = 7, dpi = 150)

pB <- fig_B()
ggsave(file.path(OUT_DIR, "B_per_word_growth.png"), pB,
       width = 11, height = 6, dpi = 150)

pC <- fig_C()
ggsave(file.path(OUT_DIR, "C_psi_vs_freq.png"), pC,
       width = 11, height = 5, dpi = 150)

cat("Saved schematic figures to", OUT_DIR, "\n")
