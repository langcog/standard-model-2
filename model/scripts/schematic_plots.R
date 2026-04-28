## Schematic plots for explaining the model and its properties.
##
## NOT FITTED — these are deterministic generations from chosen "true"
## parameters under the full model and four ablations. Useful for
## explaining the model in talks / paper figures, and for showing how
## removing each ingredient changes predicted patterns.
##
## Outputs (under model/figs/schematic/):
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
N_KIDS  <- 30           # for the population trajectories
N_ITEMS <- 200
CLASSES <- c("nouns", "verbs", "adjectives", "function_words")
N_CLASS <- length(CLASSES)
AGE_GRID <- seq(8, 36, by = 0.5)

# Parameters chosen for visual clarity, broadly inspired by fit values
PARAMS <- list(
  mu_r        = 7.34,            # log tokens/hr
  sigma_r     = 0.5,
  sigma_alpha = 1.0,             # less extreme than fit (~2) for clarity
  sigma_zeta  = 0.5,
  rho_xi_zeta = 0.5,
  s           = 4,
  delta       = 1.5,
  log_H       = log(365),
  a0          = 20,
  # Class-level mean ψ: nouns easiest, function hardest
  mu_psi_by_class = c(nouns = 5.5, verbs = 6.5, adjectives = 7.0, function_words = 9),
  tau_psi   = 0.7,
  sigma_lambda = 0.35,           # 2PL discrimination spread
  beta_p    = 1.0                # frequency coefficient (1 = pure accumulator)
)

# ---- 2. Helper: simulate a population given a parameter set -----------
simulate_world <- function(params,
                           ablation = c("full", "rasch_1pl", "no_accel",
                                        "random_psi", "no_start"),
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
  # bivariate (xi_centered, zeta) with rho
  Sigma <- matrix(c(params$sigma_alpha^2 + params$sigma_r^2,
                    params$rho_xi_zeta * sqrt(params$sigma_alpha^2 + params$sigma_r^2) * params$sigma_zeta,
                    params$rho_xi_zeta * sqrt(params$sigma_alpha^2 + params$sigma_r^2) * params$sigma_zeta,
                    params$sigma_zeta^2), 2, 2)
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

# ---- 3. Figure A: family of vocab growth curves ------------------------
fig_A <- function() {
  variants <- c("full", "rasch_1pl", "no_accel", "random_psi", "no_start")
  variant_labels <- c(
    full        = "Full model",
    rasch_1pl   = "1PL ablation (uniform discrimination)",
    no_accel    = "No acceleration  (δ = 0)",
    random_psi  = "Difficulty independent of frequency",
    no_start    = "No start time  (s = 0)"
  )

  rows <- list()
  for (v in variants) {
    world <- simulate_world(PARAMS, ablation = v)
    for (i in seq_len(world$n_kids)) {
      vocab <- predict_vocab(world, i, AGE_GRID)
      rows[[length(rows) + 1]] <- data.frame(
        variant = v, kid = i, age = AGE_GRID, vocab = vocab,
        xi_rank = rank(world$xi)[i] / world$n_kids
      )
    }
  }
  df <- do.call(rbind, rows)
  df$variant <- factor(df$variant, levels = variants,
                       labels = variant_labels[variants])

  ggplot(df, aes(age, vocab, group = kid, colour = xi_rank)) +
    geom_line(alpha = 0.65, linewidth = 0.4) +
    facet_wrap(~variant, ncol = 3) +
    scale_colour_viridis_c(option = "plasma", end = 0.9, name = "child rank\n(by ξ)") +
    coord_cartesian(ylim = c(0, N_ITEMS)) +
    labs(title = "Schematic vocabulary trajectories under the model and four ablations",
         subtitle = sprintf("%d simulated children × %d items, ages 8–36 mo. Each line is one child.",
                            N_KIDS, N_ITEMS),
         x = "Age (months)", y = "Predicted vocabulary (count of items)") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
}

# ---- 4. Figure B: per-word growth curves, full vs 1PL ------------------
fig_B <- function() {
  world_full <- simulate_world(PARAMS, ablation = "full",  seed = 99)
  world_1pl  <- simulate_world(PARAMS, ablation = "rasch_1pl", seed = 99)

  # Pick 8 items across class & log_p quartile
  set.seed(2)
  picks <- tibble::tibble(
    j = seq_len(world_full$n_items),
    class = CLASSES[world_full$cc],
    log_p = world_full$log_p,
    psi   = world_full$psi
  ) |>
    dplyr::group_by(class) |>
    dplyr::arrange(psi) |>
    dplyr::slice(c(1, dplyr::n())) |>     # easiest + hardest of each class
    dplyr::ungroup()

  # Population-average P(y=1) at each age = average over kids
  avg_curve <- function(world, j) {
    sapply(AGE_GRID, function(a) {
      ages <- rep(a, world$n_kids)
      mean(sapply(seq_len(world$n_kids),
                  function(i) predict_p(world, i, j, a)))
    })
  }

  rows <- list()
  for (k in seq_len(nrow(picks))) {
    j <- picks$j[k]
    rows[[length(rows) + 1]] <- data.frame(
      class = picks$class[k], item = sprintf("%s\n(%s)",
                                             picks$class[k],
                                             ifelse(k %% 2 == 1, "easy", "hard")),
      age = AGE_GRID,
      p_full = avg_curve(world_full, j),
      p_1pl  = avg_curve(world_1pl, j)
    )
  }
  df <- do.call(rbind, rows)
  df_long <- df |>
    tidyr::pivot_longer(c(p_full, p_1pl), names_to = "model",
                        values_to = "p") |>
    dplyr::mutate(model = ifelse(model == "p_full", "Full (2PL)", "1PL"))

  ggplot(df_long, aes(age, p, colour = model, linetype = model)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~item, ncol = 4) +
    scale_colour_manual(values = c("Full (2PL)" = "firebrick", "1PL" = "steelblue")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Per-word population growth curves, 2PL vs 1PL",
         subtitle = "8 items: easiest + hardest of each class. 2PL gives word-specific slopes; 1PL forces them equal.",
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
