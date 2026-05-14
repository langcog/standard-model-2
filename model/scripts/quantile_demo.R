## 6-panel quantile-fan demo: empirical Wordbank vs. predictions under
## six nested model variants, ALL fit to the same longitudinal bundle
## (fits/long_subset_data.rds = 200 kids x 198 items, 145k obs, 808
## longitudinal admins).
##
## Variants (1-4 reproduce the original 4-panel demo on the longitudinal
## bundle; 5-6 add per-child onset offset s_i):
##   1. pure        - pure accumulator: alpha=0, zeta=0, delta=0
##   2. alpha       - + child efficiency variation (sigma_alpha free)
##   3. kappa       - + age dynamics (sigma_zeta free, delta free)
##   4. slopes      - + both (= M_best, alpha + zeta + delta)
##   5. si_only     - alpha + per-child onset s_i (no zeta)
##   6. slopes_si   - alpha + zeta + s_i (full)
##
## All fits live on Sherlock; locally we have:
##   - fits/summaries/<tag>.draws.rds   (scalar parameter draws, ~2000)
##   - fits/summaries/<tag>_psi.csv     (per-item psi median)
##
## Empirical reference: Wordbank English (American) WS, restricted to
## the bundle's 198 items, computed across all admins via wordbankr or
## the fits/long_items.rds fallback. Same approach as the previous
## 4-panel demo.
##
## Output: outputs/figs/longitudinal/quantile_demo.png

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(wordbankr); library(quantregGrowth); library(readr)
})

OUT_DIR    <- file.path(PATHS$figs_dir, "longitudinal")
SUMMARIES  <- file.path(PATHS$fits_dir, "summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BUNDLE_PATH    <- file.path(PATHS$fits_dir, "long_subset_data.rds")
CACHE_FULL_EMP <- file.path(PATHS$fits_dir, "demo_full_empirical_long.rds")
AGE_RANGE      <- c(16, 30)         # match Wordbank WS empirical range
TAUS           <- c(0.10, 0.25, 0.50, 0.75, 0.90)
N_DRAWS_USE    <- 50                # posterior draws per variant for sim
N_KIDS_PER_AGE <- 80                # simulated kids per age x draw
SEED           <- 2026

## ---- 1. Load bundle ----------------------------------------------
bundle <- readRDS(BUNDLE_PATH)
sd_b   <- bundle$stan_data
cat(sprintf("Bundle: I=%d, J=%d items, A=%d admins, N=%d obs\n",
            sd_b$I, sd_b$J, sd_b$A, sd_b$N))

# Items in this bundle (for empirical filtering)
bundle_items <- bundle$word_info$item
stopifnot(length(bundle_items) == sd_b$J)

## ---- 2. Variant registry -----------------------------------------
VARIANTS <- list(
  pure       = list(label = "1. Pure accumulator (no α, no ζ, no s_i)",
                    tag   = "long_demo_pure"),
  alpha      = list(label = "2. + efficiency variation (α)",
                    tag   = "long_demo_alpha"),
  kappa      = list(label = "3. + scaling variation (ζ → κ)",
                    tag   = "long_demo_kappa"),
  si_only    = list(label = "4. α + per-child onset s_i",
                    tag   = "long_no_freq_si_only"),
  m_best     = list(label = "5. α + ζ (M_best)",
                    tag   = "long_no_freq_slopes"),
  slopes_si  = list(label = "6. α + ζ + s_i (independent)",
                    tag   = "long_no_freq_slopes_si"),
  slopes_si_corr = list(label = "7. α + ζ + s_i (correlated)",
                        tag   = "long_no_freq_slopes_si_corr")
)
# Drop variants whose summary files are missing (e.g. while a refit is
# still running). Lets us produce a 5-panel preview before slopes_si lands.
VARIANTS <- VARIANTS[sapply(VARIANTS, function(v) {
  ok <- file.exists(file.path(SUMMARIES, paste0(v$tag, ".draws.rds")))
  if (!ok) message(sprintf("Skipping variant %s (no draws.rds): %s",
                           v$tag, v$label))
  ok
})]

## ---- 3. Empirical reference (full Wordbank restricted to bundle items) -
build_full_empirical <- function(item_definitions) {
  if (file.exists(CACHE_FULL_EMP)) {
    cat("Loading cached empirical:", CACHE_FULL_EMP, "\n")
    return(readRDS(CACHE_FULL_EMP))
  }

  # Try wordbankr first (cross-sectional + longitudinal)
  emp <- NULL
  cat("Trying wordbankr live pull...\n")
  options(timeout = max(600, getOption("timeout", 60)))
  raw <- tryCatch(
    get_instrument_data(language = "English (American)", form = "WS",
                        administration_info = TRUE, item_info = TRUE),
    error = function(e) {
      message(sprintf("wordbankr pull failed: %s", conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(raw)) {
    emp <- raw |>
      filter(item_kind == "word",
             item_definition %in% item_definitions,
             age >= AGE_RANGE[1], age <= AGE_RANGE[2],
             !is.na(produces)) |>
      group_by(child_id, age) |>
      summarise(vocab = sum(produces, na.rm = TRUE),
                n_items = n(), .groups = "drop")
  }

  # Fallback: cached longitudinal pull
  if (is.null(emp)) {
    long_items_path <- file.path(PATHS$fits_dir, "long_items.rds")
    if (!file.exists(long_items_path)) {
      stop("wordbankr unreachable and ", long_items_path, " missing.")
    }
    cat("Falling back to cached longitudinal pull:", long_items_path, "\n")
    d_long <- readRDS(long_items_path)
    emp <- d_long |>
      filter(language == "English (American)", form == "WS",
             item %in% item_definitions,
             age >= AGE_RANGE[1], age <= AGE_RANGE[2],
             !is.na(produces)) |>
      group_by(child_id, age) |>
      summarise(vocab = sum(produces, na.rm = TRUE),
                n_items = n(), .groups = "drop")
  }

  cat(sprintf("Empirical: %d admins, median %d items per kid\n",
              nrow(emp), median(emp$n_items)))
  saveRDS(emp, CACHE_FULL_EMP)
  emp
}

emp_df <- build_full_empirical(bundle_items)
J_ITEMS <- sd_b$J  # for y-axis ceiling

## ---- 4. Per-variant simulation -----------------------------------
##
## For each variant, sample N_DRAWS_USE posterior draws of scalar params,
## and for each draw simulate N_KIDS_PER_AGE kids at each age in
## AGE_GRID. Per-kid linear predictor includes:
##   xi_i   ~ N(mu_r, sqrt(sigma_r^2 + sigma_alpha^2))
##   zeta_i ~ N(0, sigma_zeta)
##   s_i    ~ Normal_+(0, sigma_s)  (zero if sigma_s not free)
##   kappa_i = 1 + delta + zeta_i
##   eta_ij = xi_i + log_p_j + log_H + kappa_i * log((t - s - s_i)/a0) - psi_j
##
## Items: bundle's J items with empirical log_p_j and posterior median psi_j.
load_variant <- function(tag) {
  draws_path <- file.path(SUMMARIES, paste0(tag, ".draws.rds"))
  psi_path   <- file.path(SUMMARIES, paste0(tag, "_psi.csv"))

  draws <- readRDS(draws_path)
  if ("draws_df" %in% class(draws)) draws <- as.data.frame(draws)

  # M_best (long_no_freq_slopes) doesn't have sigma_s -- fill with 0.
  if (!"sigma_s" %in% names(draws)) draws$sigma_s <- 0

  # M_best psi: read from a proxy if its own _psi.csv doesn't exist
  if (!file.exists(psi_path)) {
    proxy <- file.path(SUMMARIES, "long_no_freq_slopes_si_psi.csv")
    cat(sprintf("[%s] psi not found; using slopes_si as proxy\n", tag))
    psi <- read_csv(proxy, show_col_types = FALSE)
  } else {
    psi <- read_csv(psi_path, show_col_types = FALSE)
  }

  list(draws = draws, psi = psi$psi_median)
}

simulate_variant <- function(tag, label) {
  cat(sprintf("\n=== Simulating %s ===\n", label))
  v <- load_variant(tag)
  draws <- v$draws
  psi   <- v$psi

  # Sub-sample posterior draws
  set.seed(SEED)
  draw_idx <- sort(sample(seq_len(nrow(draws)), N_DRAWS_USE))

  AGE_GRID <- seq(AGE_RANGE[1], AGE_RANGE[2], by = 0.5)
  log_p    <- sd_b$log_p
  log_H    <- sd_b$log_H
  a0       <- sd_b$a0
  mu_r     <- sd_b$mu_r
  sigma_r  <- sd_b$sigma_r

  # si_corr variants carry rho_xi_s and rho_zeta_s in their draws
  # (from log_irt_long_si_corr.stan's generated quantities). When
  # present, sample (xi, zeta, s_lat) jointly from MVN(0, Sigma) per
  # draw with the LKJ-recovered correlation matrix, then apply Tobit
  # clipping s_i = fmax(s_lat, 0). For other variants, fall back to
  # the independent-sampling path.
  is_si_corr <- all(c("rho_xi_s", "rho_zeta_s") %in% names(draws))
  if (is_si_corr) {
    cat(sprintf("[%s] si_corr variant: sampling correlated (xi, zeta, s_lat) per draw\n", tag))
  }

  rows <- vector("list", length(draw_idx) * length(AGE_GRID))
  idx <- 0
  for (k in draw_idx) {
    sa <- as.numeric(draws$sigma_alpha[k])
    sz <- as.numeric(draws$sigma_zeta[k])
    ss <- as.numeric(draws$sigma_s[k])
    delta_d <- as.numeric(draws$delta[k])
    s_d  <- as.numeric(draws$s[k])
    sigma_xi <- sqrt(sa^2 + sigma_r^2)

    if (is_si_corr) {
      # Build 3x3 correlation matrix from this draw's three rho's, then
      # Cholesky-decompose. Symmetric, diag = 1.
      r_xz <- as.numeric(draws$rho_xi_zeta[k])
      r_xs <- as.numeric(draws$rho_xi_s[k])
      r_zs <- as.numeric(draws$rho_zeta_s[k])
      Sigma_corr <- matrix(c(1,    r_xz, r_xs,
                             r_xz, 1,    r_zs,
                             r_xs, r_zs, 1), nrow = 3, byrow = TRUE)
      # Safety: nudge to PSD if near-singular due to draw noise
      L_corr <- tryCatch(chol(Sigma_corr),
                         error = function(e) {
                           Sigma_corr <- Sigma_corr + diag(1e-6, 3)
                           chol(Sigma_corr)
                         })
      L_corr <- t(L_corr)  # lower triangular
      scales <- c(sigma_xi, sz, ss)
    }

    for (a in AGE_GRID) {
      idx <- idx + 1
      set.seed(SEED + k * 1000 + as.integer(a * 10))
      if (is_si_corr) {
        # Sample standardized 3-vectors per kid, scale by L_corr and sigmas.
        Z <- matrix(rnorm(N_KIDS_PER_AGE * 3), nrow = 3)
        effs <- t(L_corr %*% Z) %*% diag(scales)  # N x 3
        xi    <- mu_r + effs[, 1]
        zeta  <- effs[, 2]
        s_lat <- effs[, 3]
        s_i   <- pmax(s_lat, 0)  # Tobit clip
      } else {
        xi   <- rnorm(N_KIDS_PER_AGE, mu_r, sigma_xi)
        zeta <- rnorm(N_KIDS_PER_AGE, 0, sz)
        # half-normal: abs of N(0, sigma_s)
        s_i  <- abs(rnorm(N_KIDS_PER_AGE, 0, ss))
      }
      kappa <- 1 + delta_d + zeta
      log_age <- log(pmax(a - s_d - s_i, 0.01) / a0)
      theta   <- xi + kappa * log_age   # per-kid latent at this age
      # All variants are no_freq (beta_c pinned at 0), so log_p drops
      # from eta. Fitted psi absorbs any frequency-related shift.
      # eta_ij = theta_i + log_H - psi_j
      base_kid <- theta + log_H
      eta <- outer(base_kid, psi, "-")
      vocab <- rowSums(plogis(eta))
      rows[[idx]] <- data.frame(variant = tag, age = a, vocab = vocab)
    }
  }
  bind_rows(rows)
}

sim_list <- lapply(names(VARIANTS), function(nm) {
  simulate_variant(VARIANTS[[nm]]$tag, VARIANTS[[nm]]$label) |>
    mutate(variant_key = nm)
})
sim_df <- bind_rows(sim_list)

## ---- 5. Quantile summaries ---------------------------------------
sim_preds <- sim_df |>
  group_by(variant_key, age) |>
  summarise(
    `0.1`  = quantile(vocab, 0.10, na.rm = TRUE),
    `0.25` = quantile(vocab, 0.25, na.rm = TRUE),
    `0.5`  = quantile(vocab, 0.50, na.rm = TRUE),
    `0.75` = quantile(vocab, 0.75, na.rm = TRUE),
    `0.9`  = quantile(vocab, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(`0.1`, `0.25`, `0.5`, `0.75`, `0.9`),
               names_to = "tau", values_to = "vocab")

# Empirical gcrq
fit_gcrq <- function(df_sub) {
  tryCatch(
    quantregGrowth::gcrq(vocab ~ ps(age, monotone = 1, lambda = 1000),
                         tau = TAUS, data = df_sub),
    error = function(e) NULL
  )
}
pred_gcrq <- function(model, age_grid) {
  if (is.null(model)) return(NULL)
  preds <- predict(model, newdata = data.frame(age = age_grid))
  out <- as.data.frame(preds)
  names(out) <- as.character(TAUS)
  out$age <- age_grid
  pivot_longer(out, cols = -age, names_to = "tau", values_to = "vocab")
}

emp_q    <- fit_gcrq(emp_df)
age_grid <- seq(AGE_RANGE[1], AGE_RANGE[2], by = 0.25)
emp_pred <- pred_gcrq(emp_q, age_grid)

## ---- 6. Plot -----------------------------------------------------
variant_levels <- sapply(VARIANTS, function(v) v$label)
variant_keys   <- names(VARIANTS)
variant_map    <- setNames(variant_levels, variant_keys)

sim_preds$panel <- factor(variant_map[sim_preds$variant_key],
                          levels = variant_levels)

# Empirical lines: one copy per panel
emp_lines_panel <- bind_rows(lapply(variant_levels, function(p) {
  emp_pred |> mutate(panel = p)
}))
emp_lines_panel$panel <- factor(emp_lines_panel$panel, levels = variant_levels)

# Stratified-by-age scatter (~80 per age bin)
set.seed(SEED)
emp_scatter <- emp_df |>
  mutate(age_int = round(age)) |>
  group_by(age_int) |>
  group_modify(~ slice_sample(.x, n = min(nrow(.x), 80))) |>
  ungroup() |>
  select(-age_int)
panel_scatter <- bind_rows(lapply(variant_levels, function(p) {
  emp_scatter |> mutate(panel = p)
}))
panel_scatter$panel <- factor(panel_scatter$panel, levels = variant_levels)
cat(sprintf("Scatter: %d points across %d panels\n",
            nrow(panel_scatter), length(variant_levels)))

p_main <- ggplot(panel_scatter, aes(age, vocab)) +
  geom_jitter(width = 0.3, alpha = 0.20, size = 0.5, colour = "grey30") +
  geom_line(data = emp_lines_panel,
            aes(age, vocab, group = tau),
            colour = "grey45", linetype = "dashed",
            linewidth = 0.5) +
  geom_line(data = sim_preds,
            aes(age, vocab, group = tau, colour = tau),
            linewidth = 1.0) +
  facet_wrap(~ panel, ncol = 4) +
  scale_colour_manual(
    values = c("0.1" = "#1f78b4", "0.25" = "#a6cee3",
               "0.5" = "#33a02c", "0.75" = "#fdbf6f",
               "0.9" = "#e31a1c"),
    name = "Percentile"
  ) +
  scale_x_continuous(breaks = c(16, 20, 24, 28)) +
  coord_cartesian(xlim = c(AGE_RANGE[1] - 0.4, AGE_RANGE[2] + 0.4),
                  ylim = c(0, J_ITEMS)) +
  labs(x = "Age (months)",
       y = sprintf("Productive vocabulary (of %d bundle items)", J_ITEMS),
       title = "Predicted vs. empirical percentile fans across model variants",
       subtitle = "Empirical (Wordbank WS, restricted to bundle items): grey dots + dashed lines.  Predicted (each variant): coloured lines.") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        legend.position = "bottom")

out_png <- file.path(OUT_DIR, "quantile_demo.png")
ggsave(out_png, p_main, width = 12, height = 7.0, dpi = 200)
cat(sprintf("\nWrote: %s\n", out_png))
