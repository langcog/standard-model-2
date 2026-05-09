## 4-panel quantile-fan demo: empirical Wordbank vs.\ predictions
## under four nested model variants.
##
## Variants (each fit to the SAME Wordbank English WS subsample):
##   demo_pure   - pure accumulator: sigma_alpha=0, sigma_zeta=0, delta=0
##   demo_alpha  - + child efficiency variation (sigma_alpha free)
##   demo_kappa  - + age dynamics (delta free, sigma_zeta free; alpha pinned)
##   demo_full   - + both (= M_best, all free)
##
## For each variant, we:
##   1. Fit the model on a subsample of Wordbank English WS data
##   2. Sample posterior trajectories for N=2000 simulated children
##   3. Compute predicted vocabulary at each (child, age) by summing
##      sigmoid(eta_ij) over CDI items
##   4. Run quantile regression (gcrq) on simulated points
##
## Compare the resulting percentile lines to the same regression run
## on the empirical Wordbank data.
##
## Output: outputs/figs/longitudinal/quantile_demo.png
##
## Run:  Rscript model/scripts/quantile_demo.R

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(wordbankr); library(quantregGrowth)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CACHE_BUNDLE <- file.path(PATHS$fits_dir, "demo_subset_data.rds")
AGE_RANGE    <- c(16, 30)
N_KIDS       <- 600
N_ITEMS      <- 200
TAUS         <- c(0.10, 0.25, 0.50, 0.75, 0.90)

# ---- 1. Build a clean cross-sectional bundle ----------------------
build_bundle <- function() {
  if (file.exists(CACHE_BUNDLE)) {
    cat("Loading cached bundle:", CACHE_BUNDLE, "\n")
    return(readRDS(CACHE_BUNDLE))
  }
  cat("Pulling Wordbank English (American) WS data...\n")
  raw <- get_instrument_data(language = "English (American)", form = "WS",
                             administration_info = TRUE, item_info = TRUE)
  df <- raw |>
    filter(item_kind == "word",
           age >= AGE_RANGE[1], age <= AGE_RANGE[2]) |>
    select(child_id, age, item_definition, lexical_category, produces) |>
    filter(!is.na(produces))

  set.seed(2026)
  # Subsample admins stratified by integer age bin
  admin_keys <- df |> distinct(child_id, age) |> mutate(age_int = round(age))
  per_bin <- ceiling(N_KIDS / length(unique(admin_keys$age_int)))
  admin_keys <- admin_keys |>
    group_by(age_int) |>
    group_modify(~ slice_sample(.x, n = min(nrow(.x), per_bin))) |>
    ungroup()
  df <- df |> semi_join(admin_keys, by = c("child_id", "age"))

  # Subsample items stratified by lexical class
  item_keys <- df |> distinct(item_definition, lexical_category)
  per_class <- ceiling(N_ITEMS / length(unique(item_keys$lexical_category)))
  item_keys <- item_keys |>
    group_by(lexical_category) |>
    group_modify(~ slice_sample(.x, n = min(nrow(.x), per_class))) |>
    ungroup()
  df <- df |> semi_join(item_keys, by = "item_definition")

  # Re-index
  admin_idx <- df |> distinct(child_id, age) |> arrange(child_id, age) |>
    mutate(ii = row_number())
  item_idx <- df |> distinct(item_definition, lexical_category) |>
    arrange(item_definition) |>
    mutate(jj = row_number())
  class_levels <- sort(unique(item_idx$lexical_category))
  item_idx$cc <- match(item_idx$lexical_category, class_levels)

  df <- df |>
    left_join(admin_idx, by = c("child_id", "age")) |>
    left_join(item_idx, by = c("item_definition", "lexical_category"))

  # Per-item log-probability proxy: use empirical produces rate among
  # 30-mo-olds, mapped to a CHILDES-like log-frequency scale. For the
  # demo (which uses no_freq variants anyway, beta_c = 0), this is just
  # a placeholder; the model doesn't use log_p when beta_c = 0.
  log_p <- rep(-8, max(item_idx$jj))

  stan_data <- list(
    N = nrow(df), I = max(df$ii), J = max(df$jj),
    C = max(item_idx$cc),
    ii = df$ii, jj = df$jj, cc = item_idx$cc,
    y = df$produces,
    age = admin_idx$age,
    log_p = log_p,
    log_H = log(365),
    a0 = 19,
    mu_r = 7.337734, sigma_r = 0.5338644,
    mu_mu_c = 8, sigma_mu_c = 3,
    s_prior_mean = 0.5, s_prior_sd = 0.05,
    delta_prior_mean = 0, delta_prior_sd = 0.5,
    sigma_lambda_prior_sd = 0.001,
    sigma_zeta_prior_sd = 0.001,
    sigma_alpha_prior_sd = 1,
    beta_c_prior_sd = 1, beta_c_prior_mean = 1,
    time_baseline = 1
  )
  bundle <- list(stan_data = stan_data, item_idx = item_idx,
                 admin_idx = admin_idx, df = df,
                 class_levels = class_levels)
  saveRDS(bundle, CACHE_BUNDLE)
  cat("Wrote:", CACHE_BUNDLE, "\n")
  bundle
}

bundle <- build_bundle()
cat(sprintf("Bundle: I=%d, J=%d, N=%d obs\n",
            bundle$stan_data$I, bundle$stan_data$J, bundle$stan_data$N))

# ---- 2. Fit each variant ------------------------------------------
VARIANTS <- c("demo_pure", "demo_alpha", "demo_kappa", "demo_full")

fit_one <- function(variant_name) {
  tag <- sprintf("wordbank_%s", variant_name)
  fit_path <- file.path(PATHS$fits_dir, sprintf("%s.rds", tag))
  if (file.exists(fit_path)) {
    cat(sprintf("\n[%s] cached, loading\n", variant_name))
    return(readRDS(fit_path))
  }
  cat(sprintf("\n===== Fitting %s =====\n", variant_name))
  overrides <- variant_hyperpriors(variant_name)
  cat("Overrides:\n"); str(overrides)
  stan_data <- modifyList(modifyList(bundle$stan_data, DEFAULT_PRIORS),
                          overrides)
  fit <- fit_variant(stan_data, tag = tag,
                     cfg = modifyList(DEFAULT_FIT_CONFIG,
                                      list(chains = 4, iter = 1000, warmup = 500)))
  fit
}

fits <- lapply(VARIANTS, fit_one)
names(fits) <- VARIANTS

# ---- 3. Sample posterior trajectories per variant -----------------
##
## For each variant, draw N_SIM kids' (xi, zeta) from posterior, then
## compute predicted vocabulary count = sum_j sigmoid(eta_ij) at each
## age in AGE_GRID using the FITTED psi_j for that variant.
suppressPackageStartupMessages(library(posterior))

simulate_for_variant <- function(fit, variant_name, N_SIM = 2000) {
  d <- as_draws_df(fit)
  # Pull a few representative draws (post-warmup, every k-th)
  draw_idx <- seq(1, nrow(d), length.out = 30) |> as.integer()
  AGE_GRID <- seq(AGE_RANGE[1], AGE_RANGE[2], by = 0.5)

  # Per-item psi: posterior median across draws.
  psi_cols <- grep("^psi\\[", names(d), value = TRUE)
  psi_med <- sapply(psi_cols, function(c) median(d[[c]]))
  J <- length(psi_med)
  cat(sprintf("[%s] J=%d items in posterior\n", variant_name, J))

  # Per-draw (sigma_alpha, sigma_zeta, delta, s, mu_r) posteriors
  rows <- list()
  for (k in seq_along(draw_idx)) {
    d_idx <- draw_idx[k]
    sa <- as.numeric(d$sigma_alpha[d_idx])
    sz <- ifelse("sigma_zeta" %in% names(d), as.numeric(d$sigma_zeta[d_idx]), 0)
    delta_d <- as.numeric(d$delta[d_idx])
    s_d <- as.numeric(d$s[d_idx])
    sigma_xi <- sqrt(sa^2 + bundle$stan_data$sigma_r^2)
    # Sample kids
    set.seed(2026 + k)
    n_per_age <- ceiling(N_SIM / length(AGE_GRID))
    for (a in AGE_GRID) {
      xi <- rnorm(n_per_age, mean = bundle$stan_data$mu_r, sd = sigma_xi)
      zeta <- rnorm(n_per_age, mean = 0, sd = sz)
      log_age <- log(max(a - s_d, 0.01) / bundle$stan_data$a0)
      kappa <- 1 + delta_d + zeta  # per-kid scaling exponent
      theta <- xi + kappa * log_age
      # Predicted vocab count: sum_j sigmoid(theta + log_p_j + log_H - psi_j)
      # NB: log_irt.stan adds log_p[j] to base with unit coefficient even
      # under no_freq variants (the Stan file doesn't have a beta_c
      # multiplier on log_p — that's only in log_irt_long.stan). So when
      # log_p[j] = -8 (placeholder constant in our bundle), every base
      # carries a -8 offset that the fitted psi_j absorbs. The simulation
      # MUST include log_p too or it over-shifts predictions toward ceiling.
      log_p_vec <- bundle$stan_data$log_p
      base_kid_age <- theta + bundle$stan_data$log_H
      eta_partial <- outer(base_kid_age, psi_med, "-")
      eta_partial <- sweep(eta_partial, 2, log_p_vec, "+")  # add log_p_j per item
      vocab <- rowSums(plogis(eta_partial))
      rows[[length(rows) + 1]] <- data.frame(
        variant = variant_name, draw = d_idx, age = a, vocab = vocab
      )
    }
  }
  bind_rows(rows)
}

cat("\n===== Simulating trajectories =====\n")
sims <- lapply(names(fits), function(nm) simulate_for_variant(fits[[nm]], nm))
sim_df <- bind_rows(sims)
saveRDS(sim_df, file.path(PATHS$fits_dir, "quantile_demo_sims.rds"))
cat("Saved sims to fits/quantile_demo_sims.rds\n")

# ---- 4. gcrq quantile regression on each variant + empirical -----
##
## Use the same bs-spline setup as the wordbank book.
fit_gcrq <- function(df_sub) {
  tryCatch({
    quantregGrowth::gcrq(vocab ~ ps(age, monotone = 1, lambda = 1000),
                         tau = TAUS, data = df_sub)
  }, error = function(e) {
    message("gcrq failed: ", conditionMessage(e))
    NULL
  })
}

pred_gcrq <- function(model, age_grid) {
  if (is.null(model)) return(NULL)
  preds <- predict(model, newdata = data.frame(age = age_grid))
  out <- as.data.frame(preds)
  names(out) <- as.character(TAUS)
  out$age <- age_grid
  pivot_longer(out, cols = -age, names_to = "tau", values_to = "vocab")
}

# Empirical: aggregate child total vocab across the bundle
emp_df <- bundle$df |>
  group_by(child_id, age) |>
  summarise(vocab = sum(produces, na.rm = TRUE), .groups = "drop")

# Run quantile regression
emp_q   <- fit_gcrq(emp_df)
sim_qs  <- lapply(VARIANTS, function(v) {
  fit_gcrq(sim_df |> filter(variant == v))
})
names(sim_qs) <- VARIANTS

# Predicted percentile curves
age_grid <- seq(AGE_RANGE[1], AGE_RANGE[2], by = 0.25)
emp_pred <- pred_gcrq(emp_q, age_grid) |> mutate(variant = "Empirical")
sim_preds <- bind_rows(lapply(VARIANTS, function(v) {
  p <- pred_gcrq(sim_qs[[v]], age_grid)
  if (!is.null(p)) p$variant <- v
  p
}))

# ---- 5. Plot ------------------------------------------------------
variant_labels <- c(
  demo_pure  = "Pure accumulator",
  demo_alpha = "+ efficiency variation",
  demo_kappa = "+ scaling variation",
  demo_full  = "+ both (M_best)"
)
sim_preds$variant <- factor(sim_preds$variant,
                            levels = VARIANTS,
                            labels = variant_labels[VARIANTS])
emp_pred$variant <- factor("Empirical")

# Build 4-panel: each panel shows empirical points + empirical quantile
# lines (faint) + variant's predicted quantile lines (bold).
panel_data <- emp_df |>
  mutate(panel_n = 1) |>
  bind_rows(emp_df |> mutate(panel_n = 2)) |>
  bind_rows(emp_df |> mutate(panel_n = 3)) |>
  bind_rows(emp_df |> mutate(panel_n = 4))

variant_panels <- data.frame(
  panel_n = 1:4,
  variant = factor(variant_labels[VARIANTS], levels = variant_labels[VARIANTS])
)
panel_data <- panel_data |> left_join(variant_panels, by = "panel_n")
sim_preds_panel <- sim_preds |>
  left_join(variant_panels |> rename(variant_match = variant),
            by = c("variant" = "variant_match")) |>
  rename(variant_x = variant)
emp_lines_panel <- bind_rows(lapply(1:4, function(p) {
  emp_pred |> mutate(panel_n = p)
})) |> left_join(variant_panels, by = "panel_n") |>
  rename(variant_x = variant.x, variant = variant.y)

p_main <- ggplot(panel_data, aes(age, vocab)) +
  # Empirical points (light)
  geom_jitter(width = 0.3, alpha = 0.18, size = 0.6, colour = "grey30") +
  # Empirical quantile lines (faint dashed)
  geom_line(data = emp_lines_panel,
            aes(age, vocab, group = tau),
            colour = "grey45", linetype = "dashed",
            linewidth = 0.5) +
  # Predicted quantile lines (bold colored)
  geom_line(data = sim_preds,
            aes(age, vocab, group = tau, colour = tau),
            linewidth = 1.0) +
  facet_wrap(~ variant, ncol = 2) +
  scale_colour_manual(
    values = c("0.1" = "#1f78b4", "0.25" = "#a6cee3",
               "0.5" = "#33a02c", "0.75" = "#fdbf6f",
               "0.9" = "#e31a1c"),
    name = "Predicted percentile"
  ) +
  coord_cartesian(ylim = c(0, max(emp_df$vocab) * 1.05)) +
  labs(x = "Age (months)", y = "Productive vocabulary (count)",
       title = "How well does each model account for the empirical fan?",
       subtitle = "Black points: Wordbank English WS empirical data. Dashed grey: empirical 10/25/50/75/90 percentile curves (gcrq). Coloured: same percentiles from posterior-simulated kids under each variant.") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        legend.position = "bottom")

out_png <- file.path(OUT_DIR, "quantile_demo.png")
ggsave(out_png, p_main, width = 12, height = 8, dpi = 150)
cat(sprintf("\nWrote: %s\n", out_png))
