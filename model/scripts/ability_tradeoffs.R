## Two diagnostic figures for the ability side of the model:
##
##   A. (s, delta) joint posterior on the free_s_slopes variant.
##      Tests whether s and delta trade off along a 1-d ridge or are
##      roughly orthogonal.
##
##   B. Per-child (xi, zeta) joint posterior, English vs. Norwegian.
##      Diagnoses the rho_xi_zeta sign flip (+0.35 vs -0.31). Color
##      kids by admin count to test whether sparse-coverage kids drive
##      the correlation; provide robust-subset r as a sanity check.
##
## Output:
##   outputs/figs/longitudinal/ability_s_delta_tradeoff.png
##   outputs/figs/longitudinal/ability_xi_zeta_languages.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(posterior); library(rstan)
})

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# A. (s, delta) joint posterior
# =====================================================================

cat("\n=== A. (s, delta) joint posterior ===\n")

extract_s_delta <- function(variant, dataset_label) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) { cat("  missing:", path, "\n"); return(NULL) }
  fit <- readRDS(path)
  d <- as_draws_df(fit)
  tibble(variant = variant, dataset = dataset_label,
         s = d$s, delta = d$delta)
}

ab_a <- bind_rows(
  extract_s_delta("long_slopes",                "english"),
  extract_s_delta("long_free_s_slopes",         "english"),
  extract_s_delta("long_slopes_norwegian",      "norwegian"),
  extract_s_delta("long_free_s_slopes_norwegian","norwegian")
) %>%
  mutate(variant_label = factor(
    case_when(grepl("free_s", variant) ~ "free s",
              TRUE ~ "lean ref (s pinned ~0.5)"),
    levels = c("lean ref (s pinned ~0.5)", "free s")),
    dataset = factor(dataset, levels = c("english", "norwegian")))

# Per-panel correlation
cors_a <- ab_a %>% group_by(variant_label, dataset) %>%
  summarise(r = cor(s, delta), n = n(), .groups = "drop") %>%
  mutate(label = sprintf("r=%.2f", r))
cat("Per-panel r(s, delta):\n"); print(as.data.frame(cors_a))

p_a <- ggplot(ab_a, aes(x = s, y = delta)) +
  geom_point(alpha = 0.10, size = 0.5, color = "steelblue") +
  geom_density_2d(color = "navy", linewidth = 0.4, bins = 8) +
  geom_text(data = cors_a, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
            inherit.aes = FALSE, size = 3.5) +
  facet_grid(dataset ~ variant_label, scales = "free") +
  labs(x = "s  (start time, months)",
       y = expression(delta ~ "  (population age rate-change exponent)"),
       title = expression("(A) Joint posterior of " ~ s ~ " and " ~ delta),
       subtitle = paste0("Lean ref pins s near 0.5; free_s lets it ",
                          "wander. The correlation in free_s ",
                          "reveals the trade-off ridge.")) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(OUT_FIGS, "ability_s_delta_tradeoff.png"),
       p_a, width = 9, height = 6, dpi = 200)
cat("Wrote ability_s_delta_tradeoff.png\n")

# =====================================================================
# B. (xi, zeta) per kid, English vs Norwegian
# =====================================================================

cat("\n=== B. (xi, zeta) by language ===\n")

extract_xi_zeta_per_kid <- function(variant, dataset_key, language_label) {
  path <- file.path(PATHS$fits_dir, sprintf("%s.rds", variant))
  if (!file.exists(path)) { cat("  missing:", path, "\n"); return(NULL) }
  fit <- readRDS(path)
  d <- as_draws_df(fit)
  bundle <- load_dataset_bundle(dataset_key)
  I <- bundle$stan_data$I
  J <- bundle$stan_data$J
  a0 <- bundle$stan_data$a0
  s_med <- median(d$s)
  xi   <- sapply(seq_len(I), function(i) median(d[[sprintf("xi[%d]", i)]]))
  zeta <- sapply(seq_len(I), function(i) median(d[[sprintf("zeta[%d]", i)]]))

  # Per-kid coverage / ceiling info
  df <- as_tibble(bundle$df)
  admin_info <- as_tibble(bundle$admin_info)
  # Per-admin total produced (max J=200)
  admin_totals <- df %>% group_by(aa) %>%
    summarise(produced = sum(produces), .groups = "drop") %>%
    left_join(admin_info %>% select(aa, ii, age), by = "aa")
  per_kid <- admin_totals %>% group_by(ii) %>%
    summarise(n_admins = n(),
              age_min  = min(age),
              age_max  = max(age),
              age_span = max(age) - min(age),
              med_age  = median(age),
              max_vocab_frac = max(produced) / J,
              .groups = "drop")

  out <- tibble(language = language_label, ii = seq_len(I), xi = xi, zeta = zeta) %>%
    left_join(per_kid, by = "ii")
  # Re-centered xi: evaluate the per-kid ability at the kid s own median
  # admin age instead of at a_0. Removes the ridge tilt that contaminates
  # marginal r(xi, zeta) when admin ages skew far from a_0.
  out <- out %>%
    mutate(log_mid = log(pmax(med_age - s_med, 0.01) / a0),
           xi_centered = xi + zeta * log_mid)
  out
}

xz <- bind_rows(
  extract_xi_zeta_per_kid("long_slopes",           "english",   "English"),
  extract_xi_zeta_per_kid("long_slopes_norwegian", "norwegian", "Norwegian")
) %>%
  mutate(language = factor(language, levels = c("English", "Norwegian")))

# Per-language correlations -- raw + re-centered
cors_b <- xz %>% group_by(language) %>%
  summarise(
    r_raw       = cor(xi, zeta),
    r_dense     = cor(xi[n_admins >= 4], zeta[n_admins >= 4]),
    r_noceil    = cor(xi[max_vocab_frac < 0.95], zeta[max_vocab_frac < 0.95]),
    r_robust    = cor(xi[n_admins >= 4 & max_vocab_frac < 0.95],
                      zeta[n_admins >= 4 & max_vocab_frac < 0.95]),
    r_centered  = cor(xi_centered, zeta),
    n_all     = n(),
    n_dense   = sum(n_admins >= 4),
    n_noceil  = sum(max_vocab_frac < 0.95),
    n_robust  = sum(n_admins >= 4 & max_vocab_frac < 0.95),
    .groups = "drop"
  )
cat("Per-language correlations:\n"); print(as.data.frame(cors_b))

# ---- Panel B1: scatter colored by admin count ---- #
xz_lab <- xz %>% mutate(label = sprintf("r = %.2f  (n=%d)",
                                        cor(xi, zeta), n()))

p_b1 <- ggplot(xz, aes(x = xi, y = zeta, color = factor(n_admins))) +
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
  geom_vline(xintercept = NULL) +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_density_2d(aes(group = 1), color = "gray30",
                  linewidth = 0.4, bins = 6, inherit.aes = TRUE,
                  show.legend = FALSE) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              color = "firebrick", linewidth = 0.6, linetype = "dashed") +
  geom_text(data = xz %>% group_by(language) %>%
              summarise(r = cor(xi, zeta), n = n(), .groups = "drop") %>%
              mutate(label = sprintf("r = %+.2f  (n=%d kids)", r, n)),
            aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
            inherit.aes = FALSE, color = "firebrick", size = 3.6) +
  facet_wrap(~language, nrow = 1) +
  scale_color_brewer(palette = "YlOrRd", name = "n admins") +
  labs(x = expression(xi[i] ~ "  (per-child intercept, posterior median)"),
       y = expression(zeta[i] ~ "  (per-child slope deviation)"),
       title = "(B) Per-child (xi, zeta): English vs. Norwegian",
       subtitle = paste0("Color = how many longitudinal admins each ",
                          "kid contributed (3-9). Density contours and ",
                          "linear fit overlaid.")) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "right")

# ---- Panel B2: robustness table + re-centered correlation ---- #
robust_long <- cors_b %>%
  mutate(n_raw = n_all, n_centered = n_all) %>%
  select(language, r_raw, r_dense, r_noceil, r_robust, r_centered,
         n_raw, n_dense, n_noceil, n_robust, n_centered) %>%
  pivot_longer(cols = -language,
               names_to = c(".value", "subset"),
               names_pattern = "(r|n)_(.*)") %>%
  mutate(subset = factor(subset,
                         levels = c("raw", "dense", "noceil", "robust", "centered"),
                         labels = c("raw (xi at a_0=20)",
                                    "n_admins >= 4",
                                    "no ceiling (max_vocab < 95%)",
                                    "both filters",
                                    "xi re-centered at per-kid median age")))

p_b2 <- ggplot(robust_long,
               aes(y = subset, x = r, fill = language, label = sprintf("r=%+.2f (n=%d)", r, n))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           alpha = 0.8) +
  geom_text(position = position_dodge(width = 0.8),
            color = "black", size = 3, hjust = -0.05) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.3) +
  scale_fill_manual(values = c("English" = "#1f77b4",
                                "Norwegian" = "#d62728")) +
  coord_cartesian(xlim = c(-0.7, 0.7)) +
  labs(x = expression("Pearson correlation " ~ r(xi[i], zeta[i])),
       y = NULL,
       title = "Sign flip diagnoses as a parameterization artifact",
       subtitle = paste0("Subset filters preserve the flip; ",
                          "re-centering xi to each kid's mid-observation ",
                          "age makes it vanish.")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

# ---- Panel B3: re-centered scatter ---- #
p_b3 <- ggplot(xz, aes(x = xi_centered, y = zeta)) +
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
  geom_point(alpha = 0.55, size = 1.4, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE,
              color = "firebrick", linewidth = 0.6, linetype = "dashed") +
  geom_text(data = xz %>% group_by(language) %>%
              summarise(r = cor(xi_centered, zeta), n = n(), .groups = "drop") %>%
              mutate(label = sprintf("r = %+.2f", r)),
            aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
            inherit.aes = FALSE, color = "firebrick", size = 3.6) +
  facet_wrap(~language, nrow = 1) +
  labs(x = expression("re-centered " ~ xi[i] ~
                      " (evaluated at per-kid median age)"),
       y = expression(zeta[i]),
       title = "Re-centered scatter: ridge tilt removed",
       subtitle = "Both languages now show similar weak positive correlation.") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# Combine
p_b <- p_b1 / p_b3 / p_b2 + plot_layout(heights = c(1.4, 1.0, 1.0))

ggsave(file.path(OUT_FIGS, "ability_xi_zeta_languages.png"),
       p_b, width = 11, height = 9, dpi = 200)
cat("Wrote ability_xi_zeta_languages.png\n")
