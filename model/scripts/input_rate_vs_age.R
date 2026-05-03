## Descriptive: does observed input rate r_i vary with age?
##
## The model assumes log_r_i is a per-child constant (no age-varying
## input rate growth). The (1+delta+zeta_i)*log_age term then absorbs
## any age-varying child-side dynamics into "efficiency growth". This
## is justifiable only if log r_i really is roughly age-invariant.
##
## Test: plot per-recording log_r_obs vs age, both for BabyView (head-cam,
## adult tokens / video duration) and SEEDLingS (LENA AWC), with per-kid
## within-subject linear fits. If slopes are clustered around zero, the
## age-invariance assumption holds.
##
## Output: outputs/figs/io/input_rate_vs_age.png

source("model/R/config.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(tidyr)
})

OUT_FIGS <- file.path(PATHS$figs_dir, "io")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

# ---- BabyView ---- #
bv <- load_dataset_bundle("babyview")
v_bv <- bv$videos %>% as_tibble() %>%
  select(child_ii, age_mo, log_r_obs) %>%
  mutate(dataset = "BabyView (head-cam, FEM)")

cat(sprintf("BabyView: %d videos, %d kids\n",
            nrow(v_bv), length(unique(v_bv$child_ii))))

# ---- SEEDLingS ---- #
sd <- load_dataset_bundle("seedlings")
v_sd <- sd$recordings %>% as_tibble() %>%
  select(child_ii, age_mo = month, log_r_obs) %>%
  mutate(dataset = "SEEDLingS (LENA AWC)")

cat(sprintf("SEEDLingS: %d recordings, %d kids\n",
            nrow(v_sd), length(unique(v_sd$child_ii))))

all_rec <- bind_rows(v_bv, v_sd)

# ---- Per-kid within-subject slopes ---- #
fit_kid_slope <- function(df) {
  if (nrow(df) < 3) return(tibble(slope = NA, intercept = NA, n = nrow(df)))
  m <- lm(log_r_obs ~ age_mo, data = df)
  tibble(slope = coef(m)[2], intercept = coef(m)[1], n = nrow(df))
}
kid_fits <- all_rec %>% group_by(dataset, child_ii) %>%
  do(fit_kid_slope(.)) %>% ungroup()

cat("\nWithin-kid slope of log_r_obs vs age (logits/month):\n")
print(kid_fits %>% group_by(dataset) %>%
      summarise(n_kids = n(),
                mean_slope = mean(slope, na.rm = TRUE),
                median_slope = median(slope, na.rm = TRUE),
                sd_slope = sd(slope, na.rm = TRUE),
                pct_slope_pos = mean(slope > 0, na.rm = TRUE) * 100,
                .groups = "drop"))

# ---- Pop-level slope (cross-sectional, ignoring kid identity) ---- #
pop_fits <- all_rec %>% group_by(dataset) %>%
  do(broom::tidy(lm(log_r_obs ~ age_mo, data = .))) %>%
  filter(term == "age_mo") %>% ungroup()
cat("\nPopulation-level (between+within) slope of log_r_obs vs age:\n")
print(pop_fits %>% select(dataset, slope_estimate = estimate, se = std.error,
                           p = p.value))

# ---- Plot ---- #
# For visualizing: jitter age slightly so all videos visible
all_rec_pl <- all_rec %>% mutate(age_jit = age_mo + runif(n(), -0.1, 0.1))

p <- ggplot(all_rec_pl,
            aes(x = age_jit, y = log_r_obs, color = factor(child_ii))) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_smooth(aes(group = child_ii), method = "lm", se = FALSE,
              linewidth = 0.4, alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "black", linewidth = 1, fill = "gray60",
              alpha = 0.25) +
  facet_wrap(~dataset, scales = "free", nrow = 1) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  labs(x = "age (months)",
       y = expression(log~r[obs] ~ "  (per recording)"),
       title = "Per-recording observed input rate vs. age",
       subtitle = paste0("Each thin colored line = one child's within-",
                          "subject regression. Black line = pooled. ",
                          "If log r is roughly age-invariant, slopes ",
                          "should cluster near zero.")) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# Histogram of within-kid slopes
p_hist <- ggplot(kid_fits, aes(x = slope, fill = dataset)) +
  geom_histogram(bins = 12, alpha = 0.8, color = "white") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "gray30") +
  facet_wrap(~dataset, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("BabyView (head-cam, FEM)" = "#1f77b4",
                                "SEEDLingS (LENA AWC)" = "#d62728"),
                     guide = "none") +
  labs(x = "within-kid slope (logits/month)",
       y = "count of kids",
       title = "Distribution of per-kid slopes",
       subtitle = "Mass concentrated near zero supports age-invariance assumption.") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

out <- p / p_hist + plot_layout(heights = c(1.4, 1))
ggsave(file.path(OUT_FIGS, "input_rate_vs_age.png"),
       out, width = 11, height = 8, dpi = 200)
cat("\nWrote input_rate_vs_age.png\n")
