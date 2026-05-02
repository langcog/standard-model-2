## Effective tokens / real tokens diagnostic.
##
## Make the (1+delta) exponent visceral: with delta ~9, the model says
## effective tokens accumulate as (t/a_0)^10. Real input grows linearly
## in t. So the model is saying "tokens get more valuable" steeply with
## age. Plot:
##
##   1. For example words at low/mid/high frequency, plot CUMULATIVE
##      REAL TOKENS heard by the population-mean child as a function of
##      age. Mark the age at which the model predicts 50% production.
##
##   2. "Real tokens at 50% acquisition" vs. age of 50% acquisition,
##      colored by log p_j. Tests Mike's intuition: do words acquired
##      early need more real tokens than words acquired later?
##
## Output: model/figs/longitudinal/effective_tokens.png

source("model/R/config.R")
source("model/R/helpers.R")
source("model/R/datasets.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(posterior); library(rstan)
})

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

# ---- Load lean reference English fit and bundle ---- #
bundle <- load_dataset_bundle("english")
fit    <- readRDS(file.path(PATHS$fits_dir, "long_slopes.rds"))
sd_    <- bundle$stan_data
d      <- as_draws_df(fit)
log_p  <- log(bundle$word_info$prob)
words  <- bundle$word_info$item

med <- function(p) median(d[[p]])
mu_r   <- sd_$mu_r
log_H  <- sd_$log_H
a0     <- sd_$a0
s_med  <- med("s")
delta  <- med("delta")
psi    <- sapply(seq_len(sd_$J),
                 function(j) median(d[[sprintf("psi[%d]", j)]]))
beta_j <- psi - log_p - log_H

cat(sprintf("Population: mu_r=%.2f, delta=%.2f, s=%.2f, a_0=%d\n",
            mu_r, delta, s_med, a0))

# Inverse Sperry-style mean input rate (tokens/hour)
mean_rate_per_hr <- exp(mu_r)
# Real tokens of word j heard by population-mean kid at age t (t > s)
real_tokens <- function(t, log_p_j) {
  pmax(t - s_med, 0.01) * exp(mu_r) * exp(log_p_j) * exp(log_H) / exp(log_H) * exp(log_H)
  # = (t-s) * mean_rate * p_j * H , in tokens
  # Carefully: mean_rate is tokens/hour, H is hours/month
  # so tokens/month = mean_rate * H * p_j; tokens by age t (months) = mean_rate*H*p_j*(t-s)
}
# cleaner:
real_tokens <- function(t, log_p_j) {
  exp(mu_r) * exp(log_H) * exp(log_p_j) * pmax(t - s_med, 0.01)
}

# Effective tokens = exp(theta_pop) at population mean child
# theta_pop(t) = mu_r + log_H + (1+delta) * log((t-s)/a_0)  + log p_j on item side
# But for "effective tokens of word j" specifically:
#   effective_j(t) = mean_rate * H * p_j * a_0 * ((t-s)/a_0)^(1+delta)
# (this is the model's effective-tokens-of-word-j proxy)
eff_tokens <- function(t, log_p_j) {
  exp(mu_r) * exp(log_H) * exp(log_p_j) * a0 *
    (pmax(t - s_med, 0.01) / a0)^(1 + delta)
}

# 50% acquisition age for population-mean kid:
#   theta_pop(t)            = beta_j
#   mu_r + (1+delta)*log((t-s)/a_0)  =  psi_j - log p_j - log H
#   solve: t-s = a_0 * exp((psi_j - log_p_j - log_H - mu_r) / (1+delta))
age_50 <- function(j) {
  s_med + a0 * exp((psi[j] - log_p[j] - log_H - mu_r) / (1 + delta))
}

# ---- Pick example words at low/mid/high frequency ---- #
# Use frequency quintiles of log_p, pick words with class = noun if possible
# for a familiar exemplar
quintile_idx <- function(q) which(rank(log_p) >= q*sd_$J & rank(log_p) < (q+0.05)*sd_$J)

# Pick 6 words spanning the age50 distribution.  Compute age50 for all
# words first, then pick representatives at quantiles of age50.
all_age50 <- sapply(seq_len(sd_$J), age_50)
all_tokens50 <- mapply(function(t, lp) exp(mu_r) * exp(log_H) * exp(lp) *
                                       max(t - s_med, 0.01),
                       all_age50, log_p)
in_range <- which(all_age50 >= 8 & all_age50 <= 35)
sorted   <- in_range[order(all_age50[in_range])]
n_in <- length(sorted)
picks <- sorted[round(seq(0.05, 0.95, length.out = 6) * n_in)]
example_j <- picks
ex <- tibble(j = example_j,
             word = words[example_j],
             log_p = log_p[example_j],
             psi = psi[example_j],
             beta = beta_j[example_j],
             age50 = all_age50[example_j]) %>%
  arrange(age50)
cat("Example words (spanning age50):\n"); print(ex)

# ---- Panel 1: cumulative tokens vs age, with 50% mark ---- #
ages <- seq(6, 32, by = 0.25)
trajectories <- bind_rows(lapply(seq_len(nrow(ex)), function(k) {
  tibble(j     = ex$j[k],
         word  = ex$word[k],
         log_p = ex$log_p[k],
         age50 = ex$age50[k],
         age   = ages,
         real_tokens = real_tokens(ages, ex$log_p[k]))
})) %>%
  mutate(word_label = sprintf("%s  (log p = %.1f, age50 = %.1f)",
                              word, log_p, age50))

# Per-word: tokens at age50
tokens_at_50 <- ex %>%
  mutate(real_tokens_at_50 = real_tokens(age50, log_p),
         word_label = sprintf("%s  (log p = %.1f, age50 = %.1f)",
                              word, log_p, age50))

p1 <- ggplot(trajectories, aes(x = age, y = real_tokens, color = word_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(data = tokens_at_50,
             aes(x = age50, y = real_tokens_at_50, color = word_label),
             size = 2.5, shape = 21, stroke = 1.2, fill = "white") +
  scale_y_log10(labels = function(x) format(x, scientific = TRUE, digits = 1)) +
  scale_color_viridis_d(option = "plasma", end = 0.85, name = NULL) +
  coord_cartesian(xlim = c(6, 32)) +
  labs(x = "age (months)",
       y = "real tokens of word j heard by population-mean child\n(log scale; r * p_j * H * (t-s))",
       title = "(1) Real tokens of word j heard vs. age",
       subtitle = paste0("Open circle = age of 50% production (where ",
                          "theta_pop crosses beta_j).  Vast spread in ",
                          "real-tokens-at-acquisition across words.")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

# ---- Panel 2: tokens-at-50 vs age-of-50 across ALL words ---- #
# (already computed above)
all_tab <- tibble(j = seq_len(sd_$J), word = words, log_p = log_p,
                  age50 = all_age50, tokens_at_50 = all_tokens50,
                  class = bundle$class_levels[bundle$word_info$cc]) %>%
  filter(age50 >= 6, age50 <= 36)   # plot reasonable age range

cat(sprintf("\nWords in plot range: %d / %d\n",
            nrow(all_tab), sd_$J))

p2 <- ggplot(all_tab, aes(x = age50, y = tokens_at_50, color = log_p)) +
  geom_point(alpha = 0.7, size = 1.4) +
  scale_y_log10(labels = function(x) {
    ifelse(x < 1, sprintf("%.1f", x), format(x, big.mark = ",", scientific = FALSE))
  }) +
  scale_color_viridis_c(option = "plasma", end = 0.9,
                        name = expression(log~p[j])) +
  geom_hline(yintercept = c(1, 10, 100, 1000, 10000, 100000),
             color = "gray85", linewidth = 0.2) +
  labs(x = "age of 50% production (months)",
       y = "real tokens of word j heard at age50\n(log scale)",
       title = "(2) Real tokens to 50% production vs age, all words",
       subtitle = paste0("Each dot = one of J=", sd_$J,
                          " items.  Frequent words (yellow) acquired ",
                          "early with many tokens; rare words ",
                          "(purple) acquired late with few tokens. ",
                          "Red line = median by age bin.")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Annotations: at three example ages, what's the typical tokens-at-acquisition
trend_summary <- all_tab %>%
  mutate(age_bin = round(age50)) %>%
  group_by(age_bin) %>%
  summarise(n = n(), median_tokens = median(tokens_at_50),
            .groups = "drop") %>%
  filter(n >= 5)

p2_with_trend <- p2 +
  geom_line(data = trend_summary,
            aes(x = age_bin, y = median_tokens),
            color = "firebrick", linewidth = 0.8, inherit.aes = FALSE)

# ---- Combine ---- #
out <- p1 / p2_with_trend + plot_layout(heights = c(1, 1.1))
ggsave(file.path(OUT_FIGS, "effective_tokens.png"),
       out, width = 11, height = 9, dpi = 200)
cat("Wrote effective_tokens.png\n")

# Print quantitative summary
cat("\n--- Summary ---\n")
cat(sprintf("Mean rate (tokens/hr): %.0f\n", mean_rate_per_hr))
cat(sprintf("Tokens of mean-freq word heard by 30 mo:\n"))
mid_log_p <- median(log_p)
cat(sprintf("  At log p = %.2f: %.0f tokens by age 30\n",
            mid_log_p, real_tokens(30, mid_log_p)))

cat(sprintf("\nReal tokens at 50%% acquisition, percentiles across all %d words:\n",
            sd_$J))
print(quantile(all_tab$tokens_at_50,
               probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
               na.rm = TRUE), digits = 4)

cat(sprintf("\nReal tokens at 50%% acquisition vs age:\n"))
print(trend_summary %>% filter(age_bin %in% c(12, 16, 20, 24, 28, 32)),
      n = 20)
