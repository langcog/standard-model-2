## RQ2 / out-of-vocabulary generalization test.
##
## Runs the regression of per-word log-threshold (psi_j) on log(p_j) +
## lexical class, at the level of M1 (= long_m1_time_only): the spine
## stage where time has been added but frequency has NOT. At that
## stage psi_j captures raw word difficulty in token-time units, with
## no portion already absorbed by a structural log p_j coefficient.
##
## The R^2 of this regression tells us how much of word-level difficulty
## variance frequency + class explains. Critically, this is a different
## quantity from the LOO-ELPD comparison in the spine: ELPD measures
## how much frequency adds for predicting *held-out observations* of
## already-seen words (where psi_j is itself estimated), while R^2
## here measures how much frequency would help predict *new words'*
## difficulty (where psi_j is unknown).
##
## Usage:
##   Rscript model/scripts/psi_freq_regression.R [fit_path] [bundle_path]
## Defaults: fits/long_m1_time_only.rds, fits/long_subset_data.rds
##
## Outputs:
##   outputs/figs/longitudinal/psi_freq_regression.csv
##   outputs/figs/longitudinal/psi_freq_regression.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
fit_path    <- if (length(args) >= 1) args[1]
                else file.path(PATHS$fits_dir, "long_m1_time_only.rds")
bundle_path <- if (length(args) >= 2) args[2]
                else file.path(PATHS$fits_dir, "long_subset_data.rds")

OUT_FIGS <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(fit_path))
  stop(sprintf("Fit not found: %s\n  (Sync from Sherlock when the m1_time_only fit lands.)", fit_path))
if (!file.exists(bundle_path))
  stop(sprintf("Bundle not found: %s", bundle_path))

cat(sprintf("Reading fit: %s\n", fit_path))
fit    <- readRDS(fit_path)
bundle <- readRDS(bundle_path)

d <- as_draws_df(fit)

psi_cols <- grep("^psi\\[", names(d), value = TRUE)
if (length(psi_cols) == 0) stop("No psi[] columns in fit draws.")
cat(sprintf("Found %d per-word psi parameters.\n", length(psi_cols)))

# Per-word posterior summary
psi_summary <- tibble(
  jj      = as.integer(sub("psi\\[(\\d+)\\]", "\\1", psi_cols)),
  psi_mean   = sapply(psi_cols, function(p) mean(d[[p]])),
  psi_median = sapply(psi_cols, function(p) median(d[[p]])),
  psi_sd     = sapply(psi_cols, function(p) sd(d[[p]])),
  psi_q025   = sapply(psi_cols, function(p) quantile(d[[p]], 0.025, names = FALSE)),
  psi_q975   = sapply(psi_cols, function(p) quantile(d[[p]], 0.975, names = FALSE))
)

# Join with word_info from the bundle to get item, log_p, lexical class
classes <- bundle$class_levels
word_info <- bundle$word_info %>%
  mutate(log_p = log(prob),
         lexical_class = classes[cc])

dat <- psi_summary %>% left_join(word_info, by = "jj")
cat(sprintf("Joined: %d words.\n", nrow(dat)))

# ---- Regression ---- #
m_full     <- lm(psi_median ~ log_p + lexical_class, data = dat)
m_log_p    <- lm(psi_median ~ log_p, data = dat)
m_class    <- lm(psi_median ~ lexical_class, data = dat)
m_null     <- lm(psi_median ~ 1, data = dat)

cat("\n=== Regression: psi_j ~ log p_j + lexical_class ===\n")
print(summary(m_full))

cat(sprintf("\nR^2 (log_p alone)            = %.3f\n",
            summary(m_log_p)$r.squared))
cat(sprintf("R^2 (lexical_class alone)    = %.3f\n",
            summary(m_class)$r.squared))
cat(sprintf("R^2 (log_p + lexical_class)  = %.3f\n",
            summary(m_full)$r.squared))
cat(sprintf("Adj R^2 (log_p + class)      = %.3f\n",
            summary(m_full)$adj.r.squared))

# ---- Per-class breakdown ---- #
cat("\n=== Per-class regression of psi_j on log p_j ===\n")
per_class <- dat %>%
  group_by(lexical_class) %>%
  summarise(n = n(),
            slope_log_p = coef(lm(psi_median ~ log_p))[2],
            r2 = summary(lm(psi_median ~ log_p))$r.squared,
            cor_psi_logp = cor(psi_median, log_p),
            .groups = "drop") %>%
  mutate(across(c(slope_log_p, r2, cor_psi_logp), ~ round(.x, 3)))
print(per_class)

# ---- Save CSV ---- #
out_csv <- file.path(OUT_FIGS, "psi_freq_regression.csv")
results <- list(
  per_word = dat %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))),
  per_class_summary = per_class,
  overall_r2 = data.frame(
    model = c("log_p only", "class only", "log_p + class"),
    r2 = c(summary(m_log_p)$r.squared,
           summary(m_class)$r.squared,
           summary(m_full)$r.squared))
)

# Write the per-word table as the main CSV (most useful for the paper)
write.csv(results$per_word,
          file.path(OUT_FIGS, "psi_freq_regression_per_word.csv"),
          row.names = FALSE)
write.csv(results$per_class_summary,
          file.path(OUT_FIGS, "psi_freq_regression_per_class.csv"),
          row.names = FALSE)
write.csv(results$overall_r2,
          file.path(OUT_FIGS, "psi_freq_regression_r2.csv"),
          row.names = FALSE)

cat(sprintf("\nWrote 3 CSVs to %s\n", OUT_FIGS))

# ---- Plot ---- #
p <- ggplot(dat, aes(x = log_p, y = psi_median, color = lexical_class)) +
  geom_pointrange(aes(ymin = psi_q025, ymax = psi_q975),
                  size = 0.2, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  facet_wrap(~lexical_class) +
  labs(x = "log p_j (CHILDES American English)",
       y = "psi_j posterior median (M1: time only)",
       title = "Per-word log-threshold vs frequency, by lexical class",
       subtitle = sprintf("R^2 = %.2f (log_p alone), %.2f (log_p + class), n = %d words",
                          summary(m_log_p)$r.squared,
                          summary(m_full)$r.squared,
                          nrow(dat))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUT_FIGS, "psi_freq_regression.png"),
       p, width = 9, height = 7, dpi = 200)
cat(sprintf("Wrote %s\n", file.path(OUT_FIGS, "psi_freq_regression.png")))
