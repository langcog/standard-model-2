## Per-word "exposures to 50% threshold" plot, derived from the EN
## M_best fit. Each dot is one word in the bundle, positioned at
##   x = predicted age at which a typical kid reaches 50% production,
##   y = cumulative exposures of THAT word by that age.
##
## Solves for t_j from the linear predictor under a typical kid:
##   eta = xi_typ + log H + kappa_typ * log(t/a0) - delta_j = 0
##   => t_j = a0 * exp((delta_j - log H - xi_typ) / kappa_typ)
## Cumulative word-j exposures:
##   N_j = exp(xi_typ) * H * t_j * p_j         (rate in tokens/hr; H in hr/mo)
##
## Two reference lines:
##   1. Pure accumulator (kappa = 1): same math with kappa fixed at 1.
##      Shows what a non-efficiency-gaining learner would need.
##   2. Constant-threshold (rough): a fixed N per word, age would just
##      track inversely with frequency. (Plotted only if we want.)
##
## Output: outputs/figs/longitudinal/exposure_to_learn_EN.png

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(ggrepel)
})

OUT_DIR    <- file.path(PATHS$figs_dir, "longitudinal")
SUMMARIES  <- file.path(PATHS$fits_dir, "summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BUNDLE_PATH <- file.path(PATHS$fits_dir, "long_subset_data.rds")
TAG         <- "long_no_freq_slopes"
SEED        <- 20260523

bundle <- readRDS(BUNDLE_PATH)
sd_b   <- bundle$stan_data
word_info <- bundle$word_info
class_lev <- bundle$class_levels

log_H_const <- sd_b$log_H       # 5.9, hours/month exposure
a0          <- sd_b$a0          # reference age in months
mu_r        <- sd_b$mu_r        # 7.34 typical log r (tokens/hr)
sigma_r     <- sd_b$sigma_r
cat(sprintf("Bundle: log_H=%.3f a0=%g mu_r=%.3f sigma_r=%.3f\n",
            log_H_const, a0, mu_r, sigma_r))

## ---- Posterior medians for the typical kid -----------------------
draws <- readRDS(file.path(SUMMARIES, paste0(TAG, ".draws.rds")))
draws <- as.data.frame(draws)
delta_post   <- median(draws$delta)
sigma_alpha_post <- median(draws$sigma_alpha)
cat(sprintf("Posterior medians: delta=%.3f sigma_alpha=%.3f\n",
            delta_post, sigma_alpha_post))
xi_typ    <- mu_r                  # typical kid: log_alpha = 0, log_r = mu_r
kappa_typ <- 1 + delta_post        # typical kappa = 1 + delta + 0
r_typ     <- exp(xi_typ)           # tokens / hr

## ---- Item difficulties -------------------------------------------
psi <- read_csv(file.path(SUMMARIES, paste0(TAG, "_psi.csv")),
                show_col_types = FALSE) |>
  rename(delta_j = delta_j_median)
items <- word_info |>
  mutate(jj = row_number()) |>
  left_join(psi |> select(jj, delta_j), by = "jj") |>
  filter(!is.na(delta_j), prob > 0) |>
  mutate(lexical_class = factor(class_lev[cc], levels = class_lev))
cat(sprintf("Items: %d with delta_j and prob > 0\n", nrow(items)))

# Drop CDI items whose CHILDES frequency is suspiciously low. ~54 items
# got the explicit no-match floor (1.16e-7); another ~20 items have
# tokenization mismatches (CDI "teddybear" vs CHILDES "teddy bear",
# CDI "grrr" vs CHILDES "grr", multi-word phrases like "a lot",
# sound-effect-like items, etc.) that produced near-zero frequencies.
# All of these would form a low-N artifact in the plot without telling
# us anything about real word-learning. Threshold at 1e-5 (any word
# heard < 1 in 100,000 tokens in CHILDES is almost certainly a match
# bug rather than a truly rare word a kid is acquiring).
FREQ_MIN <- 1e-5
n_drop <- sum(items$prob < FREQ_MIN)
cat(sprintf("Dropped %d items with prob < %g (CHILDES no-match / tokenization mismatch)\n",
            n_drop, FREQ_MIN))
items <- items |> filter(prob >= FREQ_MIN)

## ---- Predicted age-of-50% and exposure count per word ------------
solve_t <- function(delta_j, kappa, xi) {
  a0 * exp((delta_j - log_H_const - xi) / kappa)
}
items <- items |>
  mutate(
    t_50_mbest = solve_t(delta_j, kappa_typ, xi_typ),
    t_50_pure  = solve_t(delta_j, 1, xi_typ),
    N_word_mbest = r_typ * exp(log_H_const) * t_50_mbest * prob,
    N_word_pure  = r_typ * exp(log_H_const) * t_50_pure  * prob,
    N_total_mbest = r_typ * exp(log_H_const) * t_50_mbest
  )

## ---- Print a few representative entries --------------------------
ex <- items |> arrange(delta_j) |>
  slice(c(1, 3, round(n() * 0.25), round(n() * 0.5),
          round(n() * 0.75), n() - 2, n())) |>
  select(item, lexical_class, delta_j, prob, t_50_mbest, N_word_mbest)
cat("\nExample words (easy -> hard):\n")
print(ex)

## ---- Label-pick: top + bottom item per class on log_N ------------
# Log-N because the y-axis is log scale; choosing extremes in log_N
# gives a balanced set across classes regardless of count magnitude.
labels_df <- items |>
  mutate(log_N = log10(N_word_mbest)) |>
  group_by(lexical_class) |>
  arrange(log_N) |>
  slice(c(1, n())) |>
  ungroup()

## ---- Plot: word j exposures by age of 50% ------------------------
p <- ggplot(items, aes(t_50_mbest, N_word_mbest, colour = lexical_class)) +
  geom_point(alpha = 0.55, size = 1.1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  geom_label_repel(data = labels_df,
                    aes(label = item),
                    size = 2.8, alpha = 0.95,
                    show.legend = FALSE,
                    box.padding = 0.4, max.overlaps = Inf,
                    min.segment.length = 0) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  scale_x_continuous(limits = c(10, 40)) +
  scale_colour_brewer(palette = "Dark2", name = "Class") +
  labs(x = "Age (months) when word is produced at 50%",
       y = "Cumulative exposures of THIS word at that age",
       title = "How many times has a typical kid heard each word by the time they produce it?",
       subtitle = sprintf(
         "EN M_best at I=500 (%d items; CHILDES no-match floor items excluded). Typical kid: log r=%.2f, kappa=%.2f. Lines: per-class lm fits.",
         nrow(items), mu_r, kappa_typ)) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        legend.position = "bottom")

out_png <- file.path(OUT_DIR, "exposure_to_learn_EN.png")
ggsave(out_png, p, width = 10, height = 6, dpi = 200)
cat(sprintf("\nWrote: %s\n", out_png))

# (Pure-accumulator side-by-side comparison removed -- the headline
# plot tells the story on its own.)
