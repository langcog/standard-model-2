## Chang & Bergen 2022 sigmoid-slope comparison.
##
## Both kids' and LMs' word-acquisition curves are sigmoidal in
## log-experience. The sigmoid slope on logit-per-log-experience is
## directly comparable across the two:
##
##   Kids (per-child, per-word): logit P(produces) = ... + kappa_i log t
##                               slope on ln(t) = kappa_i
##   LMs (Chang & Bergen, per-word): logit P(learned) = (x - xmid)/scal
##                                   x = log10(steps); slope on log10 = 1/scal
##
## To compare on the same axis, convert LM slopes to natural-log units:
##   slope_nat = (1/scal) / ln(10)
##
## Output:
##   outputs/figs/longitudinal/chang_bergen_slope_comparison.png
##   outputs/figs/longitudinal/chang_bergen_slope_summary.csv

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(posterior)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")
DATA_DIR <- "data/chang_bergen_2022"

# ---- Load Chang & Bergen sigmoid fits ----------------------------
LMS <- c("bert", "bilstm", "gpt2", "lstm")
lm_data <- lapply(LMS, function(lm) {
  d <- read.delim(file.path(DATA_DIR, sprintf("%s_sigmoids.txt", lm)),
                  stringsAsFactors = FALSE)
  d$model <- lm
  d
}) |> bind_rows()

cat(sprintf("Loaded %d (LM, word) fits.\n", nrow(lm_data)))
cat("Columns:", paste(names(lm_data), collapse = ", "), "\n\n")

# ---- Compute per-(LM, word) slope on natural-log experience ------
##
## Filter out degenerate fits: ParamScale very small (<0.01) or very
## large (>10) suggests the sigmoid didn't fit a meaningful transition
## (e.g., word never learned, or fit went to a numerical edge).
## Also drop where the surprisal range is tiny (max - min < 1.0 nat).
lm_data <- lm_data |>
  mutate(
    surprisal_range = ParamUpper - ParamLower,
    slope_log10 = 1 / ParamScale,
    slope_natural = slope_log10 / log(10)   # logit per natural-log experience
  ) |>
  filter(is.finite(slope_natural),
         ParamScale > 0.01, ParamScale < 10,
         surprisal_range > 1.0)

cat(sprintf("After filtering: %d (LM, word) fits.\n", nrow(lm_data)))

# Per-LM summary
lm_summary <- lm_data |>
  group_by(model) |>
  summarise(
    n_words = n(),
    median_slope = median(slope_natural),
    q025 = quantile(slope_natural, 0.025),
    q975 = quantile(slope_natural, 0.975),
    mean_slope = mean(slope_natural),
    sd_slope = sd(slope_natural)
  )
cat("\n=== Per-LM slope summary (logit per natural-log experience) ===\n")
print(as.data.frame(lm_summary), digits = 3)

# ---- Load our kid kappa distribution -----------------------------
##
## kappa_i = 1 + delta + zeta_i. Sample N_KIDS_DRAW kids from posterior
## of M_best (English long_no_freq_slopes) and Norwegian for comparison.
draws_eng <- readRDS(file.path(PATHS$fits_dir, "summaries",
                                "long_no_freq_slopes.draws.rds"))
draws_nor <- readRDS(file.path(PATHS$fits_dir, "summaries",
                                "long_no_freq_slopes_norwegian.draws.rds"))

as_num <- function(x) as.numeric(unlist(x))
sample_kappa_dist <- function(draws, label, N_kids = 5000) {
  delta <- as_num(draws$delta)
  sigma_zeta <- as_num(draws$sigma_zeta)
  set.seed(2026)
  # Sample (delta_d, sigma_zeta_d) from posterior, then per-draw sample kids
  draw_idx <- sample.int(length(delta), size = 50)
  out <- lapply(draw_idx, function(d) {
    n <- ceiling(N_kids / length(draw_idx))
    kappa_draw <- 1 + delta[d] + rnorm(n, mean = 0, sd = sigma_zeta[d])
    data.frame(model = label, slope_natural = kappa_draw)
  }) |> bind_rows()
  out
}

kid_eng <- sample_kappa_dist(draws_eng, "Kids (English, M_best)")
kid_nor <- sample_kappa_dist(draws_nor, "Kids (Norwegian, M_best)")

cat(sprintf("\nEnglish kids: median kappa = %.2f, IQR = [%.2f, %.2f]\n",
            median(kid_eng$slope_natural),
            quantile(kid_eng$slope_natural, 0.25),
            quantile(kid_eng$slope_natural, 0.75)))
cat(sprintf("Norwegian kids: median kappa = %.2f, IQR = [%.2f, %.2f]\n",
            median(kid_nor$slope_natural),
            quantile(kid_nor$slope_natural, 0.25),
            quantile(kid_nor$slope_natural, 0.75)))

# ---- Combine for plotting ----------------------------------------
lm_plot_data <- lm_data |>
  mutate(model = factor(model,
                        levels = c("bert", "gpt2", "bilstm", "lstm"),
                        labels = c("BERT", "GPT-2", "BiLSTM", "LSTM"))) |>
  select(model, slope_natural)
lm_plot_data$category <- "LMs (Chang & Bergen 2022)"

kid_plot_data <- bind_rows(kid_eng, kid_nor)
kid_plot_data$model <- factor(kid_plot_data$model)
kid_plot_data$category <- "Children (this work)"

all_data <- bind_rows(
  lm_plot_data |> rename(label = model),
  kid_plot_data |> rename(label = model)
)
all_data$label <- factor(all_data$label,
                          levels = c("Kids (English, M_best)",
                                     "Kids (Norwegian, M_best)",
                                     "BERT", "GPT-2", "BiLSTM", "LSTM"))
all_data$category <- factor(all_data$category,
                             levels = c("Children (this work)",
                                        "LMs (Chang & Bergen 2022)"))

# ---- Plot: density distributions on common axis ------------------
write.csv(all_data |> group_by(label) |>
            summarise(median = median(slope_natural),
                      q025 = quantile(slope_natural, 0.025),
                      q975 = quantile(slope_natural, 0.975),
                      n = n()),
          file.path(OUT_DIR, "chang_bergen_slope_summary.csv"),
          row.names = FALSE)

p <- ggplot(all_data, aes(slope_natural, fill = label, colour = label)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55",
             linewidth = 0.4) +
  annotate("text", x = 1.05, y = 0.5, label = "kappa = 1\n(unit accumulator)",
           hjust = 0, vjust = 0.5, size = 3.0, colour = "grey45") +
  scale_fill_manual(
    values = c("Kids (English, M_best)" = "#c41e37",
               "Kids (Norwegian, M_best)" = "#fb9a99",
               "BERT" = "#1f78b4", "GPT-2" = "#33a02c",
               "BiLSTM" = "#a6cee3", "LSTM" = "#b2df8a")
  ) +
  scale_colour_manual(
    values = c("Kids (English, M_best)" = "#c41e37",
               "Kids (Norwegian, M_best)" = "#fb9a99",
               "BERT" = "#1f78b4", "GPT-2" = "#33a02c",
               "BiLSTM" = "#a6cee3", "LSTM" = "#b2df8a")
  ) +
  coord_cartesian(xlim = c(0, 22)) +
  labs(
    x = "Slope on log(experience): logit per natural-log unit",
    y = "Density",
    fill = NULL, colour = NULL,
    title = "Per-instance scaling slopes: children vs.\ LLMs",
    subtitle = "Sigmoid slope of P(word acquired) vs.\ log(experience). Both populations on the same axis: kids (per-child kappa_i in natural-log time); LMs (per-word 1/ParamScale in log10-steps, converted to natural log)."
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"))

ggsave(file.path(OUT_DIR, "chang_bergen_slope_comparison.png"), p,
       width = 10, height = 5.8, dpi = 150)
cat(sprintf("\nWrote: %s\n",
            file.path(OUT_DIR, "chang_bergen_slope_comparison.png")))

cat("\nDone.\n")
