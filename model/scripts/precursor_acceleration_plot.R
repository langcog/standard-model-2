## Model-free acceleration figure for slide 3a.
##
## The canonical Wordbank vocabulary-vs-age plot with smoothed quantile
## lines, in the visual style Mike highlighted (vocab-eng-prod-1):
## individual children as points, gcrq quantile curves at 10/25/50/75/90.
##
## Key visual point: vocabulary growth is super-linear in age — the
## quantile curves are concave-up across the productive window. No
## model fitting; just empirical Wordbank data.
##
## Output: outputs/figs/longitudinal/precursor_acceleration.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(quantregGrowth); library(wordbankr)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")

# Pull English (American) WS — the canonical sample
cat("Pulling English (American) WS data...\n")
raw <- get_instrument_data(language = "English (American)", form = "WS",
                           administration_info = TRUE, item_info = TRUE)

# Total vocabulary per (child, age)
vocab <- raw |>
  filter(item_kind == "word") |>
  group_by(child_id, age, data_id) |>
  summarise(vocab = sum(produces, na.rm = TRUE), .groups = "drop") |>
  filter(age >= 16, age <= 30)

cat(sprintf("N admins: %d, ages %d-%d\n",
            nrow(vocab), min(vocab$age), max(vocab$age)))

# Quantile regression
TAUS <- c(0.1, 0.25, 0.5, 0.75, 0.9)
q_fit <- gcrq(vocab ~ ps(age, monotone = 1, lambda = 1000),
              tau = TAUS, data = vocab)
age_grid <- seq(16, 30, by = 0.25)
preds <- predict(q_fit, newdata = data.frame(age = age_grid))
preds_df <- as.data.frame(preds)
names(preds_df) <- as.character(TAUS)
preds_df$age <- age_grid
preds_long <- pivot_longer(preds_df, cols = -age, names_to = "tau",
                            values_to = "vocab")

p <- ggplot(vocab, aes(age, vocab)) +
  geom_jitter(width = 0.3, alpha = 0.18, size = 0.6, colour = "grey25") +
  geom_line(data = preds_long, aes(age, vocab, colour = tau, group = tau),
            linewidth = 1.1) +
  scale_colour_manual(
    values = c("0.1" = "#1f78b4", "0.25" = "#a6cee3",
               "0.5" = "#33a02c", "0.75" = "#fdbf6f", "0.9" = "#e31a1c"),
    name = "Percentile"
  ) +
  scale_x_continuous(breaks = c(16, 18, 20, 22, 24, 26, 28, 30)) +
  labs(x = "Age (months)",
       y = "Productive vocabulary (number of words)",
       title = "Acceleration: vocabulary grows super-linearly with age",
       subtitle = "Wordbank English (American) WS form. Each point = one child; lines = 10/25/50/75/90 percentiles via quantile regression.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey25"),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "precursor_acceleration.png"), p,
       width = 8, height = 5.5, dpi = 200)
cat(sprintf("Wrote: %s\n",
            file.path(OUT_DIR, "precursor_acceleration.png")))
