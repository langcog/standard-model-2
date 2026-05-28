## Aggregate per-fit ladder summaries into one table + figure.
##
## Usage:   Rscript model/scripts/aggregate_ladder.R
##
## Inputs:  fits/glmer_ladder/summary_<lang>_<model>.csv (one per fit)
## Outputs:
##   outputs/glmer_ladder/ladder_summary.csv              -- combined long table
##   outputs/glmer_ladder/ladder_deltas.csv               -- ΔAIC vs worst per language
##   outputs/figs/glmer_ladder/deltaAIC.png-- panel figure

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

OUT_SUMM <- file.path(PATHS$outputs_dir, "glmer_ladder/ladder_summary.csv")
OUT_DELT <- file.path(PATHS$outputs_dir, "glmer_ladder/ladder_deltas.csv")
OUT_PNG  <- file.path(PATHS$figs_dir, "glmer_ladder", "deltaAIC.png")
dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)

ladder_dir <- file.path(PATHS$fits_dir, "glmer_ladder")
files <- list.files(ladder_dir, pattern = "^summary_.*\\.csv$",
                     full.names = TRUE)
cat(sprintf("Found %d per-fit summary CSVs\n", length(files)))
if (length(files) == 0) stop("no summary files — did you sync from Sherlock?")

summ <- bind_rows(lapply(files, read_csv, show_col_types = FALSE))

## Model factor with the ladder ordering
MODEL_LEVELS <- c("A", "B_lin", "B_log", "C_lin", "C_log", "D_lin", "D_log")
MODEL_LABELS <- c(
  A     = "A: pure accumulator",
  B_lin = "B_lin: + κ on age",
  B_log = "B_log: + κ on log_age",
  C_lin = "C_lin: + (1|child)",
  C_log = "C_log: + (1|child)",
  D_lin = "D_lin: + (1+age|child)",
  D_log = "D_log: + (1+log_age|child)"
)
summ <- summ |>
  mutate(model = factor(model, levels = MODEL_LEVELS))

write_csv(summ, OUT_SUMM)
cat(sprintf("Wrote %s\n", OUT_SUMM))

## ΔAIC relative to the WORST model in each language (so all bars are ≤ 0
## with the best model at the most-negative; or flip sign for "better is up").
## We use "AIC improvement vs worst" = AIC_worst − AIC, so bigger = better.
deltas <- summ |>
  group_by(language) |>
  mutate(
    AIC_worst       = max(AIC),
    AIC_improvement = AIC_worst - AIC,
    AIC_best        = min(AIC),
    delta_from_best = AIC - AIC_best
  ) |>
  ungroup() |>
  arrange(language, model)

write_csv(deltas, OUT_DELT)
cat(sprintf("Wrote %s\n", OUT_DELT))

## ---- Figure: AIC improvement vs worst, per language ----
## Each row = one language; bars colored by lin/log family.
deltas <- deltas |>
  mutate(family = case_when(
    model == "A" ~ "Pure accumulator",
    grepl("lin$", model) ~ "Linear age",
    grepl("log$", model) ~ "Log age",
    TRUE ~ "?"
  ))

p <- ggplot(deltas,
             aes(x = model,
                 y = AIC_improvement,
                 fill = family)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(delta_from_best == 0,
                                "best",
                                sprintf("Δ=%+.0f", delta_from_best))),
             vjust = -0.3, size = 2.8) +
  facet_wrap(~ language, scales = "free_y") +
  scale_fill_manual(values = c(`Pure accumulator` = "grey60",
                                `Linear age`       = "#d7191c",
                                `Log age`          = "#1f78b4")) +
  scale_x_discrete(labels = MODEL_LABELS) +
  labs(x = NULL,
       y = "AIC improvement vs worst model in language",
       title = "Model-ladder comparison: longitudinal CDI data, all languages with ≥100 longitudinal kids",
       subtitle = "Higher = better fit. Linear-age and log-age branches diverge at B; log-age wins consistently.",
       fill = NULL,
       caption = "All models fit via glmer (nAGQ=0). N per language varies; see ladder_summary.csv.") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        legend.position = "top",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        strip.text = element_text(face = "bold"))

ggsave(OUT_PNG, p, width = 13, height = 9, dpi = 200)
cat(sprintf("Wrote %s\n", OUT_PNG))

## ---- Quick text summary ----
cat("\n=== Summary table ===\n")
print(deltas |>
        select(language, model, AIC, delta_from_best, n_kids, n_obs) |>
        arrange(language, model), n = 100)

cat("\n=== Best model per language ===\n")
print(deltas |>
        filter(delta_from_best == 0) |>
        select(language, model, AIC, n_kids, n_obs))
