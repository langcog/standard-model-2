## Pre-model input-vs-uptake check for the four IO datasets.
##
## Motivation: the anchored fits return pi_alpha ~ 0.97-0.99 for the
## sparse-LENA datasets (AM2018, FMW2013), which is suspicious. ID
## matching is verified correct and the raw input genuinely varies
## (per-kid log-AWC SD 0.37-0.58). The high pi_alpha is a sparse-
## replication artifact: with 1-2 recordings/kid the model can't pin
## sigma_r, so it collapses. This script asks the question MODEL-FREE:
## does a child's observed input rate predict their vocabulary
## trajectory?
##
## For each dataset: bin kids by their mean observed log input rate
## (terciles), then plot vocab (proportion of the form's items, to
## handle WG/WS denominators) vs age coloured by input tercile, with a
## per-tercile smoother. If higher-input kids sit higher, input
## predicts uptake in the raw data even where the model can't claim it.
##
## Output: outputs/figs/io/input_uptake_premodel.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

OUT_PNG <- file.path(PATHS$figs_dir, "io", "input_uptake_premodel.png")
dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)

# bundle, id column, and where the per-recording input lives
SPECS <- list(
  BabyView  = list(bundle = "babyview_subset_data.rds",  id = "subject_id",
                    rec = "videos"),
  SEEDLingS = list(bundle = "seedlings_subset_data.rds",  id = "subject_id",
                    rec = "recordings"),
  AM2018    = list(bundle = "io_am2018_subset_data.rds",  id = "subject_id",
                    rec = "recordings"),
  FMW2013   = list(bundle = "io_fmw2013_subset_data.rds", id = "subject_id",
                    rec = "recordings")
)

one <- function(label, spec) {
  b <- readRDS(file.path(PATHS$fits_dir, spec$bundle))
  df <- b$df
  # The df carries an integer child index `ii`; recordings carry the
  # matching `child_ii`. Join on the integer index (robust — BabyView's
  # videos$subject_id "S00360001" differs from df$subject_id formatting).
  stopifnot("ii" %in% names(df))

  rec <- b[[spec$rec]]
  stopifnot("child_ii" %in% names(rec))
  kid_in <- rec |> group_by(ii = child_ii) |>
    summarise(mean_log_r = mean(log_r_obs, na.rm = TRUE), .groups = "drop")

  per_admin <- df |>
    distinct(ii, age, form, item, .keep_all = TRUE) |>
    group_by(ii, age, form) |>
    summarise(vocab = sum(produces, na.rm = TRUE), n = n(),
              prop = vocab / n, .groups = "drop")

  d <- per_admin |> inner_join(kid_in, by = "ii")
  nq <- if (n_distinct(d$ii) >= 30) 3 else 2
  brk <- quantile(kid_in$mean_log_r, seq(0, 1, length.out = nq + 1), na.rm = TRUE)
  d$input_bin <- cut(d$mean_log_r, breaks = brk, include.lowest = TRUE,
                      labels = if (nq == 3) c("low","mid","high")
                               else c("low","high"))
  d$dataset <- label
  d
}

dat <- bind_rows(Map(one, names(SPECS), SPECS))
dat$dataset <- factor(dat$dataset, levels = names(SPECS))

# report the raw input-bin vocab gap per dataset (model-free)
cat("Per-dataset: mean vocab proportion by input tercile (model-free)\n")
print(dat |> group_by(dataset, input_bin) |>
        summarise(n_kids = n_distinct(ii), mean_prop = mean(prop),
                  .groups = "drop"))

BIN_COL <- c(low = "#2c7bb6", mid = "#fdae61", high = "#d7191c")

p <- ggplot(dat, aes(age, prop, colour = input_bin, fill = input_bin)) +
  geom_point(alpha = 0.30, size = 0.9) +
  geom_smooth(method = "loess", se = TRUE, span = 1.0, linewidth = 1.0,
              alpha = 0.15) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = BIN_COL, name = "Observed input") +
  scale_fill_manual(values = BIN_COL, name = "Observed input") +
  labs(x = "Age (months)",
       y = "Proportion of form's items produced",
       title = "Pre-model: does observed input rate predict vocabulary uptake?",
       subtitle = "Kids binned by mean observed log input (LENA AWC/hr or head-cam tokens). If input predicts uptake, high (red) sits above low (blue).") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

ggsave(OUT_PNG, p, width = 14, height = 4.5, dpi = 200)
cat(sprintf("\nWrote %s\n", OUT_PNG))
