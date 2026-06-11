## Feng et al. (2026) CHILDES-trained GPT-2 vs. children vs.
## Chang & Bergen (2022) per-word sigmoid slopes.
##
## Mirrors `studies/llm/chang_bergen_comparison.R`. Reads:
##   data/chang_bergen_2022/{bert,gpt2,bilstm,lstm}_sigmoids.txt   (BookCorpus-trained)
##   fits/llm/sigmoids/<model>_sigmoids.txt                            (CHILDES-trained, our pipeline)
##   fits/summaries/long_no_freq_slopes.draws.rds                   (kid kappa, English)
##   fits/summaries/long_no_freq_slopes_norwegian.draws.rds         (kid kappa, Norwegian)
##
## Produces:
##   figs/longitudinal/feng_slope_comparison.png            (kids + Feng only)
##   figs/longitudinal/feng_chang_bergen_slope_comparison.png  (kids + Feng + C&B)
##   figs/longitudinal/feng_slope_summary.csv

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(posterior)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")
CB_DIR  <- "data/chang_bergen_2022"
FENG_DIR <- "fits/llm/sigmoids"

# ---- Load Chang & Bergen sigmoid fits ----------------------------
CB_MODELS <- c("bert", "bilstm", "gpt2", "lstm")
cb_data <- lapply(CB_MODELS, function(lm) {
  d <- read.delim(file.path(CB_DIR, sprintf("%s_sigmoids.txt", lm)),
                  stringsAsFactors = FALSE)
  d$model <- lm
  d$source <- "C&B BookCorpus"
  d
}) |> bind_rows()

# ---- Load Feng et al. (this work) sigmoid fits -------------------
feng_files <- list.files(FENG_DIR, pattern = "_sigmoids\\.txt$", full.names = TRUE)
if (length(feng_files) == 0) {
  stop("No Feng sigmoid files under ", FENG_DIR,
       ". Run studies/llm/fit_per_word_sigmoid.py first.")
}
feng_data <- lapply(feng_files, function(p) {
  d <- read.delim(p, stringsAsFactors = FALSE)
  d$model <- sub("_sigmoids\\.txt$", "", basename(p))
  d$source <- "Feng CHILDES"
  d
}) |> bind_rows()

cat(sprintf("C&B fits: %d, Feng fits: %d\n", nrow(cb_data), nrow(feng_data)))

# ---- Slope on natural-log experience (logit per nat-log step) ----
compute_slopes <- function(d) {
  d |>
    mutate(
      surprisal_range = ParamUpper - ParamLower,
      slope_log10 = 1 / ParamScale,
      slope_natural = slope_log10 / log(10)
    ) |>
    filter(is.finite(slope_natural),
           ParamScale > 0.01, ParamScale < 10,
           surprisal_range > 1.0)
}
cb_data <- compute_slopes(cb_data)
feng_data <- compute_slopes(feng_data)

cat(sprintf("After filtering — C&B: %d, Feng: %d\n",
            nrow(cb_data), nrow(feng_data)))

# Per-model summary
lm_summary <- function(d, lab) {
  d |>
    group_by(model) |>
    summarise(
      source = first(source),
      n_words = n(),
      median_slope = median(slope_natural),
      q025 = quantile(slope_natural, 0.025),
      q975 = quantile(slope_natural, 0.975),
      mean_slope = mean(slope_natural),
      sd_slope = sd(slope_natural)
    )
}
cat("\n=== Chang & Bergen models ===\n")
print(as.data.frame(lm_summary(cb_data, "cb")), digits = 3)
cat("\n=== Feng CHILDES models ===\n")
print(as.data.frame(lm_summary(feng_data, "feng")), digits = 3)

# ---- Load our kid kappa distribution -----------------------------
# Prefer the posterior draws; fall back to the published kappa_pop/sigma_zeta
# summary from outputs/chang_bergen_derivation.tex when draws aren't present.

draws_eng_path <- file.path(PATHS$fits_dir, "summaries",
                             "long_no_freq_slopes.draws.rds")
draws_nor_path <- file.path(PATHS$fits_dir, "summaries",
                             "long_no_freq_slopes_norwegian.draws.rds")

as_num <- function(x) as.numeric(unlist(x))

sample_kappa_from_draws <- function(draws, label, N_kids = 5000) {
  delta <- as_num(draws$delta)
  sigma_zeta <- as_num(draws$sigma_zeta)
  set.seed(2026)
  draw_idx <- sample.int(length(delta), size = 50)
  lapply(draw_idx, function(d) {
    n <- ceiling(N_kids / length(draw_idx))
    kappa_draw <- 1 + delta[d] + rnorm(n, mean = 0, sd = sigma_zeta[d])
    data.frame(model = label, slope_natural = kappa_draw)
  }) |> bind_rows()
}

sample_kappa_from_summary <- function(kappa_pop, sigma_zeta, label, N_kids = 5000) {
  set.seed(2026)
  data.frame(model = label,
             slope_natural = kappa_pop + rnorm(N_kids, 0, sigma_zeta))
}

if (file.exists(draws_eng_path)) {
  kid_eng <- sample_kappa_from_draws(readRDS(draws_eng_path),
                                      "Kids (English, M_best)")
} else {
  message("NOTE: ", draws_eng_path, " not found; using published M_best summary",
          " (kappa_pop=10.3, sigma_zeta=3.47) instead of posterior draws.")
  kid_eng <- sample_kappa_from_summary(10.3, 3.47, "Kids (English, M_best)")
}
if (file.exists(draws_nor_path)) {
  kid_nor <- sample_kappa_from_draws(readRDS(draws_nor_path),
                                      "Kids (Norwegian, M_best)")
} else {
  message("NOTE: ", draws_nor_path, " not found; using published M_best summary",
          " (kappa_pop=12.5, sigma_zeta=3.50) instead of posterior draws.")
  kid_nor <- sample_kappa_from_summary(12.5, 3.50, "Kids (Norwegian, M_best)")
}

# ---- Common plotting helper --------------------------------------
pretty_lm <- function(d) {
  d |>
    mutate(label = case_when(
      source == "C&B BookCorpus" & model == "bert"   ~ "BERT (BookCorpus)",
      source == "C&B BookCorpus" & model == "gpt2"   ~ "GPT-2 (BookCorpus)",
      source == "C&B BookCorpus" & model == "bilstm" ~ "BiLSTM (BookCorpus)",
      source == "C&B BookCorpus" & model == "lstm"   ~ "LSTM (BookCorpus)",
      source == "Feng CHILDES"                       ~ sprintf("%s (CHILDES)", model),
      TRUE ~ paste(source, model)
    )) |>
    select(label, source, slope_natural)
}

cb_plot   <- pretty_lm(cb_data)
feng_plot <- pretty_lm(feng_data)

kid_plot <- bind_rows(kid_eng, kid_nor) |>
  rename(label = model) |>
  mutate(source = "Children")

# ---- Output 1: Feng-only (kids + CHILDES LMs) --------------------
feng_only <- bind_rows(
  kid_plot,
  feng_plot
)
feng_only$source <- factor(feng_only$source,
                            levels = c("Children", "Feng CHILDES"))

p_feng <- ggplot(feng_only, aes(slope_natural, fill = label, colour = label)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55",
             linewidth = 0.4) +
  annotate("text", x = 1.05, y = 0.5, label = "kappa = 1\n(unit accumulator)",
           hjust = 0, vjust = 0.5, size = 3.0, colour = "grey45") +
  coord_cartesian(xlim = c(0, 22)) +
  labs(
    x = "Slope on log(experience): logit per natural-log unit",
    y = "Density",
    fill = NULL, colour = NULL,
    title = "Per-instance scaling slopes: children vs. CHILDES-trained GPT-2",
    subtitle = "LMs match the input distribution children receive. Does the kid-vs-LM gap persist?"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"))
ggsave(file.path(OUT_DIR, "feng_slope_comparison.png"), p_feng,
       width = 9, height = 4.5, dpi = 200)

# ---- Output 2: combined (kids + Feng + C&B) ----------------------
combined <- bind_rows(kid_plot, feng_plot, cb_plot)
combined$source <- factor(combined$source,
                           levels = c("Children", "Feng CHILDES", "C&B BookCorpus"))
p_combined <- ggplot(combined, aes(slope_natural, fill = source, colour = source)) +
  geom_density(alpha = 0.4, linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55",
             linewidth = 0.4) +
  annotate("text", x = 1.05, y = 0.5, label = "kappa = 1\n(unit accumulator)",
           hjust = 0, vjust = 0.5, size = 3.0, colour = "grey45") +
  scale_fill_manual(values = c("Children" = "#c41e37",
                                "Feng CHILDES" = "#1f78b4",
                                "C&B BookCorpus" = "#33a02c")) +
  scale_colour_manual(values = c("Children" = "#c41e37",
                                  "Feng CHILDES" = "#1f78b4",
                                  "C&B BookCorpus" = "#33a02c")) +
  coord_cartesian(xlim = c(0, 22)) +
  labs(
    x = "Slope on log(experience): logit per natural-log unit",
    y = "Density",
    fill = NULL, colour = NULL,
    title = "Per-instance scaling slopes: children, CHILDES LMs, BookCorpus LMs",
    subtitle = "Per-word 1/ParamScale (logit per natural-log experience) vs. per-child kappa_i."
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"))
ggsave(file.path(OUT_DIR, "feng_chang_bergen_slope_comparison.png"), p_combined,
       width = 9, height = 4.5, dpi = 200)

# ---- Summary CSV ------------------------------------------------
summary_csv <- combined |>
  group_by(label, source) |>
  summarise(median = median(slope_natural),
            q25 = quantile(slope_natural, 0.25),
            q75 = quantile(slope_natural, 0.75),
            q025 = quantile(slope_natural, 0.025),
            q975 = quantile(slope_natural, 0.975),
            n = n(),
            .groups = "drop")
write.csv(summary_csv, file.path(OUT_DIR, "feng_slope_summary.csv"),
          row.names = FALSE)
cat("\nWrote:\n  feng_slope_comparison.png\n  feng_chang_bergen_slope_comparison.png\n  feng_slope_summary.csv\n")
