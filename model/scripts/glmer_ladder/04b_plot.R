## glmer ladder — PLOT step (fast).
##
## Reads the cache written by 04a_simulate.R and builds the figures.
## Iterate on aesthetics here without re-running the simulation.
##
## Input:   outputs/glmer_ladder/sim_cache.rds
## Outputs: outputs/figs/longitudinal/glmer_ladder_main.png  (4 langs)
##          outputs/figs/longitudinal/glmer_ladder_mega.png  (6 langs)

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})

CACHE   <- file.path(PATHS$outputs_dir, "glmer_ladder", "sim_cache.rds")
FIG_DIR <- file.path(PATHS$figs_dir, "longitudinal")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(CACHE)) {
  stop("missing sim cache — run 04a_simulate.R first: ", CACHE)
}
cache <- readRDS(CACHE)
qtiles      <- cache$qtiles
emp         <- cache$emp
summ        <- cache$summ
LANGS       <- cache$LANGS
LANG_LABELS <- cache$LANG_LABELS
MODELS      <- cache$MODELS
N_SIM_KIDS  <- cache$N_SIM_KIDS

WORDBANK_PALETTE <- c("0.1"  = "#1f78b4",   # dark blue
                       "0.25" = "#a6cee3",   # light blue
                       "0.5"  = "#33a02c",   # green
                       "0.75" = "#fdbf6f",   # gold
                       "0.9"  = "#e31a1c")   # red

qtiles <- qtiles |>
  mutate(qprob_f = factor(qprob, levels = c(0.1, 0.25, 0.5, 0.75, 0.9),
                           labels = c("0.1", "0.25", "0.5", "0.75", "0.9")))

## Friendly column labels for the conceptual-ladder (log-only) figure.
MODEL_LABELS <- c(
  A     = "A: unit accumulator (κ=1)",
  B_lin = "B: + acceleration κ (linear age)",
  B_log = "B: + acceleration κ",
  C_lin = "C: + child intercept (linear age)",
  C_log = "C: + child intercept",
  D_lin = "D: + child slope (linear age)",
  D_log = "D: + child slope"
)

## ---- Build the plot for a given language × model subset --------------
build_ladder_plot <- function(langs_subset, out_png, title,
                                models_subset = MODELS,
                                relabel_models = FALSE,
                                width = 18, height = 16) {
  qt <- qtiles |>
    filter(lang_slug %in% langs_subset, model %in% models_subset) |>
    mutate(language = factor(LANG_LABELS[lang_slug],
                              levels = LANG_LABELS[langs_subset]),
            model    = factor(model, levels = models_subset))
  el <- emp |>
    filter(lang_slug %in% langs_subset) |>
    mutate(language = factor(LANG_LABELS[lang_slug],
                              levels = LANG_LABELS[langs_subset]))
  al <- summ |>
    filter(lang_slug %in% langs_subset) |>
    group_by(language, lang_slug) |>
    mutate(AIC_best = min(AIC), dAIC = AIC - AIC_best) |>   # ΔAIC vs best of ALL 7 models
    ungroup() |>
    filter(model %in% models_subset) |>
    mutate(language = factor(language,
                              levels = LANG_LABELS[langs_subset]),
            model    = factor(model, levels = models_subset),
            label    = ifelse(dAIC == 0, "best",
                              sprintf("Δ%+.0fk", dAIC / 1000)))

  p <- ggplot() +
    geom_line(data = el,
               aes(x = age, y = vocab, group = child_id),
               colour = "grey25", alpha = 0.15, linewidth = 0.25) +
    geom_point(data = el,
                aes(x = age, y = vocab),
                colour = "grey25", alpha = 0.2, size = 0.25) +
    geom_line(data = qt |> filter(model %in% c("C_lin", "C_log",
                                                 "D_lin", "D_log")),
               aes(x = age, y = vocab,
                   colour = qprob_f, group = qprob_f,
                   linewidth = ifelse(as.character(qprob_f) == "0.5",
                                       0.9, 0.55))) +
    geom_line(data = qt |> filter(!model %in% c("C_lin", "C_log",
                                                  "D_lin", "D_log"),
                                    qprob_f == "0.5"),
               aes(x = age, y = vocab, colour = qprob_f),
               linewidth = 0.95) +
    facet_grid(language ~ model, scales = "free_y", space = "fixed",
                labeller = if (relabel_models)
                  labeller(model = MODEL_LABELS) else "label_value") +
    geom_text(data = al,
               aes(x = -Inf, y = Inf, label = label),
               hjust = -0.15, vjust = 1.5, size = 2.6, colour = "grey25") +
    scale_colour_manual(values = WORDBANK_PALETTE, name = "Percentile") +
    scale_linewidth_identity() +
    labs(x = "Age (months)",
         y = "Productive vocabulary",
         title = title,
         subtitle = sprintf("Lines = 10/25/50/75/90 quantiles across %d simulated kids (BLUP bootstrap from each model's child-RE distribution). Grey lines = per-kid longitudinal trajectories. Empirical and predictions both restricted to the largest form per language. Corner label = ΔAIC vs best of all 7 models in that language.",
                            N_SIM_KIDS)) +
    theme_minimal(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey25"),
          strip.text    = element_text(size = 8, face = "bold"),
          legend.position = "top",
          panel.spacing = unit(0.4, "lines"))

  ggsave(out_png, p, width = width, height = height, dpi = 180)
  cat(sprintf("Wrote %s\n", out_png))
}

## Main-text version: 4 well-powered languages, LOG-only ladder
## (A → B → C → D). Linear vs log is a fine-structure AIC difference
## that's visually imperceptible over the CDI age window, so we keep
## the linear columns for the supplement only.
LANGS_MAIN  <- c("english_american", "norwegian",
                 "french_quebecois", "japanese")
MODELS_LOG  <- c("A", "B_log", "C_log", "D_log")
build_ladder_plot(
  LANGS_MAIN,
  file.path(FIG_DIR, "glmer_ladder_main.png"),
  title = "Model ladder: vocabulary growth across four longitudinal CDI samples",
  models_subset = MODELS_LOG,
  relabel_models = TRUE,
  width = 13, height = 11
)

## Supplementary version 1: all 6 languages, LOG-only ladder
build_ladder_plot(
  LANGS,
  file.path(FIG_DIR, "glmer_ladder_supp_log.png"),
  title = "Model ladder (log-age only): six longitudinal CDI samples",
  models_subset = MODELS_LOG,
  relabel_models = TRUE,
  width = 13, height = 16
)

## Supplementary version 2: full 7-model grid (lin + log) — shows the
## linear-vs-log comparison that motivates using log throughout.
build_ladder_plot(
  LANGS,
  file.path(FIG_DIR, "glmer_ladder_mega.png"),
  title = "Full model grid: linear vs log age × 6 languages",
  width = 18, height = 16
)
