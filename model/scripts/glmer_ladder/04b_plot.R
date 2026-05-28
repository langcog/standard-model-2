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

## ---- Build the plot for a given language subset ----------------------
build_ladder_plot <- function(langs_subset, out_png, title_n_langs,
                                width = 18, height = 16) {
  qt <- qtiles |>
    filter(lang_slug %in% langs_subset) |>
    mutate(language = factor(LANG_LABELS[lang_slug],
                              levels = LANG_LABELS[langs_subset]),
            model    = factor(model, levels = MODELS))
  el <- emp |>
    filter(lang_slug %in% langs_subset) |>
    mutate(language = factor(LANG_LABELS[lang_slug],
                              levels = LANG_LABELS[langs_subset]))
  al <- summ |>
    filter(lang_slug %in% langs_subset) |>
    group_by(language, lang_slug) |>
    mutate(AIC_best = min(AIC), dAIC = AIC - AIC_best) |>
    ungroup() |>
    mutate(language = factor(language,
                              levels = LANG_LABELS[langs_subset]),
            model    = factor(model, levels = MODELS),
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
    facet_grid(language ~ model, scales = "free_y", space = "fixed") +
    geom_text(data = al,
               aes(x = -Inf, y = Inf, label = label),
               hjust = -0.15, vjust = 1.5, size = 2.6, colour = "grey25") +
    scale_colour_manual(values = WORDBANK_PALETTE, name = "Percentile") +
    scale_linewidth_identity() +
    labs(x = "Age (months)",
         y = "Productive vocabulary",
         title = sprintf("glmer ladder: predictions vs empirical vocab(age) — %d languages × 7 models",
                          title_n_langs),
         subtitle = sprintf("Lines = 10/25/50/75/90 quantiles across %d simulated kids (BLUP bootstrap from each model's child-RE distribution). Grey lines = per-kid longitudinal trajectories. Empirical and predictions both restricted to the largest form per language. Corner label = ΔAIC vs best within language.",
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

## Main-text version: 4 well-powered languages
LANGS_MAIN <- c("english_american", "norwegian",
                "french_quebecois", "japanese")
build_ladder_plot(
  LANGS_MAIN,
  file.path(FIG_DIR, "glmer_ladder_main.png"),
  title_n_langs = length(LANGS_MAIN),
  width = 18, height = 11
)

## Supplementary version: all 6 languages
build_ladder_plot(
  LANGS,
  file.path(FIG_DIR, "glmer_ladder_mega.png"),
  title_n_langs = length(LANGS),
  width = 18, height = 16
)
