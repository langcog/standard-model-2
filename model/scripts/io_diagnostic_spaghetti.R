## Diagnostic spaghetti for the anchored IO fits.
##
## Motivation: AM2018 came back with kappa_pop ~ 1.9 (vs ~8 for BabyView
## / SEEDLingS). AM2018 is the only IO dataset with a WG->WS form
## transition inside the fitted ages (WG 13-18, WS 20-27). Question: do
## the same kids' WG and WS admins line up, or is there a level jump at
## the form boundary?
##
## Two rows:
##   raw count   — vocab = sum(produces); WG and WS have DIFFERENT item
##                 denominators, so a jump here can be pure denominator.
##   proportion  — vocab / n_items_on_form; comparable across forms, so a
##                 jump HERE is a real form-effect problem.
## Cols = BabyView | SEEDLingS | AM2018. Per-kid spaghetti, coloured by
## form, with the per-form item count annotated.
##
## Output: outputs/figs/io/io_diagnostic_spaghetti.png

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT_PNG <- file.path(PATHS$figs_dir, "io", "io_diagnostic_spaghetti.png")

DATASETS <- list(
  BabyView  = "babyview_subset_data.rds",
  SEEDLingS = "seedlings_subset_data.rds",
  AM2018    = "io_am2018_subset_data.rds",
  FMW2013   = "io_fmw2013_subset_data.rds"
)

emp_for <- function(label, path) {
  b <- readRDS(file.path(PATHS$fits_dir, path))
  d <- b$df
  idcol <- if ("child_id" %in% names(d)) "child_id" else "subject_id"
  d$ID <- as.character(d[[idcol]])
  # per (kid, age, form): vocab + n_items scored
  per_admin <- d |>
    distinct(ID, age, form, item, .keep_all = TRUE) |>
    group_by(ID, age, form) |>
    summarise(vocab = sum(produces, na.rm = TRUE),
              n_items = n(), .groups = "drop") |>
    mutate(prop = vocab / n_items, dataset = label)
  per_admin
}

emp <- bind_rows(Map(emp_for, names(DATASETS), DATASETS))
emp$dataset <- factor(emp$dataset, levels = names(DATASETS))

# per-form item counts for annotation
form_n <- emp |> group_by(dataset, form) |>
  summarise(n_items = max(n_items), .groups = "drop")
cat("Per-form item counts:\n"); print(form_n)

FORM_COL <- c(WG = "#1f78b4", WS = "#e31a1c", WGProd = "#1f78b4",
              WGProdShort = "#a6cee3", WSShort = "#fb9a99")

mk <- function(yvar, ylab) {
  ggplot(emp, aes(x = age, y = .data[[yvar]], group = ID)) +
    geom_line(alpha = 0.25, linewidth = 0.3, colour = "grey40") +
    geom_point(aes(colour = form), alpha = 0.7, size = 1.1) +
    facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = FORM_COL, name = "Form") +
    labs(x = "Age (months)", y = ylab) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top", strip.text = element_text(face = "bold"))
}

p_count <- mk("vocab", "Vocabulary count (raw)") +
  labs(title = "IO diagnostic: per-kid spaghetti by form",
       subtitle = "Lines connect a child's admins across age. Raw count: WG/WS have different item denominators (see proportion row).")
p_prop <- mk("prop", "Proportion of form's items produced")

p <- p_count / p_prop
ggsave(OUT_PNG, p, width = 16, height = 8, dpi = 200)
cat(sprintf("Wrote %s\n", OUT_PNG))
