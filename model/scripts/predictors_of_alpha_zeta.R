## What predicts efficiency level (ξ_i ≈ log α_i + const) and change in
## efficiency (ζ_i)?
##
## Background:
##   Per the §28 GLMER equivalency, the random effects from
##   glmer(produces ~ 1 + log_age + (1 + log_age | child) + (1 | word),
##         family = binomial)
##   reproduce the Stan posterior medians of ξ_i (kid intercept = log r +
##   log α) and ζ_i (kid slope deviation = per-child κ_i deviation from
##   κ_pop). The two random effects are exactly the two latent dimensions
##   M_best uses to characterize a child.
##
##   Stan has no per-child σ_r decomposition, so we cannot separate
##   log r_i from log α_i WITHIN a fit. The data-identified
##   per-child quantity is ξ_i. Per-child α_i = (σ_α² / σ_ξ²) × (ξ_i − μ_r)
##   is just a re-scaled, re-centered ξ_i — perfect rank correlation, so
##   standardized regressions on demographics yield identical conclusions
##   whether we call the predictor "ξ_i" or "log α_i".
##
##   We use ξ_i (level) and ζ_i (change) directly.
##
## Wordbank demographics we ask:
##   - sex (Female / Male)
##   - birth_order (First / Later, dichotomized due to small higher cells)
##   - caregiver_education (Some College or less / College / Graduate)
##
## Data coverage (longitudinal bundles, I=500 each):
##   EN: sex 389/500, birth_order 112/500, caregiver_ed 138/500
##       (Smith and Thal subsets have no demographics)
##   NO: sex 500/500, birth_order 469/500, caregiver_ed 500/500
##       (single dataset: Kristoffersen — clean)
##
## Output:
##   outputs/predictors_alpha_zeta.csv  -- coefficient table
##   figures/longitudinal/predictors_alpha_zeta.png  -- coefficient plot

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(lme4); library(broom.mixed); library(wordbankr); library(knitr)
  library(forcats); library(patchwork)
})

OUT_CSV <- file.path(PATHS$outputs_dir, "predictors_alpha_zeta.csv")
OUT_PNG <- file.path(PATHS$figs_dir, "longitudinal",
                      "predictors_alpha_zeta.png")
dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)

## ---- Step 1: extract per-child BLUPs from each glmer fit -------
extract_blups <- function(fit_path, language) {
  m <- readRDS(fit_path)
  re <- ranef(m)$child
  tibble(
    child_id = as.integer(rownames(re)),
    xi       = re[, "(Intercept)"],
    zeta     = re[, "log_age"],
    language = language
  )
}

blups_en <- extract_blups(file.path(PATHS$fits_dir, "glmer_mbest_en.rds"),
                          "English")
no_path  <- file.path(PATHS$fits_dir, "glmer_mbest_no.rds")
if (file.exists(no_path)) {
  blups_no <- extract_blups(no_path, "Norwegian")
} else {
  warning("NO glmer not yet fit — skipping NO. Run glmer_mbest_norwegian.R first.")
  blups_no <- tibble()
}
cat(sprintf("EN BLUPs: %d kids; NO BLUPs: %d kids\n",
            nrow(blups_en), nrow(blups_no)))

## ---- Step 2: pull demographics -------------------------------
get_demo <- function(language_label) {
  forms <- c("WS", "WG")
  ad <- bind_rows(lapply(forms, function(f) {
    tryCatch(
      get_administration_data(language = language_label, form = f,
                              include_demographic_info = TRUE,
                              include_birth_info = TRUE),
      error = function(e) NULL
    )
  }))
  # take ANY non-NA demographic across this child's admins
  ad |>
    group_by(child_id) |>
    summarise(
      sex          = sex[!is.na(sex)][1],
      birth_order  = birth_order[!is.na(birth_order)][1],
      caregiver_ed = caregiver_education[!is.na(caregiver_education)][1],
      dataset      = dataset_name[1],
      .groups = "drop"
    )
}

cat("Pulling EN demographics...\n")
demo_en <- get_demo("English (American)")
cat("Pulling NO demographics...\n")
demo_no <- get_demo("Norwegian")

## ---- Step 3: join + recode -----------------------------------
recode_demo <- function(df) {
  df |>
    mutate(
      sex = factor(sex, levels = c("Female", "Male")),
      # birth order: First vs Later (combine 2+ due to thin cells)
      birth_order_2 = case_when(
        birth_order == "First" ~ "First",
        birth_order %in% c("Second", "Third", "Fourth", "Fifth",
                            "Sixth", "Seventh", "Eighth") ~ "Later",
        TRUE ~ NA_character_
      ),
      birth_order_2 = factor(birth_order_2, levels = c("First", "Later")),
      # caregiver_ed: 3-level (≤Some College / College / Graduate+)
      caregiver_ed_3 = case_when(
        caregiver_ed %in% c("None", "Primary", "Some Secondary",
                             "Secondary", "Some College") ~ "≤Some College",
        caregiver_ed == "College" ~ "College",
        caregiver_ed %in% c("Some Graduate", "Graduate") ~ "Graduate+",
        TRUE ~ NA_character_
      ),
      caregiver_ed_3 = factor(caregiver_ed_3,
                               levels = c("≤Some College",
                                           "College",
                                           "Graduate+"))
    )
}

en <- blups_en |> left_join(demo_en, by = "child_id") |> recode_demo()
if (nrow(blups_no) > 0) {
  no <- blups_no |> left_join(demo_no, by = "child_id") |> recode_demo()
} else {
  no <- tibble()
}

cat("\n=== EN coverage ===\n")
cat("  sex:", sum(!is.na(en$sex)), "/", nrow(en), "\n")
cat("  birth_order_2:", sum(!is.na(en$birth_order_2)), "/", nrow(en), "\n")
cat("  caregiver_ed_3:", sum(!is.na(en$caregiver_ed_3)), "/", nrow(en), "\n")
if (nrow(no) > 0) {
  cat("\n=== NO coverage ===\n")
  cat("  sex:", sum(!is.na(no$sex)), "/", nrow(no), "\n")
  cat("  birth_order_2:", sum(!is.na(no$birth_order_2)), "/", nrow(no), "\n")
  cat("  caregiver_ed_3:", sum(!is.na(no$caregiver_ed_3)), "/", nrow(no), "\n")
}

## ---- Step 4: regression analyses -----------------------------
## For each (language, outcome), fit one OLS with all three predictors.
## We standardize ξ and ζ to z-scores so coefficients are in SDs of the
## outcome — directly comparable across panels and languages.

fit_one <- function(data, outcome) {
  d <- data |>
    mutate(y = as.numeric(scale(.data[[outcome]]))) |>
    filter(!is.na(y), !is.na(sex), !is.na(birth_order_2),
            !is.na(caregiver_ed_3))
  if (nrow(d) < 20) return(NULL)
  m <- lm(y ~ sex + birth_order_2 + caregiver_ed_3, data = d)
  tidy(m, conf.int = TRUE) |>
    mutate(outcome = outcome, n = nrow(d))
}

combos <- expand.grid(
  language = c("English", "Norwegian"),
  outcome  = c("xi", "zeta"),
  stringsAsFactors = FALSE
)
results <- bind_rows(lapply(seq_len(nrow(combos)), function(i) {
  lang <- combos$language[i]
  out  <- combos$outcome[i]
  d <- if (lang == "English") en else no
  if (nrow(d) == 0) return(NULL)
  res <- fit_one(d, out)
  if (is.null(res)) return(NULL)
  res |> mutate(language = lang)
}))

cat("\n=== regression results (z-scored outcomes) ===\n")
print(results |>
        select(language, outcome, term, estimate, std.error,
                conf.low, conf.high, p.value, n) |>
        mutate(across(c(estimate, std.error, conf.low, conf.high, p.value),
                       \(x) round(x, 3))))

write_csv(results, OUT_CSV)
cat(sprintf("\nWrote %s\n", OUT_CSV))

## ---- Step 5: also fit MARGINAL (one-predictor) models -------
## EN coverage is limited; if we condition on having ALL three
## predictors observed, sample shrinks to ~80. So also report
## marginal effects (one predictor at a time) to use the full
## available sample for each variable.
fit_marginal <- function(data, outcome, predictor) {
  d <- data |>
    mutate(y = as.numeric(scale(.data[[outcome]]))) |>
    filter(!is.na(y), !is.na(.data[[predictor]]))
  if (nrow(d) < 20) return(NULL)
  m <- lm(as.formula(paste("y ~", predictor)), data = d)
  tidy(m, conf.int = TRUE) |>
    mutate(outcome = outcome, predictor = predictor, n = nrow(d))
}

predictors <- c("sex", "birth_order_2", "caregiver_ed_3")
marginals <- bind_rows(lapply(seq_len(nrow(combos)), function(i) {
  lang <- combos$language[i]
  out  <- combos$outcome[i]
  d <- if (lang == "English") en else no
  if (nrow(d) == 0) return(NULL)
  bind_rows(lapply(predictors, function(p) {
    res <- fit_marginal(d, out, p)
    if (is.null(res)) return(NULL)
    res |> mutate(language = lang)
  }))
}))

cat("\n=== marginal (single-predictor) z-scored coefficients ===\n")
print(marginals |>
        filter(term != "(Intercept)") |>
        select(language, outcome, predictor, term, estimate,
                conf.low, conf.high, p.value, n) |>
        mutate(across(c(estimate, conf.low, conf.high, p.value),
                       \(x) round(x, 3))))

write_csv(marginals,
          file.path(PATHS$outputs_dir, "predictors_alpha_zeta_marginal.csv"))

## ---- Step 6: coefficient plot --------------------------------
## Plot marginal coefficients (full-sample) since EN multivariate
## sample is tiny.
plot_data <- marginals |>
  filter(term != "(Intercept)") |>
  mutate(
    term_label = recode(term,
      sexMale                   = "Male\n(vs Female)",
      birth_order_2Later        = "Later-born\n(vs First)",
      `caregiver_ed_3College`   = "College\n(vs ≤Some College)",
      `caregiver_ed_3Graduate+` = "Graduate+\n(vs ≤Some College)"
    ),
    outcome_label = recode(outcome,
      xi   = "Efficiency level (ξ)",
      zeta = "Change in efficiency (ζ)"
    ),
    sig = ifelse(p.value < 0.05, "p<.05", "ns")
  )

p <- ggplot(plot_data,
             aes(x = estimate, y = term_label, colour = language,
                 shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                  height = 0, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  facet_wrap(~ outcome_label, ncol = 2) +
  scale_colour_manual(values = c(English = "#d7191c",
                                  Norwegian = "#1f78b4")) +
  scale_shape_manual(values = c(`p<.05` = 16, ns = 1)) +
  labs(x = "Standardized coefficient (SDs of outcome)",
       y = NULL,
       title = "What predicts efficiency (ξ) vs change in efficiency (ζ)?",
       subtitle = "Marginal OLS, z-scored outcome. EN: per-predictor n in {138-389}; NO: per-predictor n in {469-500}.",
       colour = NULL, shape = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        legend.position = "top")

ggsave(OUT_PNG, p, width = 10, height = 5, dpi = 200)
cat(sprintf("Wrote %s\n", OUT_PNG))

## ---- Step 7: 2x3 panel plot ----------------------------------
## Rows: ξ (efficiency level), ζ (change in efficiency)
## Cols: Sex, Birth order, Caregiver education
## Both languages overlaid: jittered points + per-language category
## means with error bars + connecting line.

both <- bind_rows(
  en |> select(child_id, language, xi, zeta, sex,
                birth_order_2, caregiver_ed_3),
  if (nrow(no) > 0) no |>
    select(child_id, language, xi, zeta, sex,
            birth_order_2, caregiver_ed_3) else NULL
)

## tidy-long: one row per (kid, outcome, predictor) for facet_grid
panel_long <- bind_rows(
  both |> filter(!is.na(sex)) |>
    transmute(language, outcome = "xi",   x = sex,
              x_num = as.numeric(sex), y = xi,
              predictor = "Sex"),
  both |> filter(!is.na(birth_order_2)) |>
    transmute(language, outcome = "xi",   x = birth_order_2,
              x_num = as.numeric(birth_order_2), y = xi,
              predictor = "Birth order"),
  both |> filter(!is.na(caregiver_ed_3)) |>
    transmute(language, outcome = "xi",   x = caregiver_ed_3,
              x_num = as.numeric(caregiver_ed_3), y = xi,
              predictor = "Caregiver education"),
  both |> filter(!is.na(sex)) |>
    transmute(language, outcome = "zeta", x = sex,
              x_num = as.numeric(sex), y = zeta,
              predictor = "Sex"),
  both |> filter(!is.na(birth_order_2)) |>
    transmute(language, outcome = "zeta", x = birth_order_2,
              x_num = as.numeric(birth_order_2), y = zeta,
              predictor = "Birth order"),
  both |> filter(!is.na(caregiver_ed_3)) |>
    transmute(language, outcome = "zeta", x = caregiver_ed_3,
              x_num = as.numeric(caregiver_ed_3), y = zeta,
              predictor = "Caregiver education")
) |>
  mutate(
    outcome   = factor(outcome,
                        levels = c("xi", "zeta"),
                        labels = c("Efficiency level (ξ)",
                                    "Change in efficiency (ζ)")),
    predictor = factor(predictor,
                        levels = c("Sex", "Birth order",
                                    "Caregiver education"))
  )

## per-language category means + SE for the overlay points/lines
panel_means <- panel_long |>
  group_by(language, outcome, predictor, x, x_num) |>
  summarise(mean_y = mean(y),
            se     = sd(y) / sqrt(n()),
            n      = n(),
            .groups = "drop")

p_panel <- ggplot(panel_long, aes(x = x_num, y = y, colour = language)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.15, size = 0.6) +
  geom_line(data = panel_means,
             aes(x = x_num, y = mean_y, group = language),
             linewidth = 0.9) +
  geom_errorbar(data = panel_means,
                 aes(x = x_num, y = mean_y,
                     ymin = mean_y - se, ymax = mean_y + se),
                 width = 0.12, linewidth = 0.7) +
  geom_point(data = panel_means,
              aes(x = x_num, y = mean_y),
              size = 3) +
  facet_grid(outcome ~ predictor, scales = "free_x") +
  scale_colour_manual(values = c(English   = "#d7191c",
                                  Norwegian = "#1f78b4")) +
  scale_x_continuous(
    breaks = function(limits) {
      # tick at every integer in range
      seq(ceiling(limits[1]), floor(limits[2]))
    },
    labels = function(breaks) {
      # need to translate numeric breaks back to factor labels for the
      # CURRENT panel. We can't easily get the facet here, so we rely
      # on a custom labeller below.
      breaks
    }
  ) +
  labs(x = NULL,
       y = "Per-child BLUP (mean-centered logit units)",
       title = "What predicts efficiency (ξ) vs change in efficiency (ζ)?",
       subtitle = "Per-child glmer random effects (= Stan posterior medians, §28). Jittered points = kids; dark points + line = per-language category mean ± 1 SE.",
       colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        strip.text    = element_text(face = "bold"),
        legend.position = "top",
        panel.spacing = unit(1.0, "lines"))

## Custom labels per facet column. We use a labeller-like trick:
## the simplest robust approach is to write x labels directly via
## a per-panel scale, but ggplot facet_grid uses a shared scale.
## Instead, set x as a factor with all 5 levels and let ggplot drop
## unused ones. Cleaner:

panel_long2 <- panel_long |>
  mutate(
    x_label = case_when(
      predictor == "Sex" ~ as.character(x),
      predictor == "Birth order" ~ as.character(x),
      predictor == "Caregiver education" ~ as.character(x)
    ),
    x_label = factor(x_label,
                      levels = c("Female", "Male",
                                  "First", "Later",
                                  "≤Some College", "College", "Graduate+"))
  )
panel_means2 <- panel_means |>
  mutate(
    x_label = factor(as.character(x),
                      levels = c("Female", "Male",
                                  "First", "Later",
                                  "≤Some College", "College", "Graduate+"))
  )

## Build annotation table: per (outcome, predictor, language), report
## the slope of the "trend" (= effect of the last vs first category in
## standardized terms). Pull from `marginals`.
ann <- marginals |>
  filter(term != "(Intercept)") |>
  mutate(
    predictor_lab = recode(predictor,
                            sex            = "Sex",
                            birth_order_2  = "Birth order",
                            caregiver_ed_3 = "Caregiver education"),
    outcome_lab   = recode(outcome,
                            xi   = "Efficiency level (ξ)",
                            zeta = "Change in efficiency (ζ)"),
    sig = case_when(p.value < 0.001 ~ "***",
                    p.value < 0.01  ~ "**",
                    p.value < 0.05  ~ "*",
                    TRUE            ~ "ns"),
    label = sprintf("β=%+.2f %s", estimate, sig)
  ) |>
  # for caregiver_ed, pick only the largest of College/Graduate+ effect
  # to keep annotations compact (one per language per panel)
  group_by(language, outcome, predictor) |>
  slice_max(abs(estimate), n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    outcome   = factor(outcome,
                        levels = c("xi", "zeta"),
                        labels = c("Efficiency level (ξ)",
                                    "Change in efficiency (ζ)")),
    predictor = factor(predictor_lab,
                        levels = c("Sex", "Birth order",
                                    "Caregiver education"))
  )

# y-position per facet row: top of the ξ panel = max ξ, bottom of ζ = min ζ.
# Easier: use facet-specific y-positions via in-aes data; we let ggplot
# auto-place by language with a vertical offset.
ann_pos <- panel_long2 |>
  group_by(outcome, predictor) |>
  summarise(y_top = quantile(y, 0.98),
            y_bot = quantile(y, 0.02), .groups = "drop")
ann <- ann |>
  left_join(ann_pos, by = c("outcome", "predictor")) |>
  mutate(y_ann = ifelse(language == "English", y_top, y_top * 0.78))

p_panel <- ggplot(panel_long2,
                   aes(x = x_label, y = y, colour = language)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.12, size = 0.55) +
  geom_line(data = panel_means2,
             aes(x = x_label, y = mean_y, group = language),
             linewidth = 0.9) +
  geom_errorbar(data = panel_means2,
                 aes(x = x_label, y = mean_y,
                     ymin = mean_y - se, ymax = mean_y + se),
                 width = 0.12, linewidth = 0.7) +
  geom_point(data = panel_means2,
              aes(x = x_label, y = mean_y),
              size = 3) +
  geom_text(data = ann,
             aes(x = 1, y = y_ann, label = label, colour = language),
             hjust = 0, size = 3.2, fontface = "bold",
             show.legend = FALSE) +
  facet_grid(outcome ~ predictor,
              scales = "free", space = "free_x") +
  scale_colour_manual(values = c(English   = "#d7191c",
                                  Norwegian = "#1f78b4")) +
  labs(x = NULL,
       y = "Per-child BLUP (mean-centered logit units)",
       title = "What predicts efficiency (ξ) vs change in efficiency (ζ)?",
       subtitle = "Per-child glmer random effects (= Stan posterior medians, §28). Points = kids; lines = per-language category mean ± 1 SE. β = standardized marginal coefficient (SDs of outcome). *p<.05, **p<.01, ***p<.001.",
       colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey25"),
        strip.text    = element_text(face = "bold"),
        legend.position = "top",
        axis.text.x   = element_text(angle = 0, hjust = 0.5, size = 9),
        panel.spacing = unit(1.0, "lines"))

panel_png <- file.path(PATHS$figs_dir, "longitudinal",
                        "predictors_alpha_zeta_panel.png")
ggsave(panel_png, p_panel, width = 11, height = 6.5, dpi = 200)
cat(sprintf("Wrote %s\n", panel_png))

cat("\nDone.\n")
