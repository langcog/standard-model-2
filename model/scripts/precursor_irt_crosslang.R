## Cross-language precursor analysis: model-free 1-PL Rasch on Wordbank
## WS data, replicating the English precursor across as many languages
## as Wordbank supports.
##
## For each language with a productive (WS) form and >= 200 admins:
##   1. Pull cross-sectional admin × item produces matrix from wordbankr
##   2. Filter to ages 12-30 (covers the productive-vocab window)
##   3. Fit y_aj ~ 1 + (1|aa) + (1|jj) via glmer  [1-PL Rasch IRT]
##   4. Extract per-admin theta, compute:
##        - pop age slope (logits/mo)
##        - within-age SD of theta (residual after age trend)
##        - SD vs. age (per integer age bin)
##   5. Save per-language scalars + per-age-bin SD as CSVs
##
## Outputs:
##   outputs/figs/longitudinal/precursor_crosslang_summary.csv
##   outputs/figs/longitudinal/precursor_crosslang_sd_by_age.csv
##   outputs/figs/longitudinal/precursor_crosslang_panel.png
##   outputs/figs/longitudinal/precursor_crosslang_summary.png
##
## Run:  Rscript model/scripts/precursor_irt_crosslang.R
##
## Per-language fits are cached at fits/precursor/<lang_slug>.rds
## so re-running picks up where it left off.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(lme4); library(wordbankr)
})

OUT_DIR <- file.path(PATHS$figs_dir, "longitudinal")
CACHE_DIR <- file.path(PATHS$fits_dir, "precursor")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Discover languages with WS form ---------------------------------
instruments <- get_instruments()
ws_languages <- instruments |>
  filter(form == "WS") |>
  pull(language) |>
  unique() |>
  sort()
cat("Languages with WS:", length(ws_languages), "—",
    paste(head(ws_languages, 10), collapse = ", "), "...\n")

# Performance ceiling: subsample admins per language so glmer finishes
# in tractable time. The downstream estimand is the SD-by-age CURVE,
# which is well-preserved at moderate sample sizes; we don't need every
# admin. With MAX_ADMINS=500 stratified by age, glmer fits in 1-3 min
# per language vs. tens of minutes at 1500.
MAX_ADMINS <- 500    # per language
AGE_RANGE  <- c(12, 30)
MIN_OBS_PER_AGE <- 5  # for SD-by-age computation

# ---- Per-language fit driver -----------------------------------------
fit_one_language <- function(lang) {
  slug <- gsub("[^A-Za-z0-9]+", "_", lang)
  cache_file <- file.path(CACHE_DIR, sprintf("%s.rds", slug))
  if (file.exists(cache_file)) {
    cat(sprintf("[%s] cached, loading\n", lang))
    return(readRDS(cache_file))
  }
  cat(sprintf("\n[%s] pulling data ...\n", lang))

  # Single-shot pull: get_instrument_data returns one row per (admin x item)
  # for the entire WS form in this language, with admin metadata attached.
  raw <- tryCatch(
    get_instrument_data(language = lang, form = "WS",
                        administration_info = TRUE, item_info = TRUE),
    error = function(e) {
      cat(sprintf("[%s] pull failed: %s\n", lang, conditionMessage(e)))
      NULL
    })
  if (is.null(raw) || nrow(raw) == 0) return(NULL)

  # Filter: words only, age in range. Wordbank's `produces` column is
  # already binary 0/1 (or NA); use it directly.
  if (!"produces" %in% names(raw)) {
    raw$produces <- as.integer(raw$value == "produces")
  }
  df <- raw |>
    filter(item_kind == "word",
           age >= AGE_RANGE[1], age <= AGE_RANGE[2]) |>
    select(child_id, age, item_definition, produces) |>
    filter(!is.na(produces))
  if (nrow(df) == 0) return(NULL)

  n_admins_all <- df |> distinct(child_id, age) |> nrow()
  if (n_admins_all < 200) {
    cat(sprintf("[%s] only %d (child x age) admins in age range; skipping\n",
                lang, n_admins_all))
    return(NULL)
  }

  # Subsample admins per age bin to keep glmer tractable
  admin_keys <- df |> distinct(child_id, age) |>
    mutate(age_int = round(age))
  if (nrow(admin_keys) > MAX_ADMINS) {
    set.seed(2026)
    per_bin <- ceiling(MAX_ADMINS / length(unique(admin_keys$age_int)))
    admin_keys <- admin_keys |>
      group_by(age_int) |>
      group_modify(~ slice_sample(.x, n = min(nrow(.x), per_bin))) |>
      ungroup()
    cat(sprintf("[%s] subsampled %d -> %d admins (per_bin=%d)\n",
                lang, n_admins_all, nrow(admin_keys), per_bin))
    df <- df |> semi_join(admin_keys, by = c("child_id", "age"))
  }

  # Re-index aa (admin) and jj (item)
  admin_idx <- df |> distinct(child_id, age) |>
    arrange(child_id, age) |> mutate(aa = row_number())
  item_idx  <- df |> distinct(item_definition) |>
    arrange(item_definition) |> mutate(jj = row_number())
  df <- df |>
    left_join(admin_idx, by = c("child_id", "age")) |>
    left_join(item_idx, by = "item_definition")
  admin_meta <- admin_idx

  cat(sprintf("[%s] fitting glmer on %d obs (%d admins x %d items) ...\n",
              lang, nrow(df), max(df$aa), max(df$jj)))
  t0 <- Sys.time()
  fit <- tryCatch(
    glmer(produces ~ 1 + (1 | aa) + (1 | jj),
          data = df, family = binomial,
          control = glmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 2e5))),
    error = function(e) {
      cat(sprintf("[%s] glmer failed: %s\n", lang, conditionMessage(e)))
      NULL
    })
  if (is.null(fit)) return(NULL)
  cat(sprintf("[%s] glmer fit time: %.1fs\n",
              lang, as.numeric(Sys.time() - t0, units = "secs")))

  # Extract per-admin theta = global intercept + ranef(aa)
  re <- ranef(fit)$aa
  theta_a <- tibble::tibble(
    aa = as.integer(rownames(re)),
    theta = re[, "(Intercept)"] + fixef(fit)["(Intercept)"]
  ) |>
    inner_join(admin_meta |> select(aa, age, child_id), by = "aa")

  # Population age slope + within-age SD
  slope_fit <- lm(theta ~ age, data = theta_a)
  pop_slope <- as.numeric(coef(slope_fit)["age"])
  sigma_theta_pooled <- sd(theta_a$theta)
  sigma_theta_within <- sd(residuals(slope_fit))

  age_summary <- theta_a |>
    mutate(age_int = round(age)) |>
    group_by(age_int) |>
    summarise(n = n(), mean_theta = mean(theta),
              sd_theta = sd(theta), .groups = "drop") |>
    filter(n >= MIN_OBS_PER_AGE)

  out <- list(
    language = lang,
    n_admins = nrow(theta_a),
    n_items  = max(df$jj),
    n_obs    = nrow(df),
    age_range = c(min(theta_a$age), max(theta_a$age)),
    pop_slope = pop_slope,
    sigma_theta_pooled = sigma_theta_pooled,
    sigma_theta_within = sigma_theta_within,
    months_per_sd_within = sigma_theta_within / pop_slope,
    age_summary = age_summary
  )
  saveRDS(out, cache_file)
  cat(sprintf("[%s] DONE: slope=%.3f, sigma_within=%.2f, months/SD=%.1f\n",
              lang, pop_slope, sigma_theta_within,
              sigma_theta_within / pop_slope))
  out
}

# ---- Run all languages -----------------------------------------------
results <- list()
for (lang in ws_languages) {
  r <- tryCatch(fit_one_language(lang),
                error = function(e) {
                  cat(sprintf("[%s] FAILED: %s\n", lang, conditionMessage(e)))
                  NULL
                })
  if (!is.null(r)) results[[lang]] <- r
}

cat(sprintf("\nFitted %d / %d languages.\n",
            length(results), length(ws_languages)))

# ---- Aggregate scalars table ----------------------------------------
summary_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(language = r$language,
             n_admins = r$n_admins, n_items = r$n_items,
             age_min = r$age_range[1], age_max = r$age_range[2],
             pop_slope = r$pop_slope,
             sigma_theta_pooled = r$sigma_theta_pooled,
             sigma_theta_within = r$sigma_theta_within,
             months_per_sd_within = r$months_per_sd_within,
             stringsAsFactors = FALSE)
})) |> arrange(desc(n_admins))
write.csv(summary_df,
          file.path(OUT_DIR, "precursor_crosslang_summary.csv"),
          row.names = FALSE)
cat(sprintf("Wrote: %s\n",
            file.path(OUT_DIR, "precursor_crosslang_summary.csv")))
cat("\n=== Per-language scalars ===\n")
print(summary_df, digits = 3, row.names = FALSE)

# ---- Aggregate per-age-bin SD table ---------------------------------
sd_by_age <- do.call(rbind, lapply(results, function(r) {
  data.frame(language = r$language, age_int = r$age_summary$age_int,
             n = r$age_summary$n,
             mean_theta = r$age_summary$mean_theta,
             sd_theta = r$age_summary$sd_theta,
             stringsAsFactors = FALSE)
}))
write.csv(sd_by_age,
          file.path(OUT_DIR, "precursor_crosslang_sd_by_age.csv"),
          row.names = FALSE)
cat(sprintf("Wrote: %s\n",
            file.path(OUT_DIR, "precursor_crosslang_sd_by_age.csv")))

# ---- Plots ----------------------------------------------------------
plot_df <- sd_by_age |>
  group_by(language) |>
  filter(n() >= 4) |>            # at least 4 age bins to show a curve
  ungroup() |>
  mutate(language = factor(language,
                           levels = summary_df$language))

p_panel <- ggplot(plot_df, aes(age_int, sd_theta, group = language)) +
  geom_line(linewidth = 0.5, colour = "grey25", alpha = 0.7) +
  geom_point(size = 0.9, colour = "grey25") +
  facet_wrap(~ language, ncol = 5) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = "Age (months)",
       y = expression("Within-age SD of " * theta[admin] * " (logits)"),
       title = "Cross-language replication: within-age SD of ability grows with age",
       subtitle = "1-PL Rasch IRT on Wordbank cross-sectional WS data; per-admin theta extracted; SD computed per integer age bin (n>=5).") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "precursor_crosslang_panel.png"),
       p_panel, width = 13, height = 8, dpi = 150)
cat(sprintf("Wrote: %s\n",
            file.path(OUT_DIR, "precursor_crosslang_panel.png")))

# Summary scatter: pop slope vs months-per-SD by language
p_summary <- ggplot(summary_df,
                    aes(pop_slope, months_per_sd_within, label = language)) +
  geom_point(aes(size = n_admins), alpha = 0.7, colour = "steelblue") +
  geom_text(size = 2.8, hjust = -0.15, vjust = 0.4) +
  scale_size_continuous(range = c(2, 8), name = "N admins") +
  labs(x = "Population age slope (logits / month)",
       y = "Within-age SD of theta / pop slope (months/SD)",
       title = "Variability in 'months of development' equivalents, by language",
       subtitle = "How many months of typical development +/- 1 SD across kids corresponds to.") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "precursor_crosslang_summary.png"),
       p_summary, width = 10, height = 6.5, dpi = 150)
cat(sprintf("Wrote: %s\n",
            file.path(OUT_DIR, "precursor_crosslang_summary.png")))

cat("\nDone.\n")
