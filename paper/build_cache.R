## build_cache.R — populates outputs/paper/cache/ with the small RDS
## summaries that standard_model.qmd loads in its setup chunk.
##
## Run from project root:
##   Rscript outputs/paper/build_cache.R
##
## Re-run when:
##   - underlying fits change
##   - dataset survey changes
##   - we add languages or studies

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(here)
  library(lme4); library(broom.mixed); library(posterior); library(wordbankr)
})

PAPER_LANGS <- c("english_american", "norwegian",
                 "french_quebecois", "japanese")
LANG_LABELS <- c("english_american" = "English (American)",
                 "norwegian"        = "Norwegian",
                 "french_quebecois" = "French (Quebecois)",
                 "japanese"         = "Japanese")
WORDBANKR_LABELS <- c("english_american" = "English (American)",
                      "norwegian"        = "Norwegian",
                      "french_quebecois" = "French (Quebecois)",
                      "japanese"         = "Japanese")

CACHE <- here("paper", "cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

## ---- 1. dataset table values from the language survey -----------
survey <- read_csv(here("outputs", "glmer_ladder", "00_language_survey.csv"),
                   show_col_types = FALSE)

wb_rows <- survey |>
  filter(language %in% LANG_LABELS) |>
  mutate(
    panel       = "Wordbank longitudinal",
    dataset     = language,
    n_kids      = long_kids,
    admins_per_kid_mean = round(admins / total_kids, 1),
    age_range   = paste0(age_min, "–", age_max, " mo"),
    input_info  = "—",
    citation    = "[@CITE]"
  ) |>
  select(panel, dataset, citation, n_kids,
         admins_per_kid = admins_per_kid_mean, age_range, input_info)

io_rows <- tibble::tribble(
  ~panel,                ~dataset,    ~citation, ~n_kids, ~admins_per_kid, ~age_range,  ~input_info,
  "Input-observed (IO)", "BabyView",  "[@CITE]", 22L,     5.0,             "6–32 mo",   "Head-cam (n = 5,739 videos)",
  "Input-observed (IO)", "SEEDLingS", "[@CITE]", 44L,     12.0,            "6–17 mo",   "LENA day (n = 525 recordings)",
  "Input-observed (IO)", "AM2018",    "[@CITE]", 66L,     3.0,             "16–24 mo",  "LENA (n = 126 recordings)",
  "Input-observed (IO)", "FMW2013",   "[@CITE]", 51L,     3.0,             "18–30 mo",  "LENA (n = 51 recordings)"
)

table1 <- bind_rows(wb_rows, io_rows)
saveRDS(table1, file.path(CACHE, "table1_datasets.rds"))
cat(sprintf("Wrote %s (%d rows)\n",
            file.path(CACHE, "table1_datasets.rds"), nrow(table1)))

## ---- 2. glmer-ladder predictions + empirical points -------------
sc <- readRDS(here("outputs", "glmer_ladder", "sim_cache.rds"))

# keep only the four paper languages and the log-age models
qt <- sc$qtiles |>
  filter(lang_slug %in% PAPER_LANGS,
         model %in% c("A", "B_log", "C_log", "D_log"))
ep <- sc$emp |>
  filter(lang_slug %in% PAPER_LANGS)

saveRDS(list(qtiles = qt, emp = ep, langs = PAPER_LANGS,
             lang_labels = LANG_LABELS[PAPER_LANGS]),
        file.path(CACHE, "fig2_glmer_ladder.rds"))
cat(sprintf("Wrote %s (qtiles=%d, emp=%d rows)\n",
            file.path(CACHE, "fig2_glmer_ladder.rds"),
            nrow(qt), nrow(ep)))

## ---- 3. BLUPs from D_log fits + Wordbank demographics ----------
extract_blups <- function(lang_slug) {
  fp <- here("fits", "glmer_ladder",
             sprintf("fit_%s_D_log.rds", lang_slug))
  if (!file.exists(fp)) return(NULL)
  m  <- readRDS(fp)
  re <- ranef(m)$child
  if (is.null(re)) re <- ranef(m)[[grep("child|kid", names(ranef(m)),
                                         ignore.case = TRUE)[1]]]
  slope_col <- intersect(c("log_age", "la"), colnames(re))[1]
  tibble::tibble(
    child_id = as.integer(rownames(re)),
    xi       = re[, "(Intercept)"],
    zeta     = re[, slope_col],
    lang_slug = lang_slug
  )
}

cat("Extracting BLUPs from D_log fits ...\n")
blups <- bind_rows(lapply(PAPER_LANGS, function(l) {
  b <- extract_blups(l); if (is.null(b)) cat(" - skipping", l, "\n")
  b
}))
cat(sprintf("  total BLUPs: %d kids across %d languages\n",
            nrow(blups), length(unique(blups$lang_slug))))

cat("Pulling Wordbank demographics ...\n")
get_demo <- function(lang_label) {
  forms <- c("WS", "WG")
  ad <- bind_rows(lapply(forms, function(f) {
    tryCatch(
      get_administration_data(language = lang_label, form = f,
                              include_demographic_info = TRUE,
                              include_birth_info = TRUE),
      error = function(e) NULL)
  }))
  if (is.null(ad) || nrow(ad) == 0) return(NULL)
  ad |>
    group_by(child_id) |>
    summarise(
      sex          = sex[!is.na(sex)][1],
      birth_order  = birth_order[!is.na(birth_order)][1],
      caregiver_ed = caregiver_education[!is.na(caregiver_education)][1],
      .groups = "drop"
    )
}

demos <- bind_rows(lapply(PAPER_LANGS, function(l) {
  d <- get_demo(WORDBANKR_LABELS[l])
  if (is.null(d)) NULL else mutate(d, lang_slug = l)
}))
cat(sprintf("  demographics: %d rows\n", nrow(demos)))

# join + recode
recode_demo <- function(df) {
  df |>
    mutate(
      sex = factor(sex, levels = c("Female", "Male")),
      birth_order_2 = case_when(
        birth_order == "First" ~ "First",
        birth_order %in% c("Second", "Third", "Fourth", "Fifth",
                            "Sixth", "Seventh", "Eighth") ~ "Later",
        TRUE ~ NA_character_
      ) |> factor(levels = c("First", "Later")),
      caregiver_ed_3 = case_when(
        caregiver_ed %in% c("None", "Primary", "Some Secondary",
                             "Secondary", "Some College") ~ "<=Some College",
        caregiver_ed == "College" ~ "College",
        caregiver_ed %in% c("Some Graduate", "Graduate") ~ "Graduate+",
        TRUE ~ NA_character_
      ) |> factor(levels = c("<=Some College", "College", "Graduate+"))
    )
}
joined <- blups |>
  left_join(demos, by = c("lang_slug", "child_id")) |>
  recode_demo()
saveRDS(joined, file.path(CACHE, "blups_demographics.rds"))
cat(sprintf("Wrote %s (%d rows)\n",
            file.path(CACHE, "blups_demographics.rds"), nrow(joined)))

## ---- 4. IO summary (wide-delta) for figure 5 -------------------
b <- readRDS(here("fits", "io_pooled_subset_data.rds"))
fit <- readRDS(here("fits", "io_pooled_widedelta.rds"))
fitg <- readRDS(here("fits", "io_pooled_gamma_widedelta_add.rds"))
d <- fit$draws(format = "df"); dg <- fitg$draws(format = "df")
io_summary <- list(
  obs = b$df |>
    group_by(study, ii, aa, age) |>
    summarise(prop = mean(produces), .groups = "drop"),
  params = list(
    sigma_r            = median(d$sigma_r),
    sigma_alpha        = median(d$sigma_alpha),
    sigma_zeta_base    = median(d$sigma_zeta),
    sigma_zeta_add     = median(dg$sigma_zeta),
    gamma_add          = median(dg$gamma),
    kappa_pop          = median(d$delta) + 1
  )
)
saveRDS(io_summary, file.path(CACHE, "fig5_io_summary.rds"))
cat(sprintf("Wrote %s\n", file.path(CACHE, "fig5_io_summary.rds")))

cat("\nAll caches built.\n")
