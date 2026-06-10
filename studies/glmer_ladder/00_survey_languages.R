## Survey wordbank for languages with enough longitudinal data to run
## the model-ladder comparison. We want, per language:
##   * total kids
##   * kids with ≥ 2 admins (any form)
##   * total production observations (across forms)
##   * age range
##
## Threshold: ≥ 100 longitudinal kids (Mike's spec). We combine WS and
## WG production data at the item level — same word ("dog") gets the
## same item RE regardless of which form it came in on.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(wordbankr); library(dplyr); library(tidyr); library(readr)
})

OUT <- file.path(PATHS$outputs_dir, "glmer_ladder/00_language_survey.csv")

langs <- get_administration_data() |>
  distinct(language) |>
  pull(language)
cat("Wordbank has", length(langs), "languages. Surveying each...\n\n")

survey_one <- function(lang) {
  ## NOTE: don't filter on production — many languages have production
  ## as NA in the admin table (item-level pulls are required). We just
  ## count admins; whether production data is available is checked
  ## downstream when we actually extract item data.
  ad <- tryCatch(
    get_administration_data(language = lang) |>
      mutate(form_simple = case_when(
        grepl("WS$|WordsAndSentences", form) ~ "WS",
        grepl("WG$|WordsAndGestures", form)  ~ "WG",
        TRUE                                  ~ form
      )),
    error = function(e) NULL
  )
  if (is.null(ad) || nrow(ad) == 0) return(NULL)
  total_kids <- length(unique(ad$child_id))
  admins_per_kid <- ad |> count(child_id) |> pull(n)
  long_kids <- sum(admins_per_kid >= 2)
  age_range <- range(ad$age, na.rm = TRUE)
  forms_present <- ad |> count(form_simple) |>
    transmute(s = sprintf("%s(%d)", form_simple, n)) |>
    pull(s) |> paste(collapse = " ")
  # rough production obs estimate: each admin's item count would require
  # pulling item-level data, which is expensive across all langs. So we
  # report admins instead and item counts only for the ones we'll actually
  # run.
  tibble(language   = lang,
          total_kids = total_kids,
          long_kids  = long_kids,
          admins     = nrow(ad),
          forms      = forms_present,
          age_min    = age_range[1],
          age_max    = age_range[2])
}

results <- bind_rows(lapply(langs, function(l) {
  cat(sprintf("  %s ... ", l))
  r <- survey_one(l)
  if (is.null(r)) { cat("skip\n"); return(NULL) }
  cat(sprintf("kids=%d long=%d admins=%d\n",
              r$total_kids, r$long_kids, r$admins))
  r
}))

results <- results |> arrange(desc(long_kids))
cat("\n=== Languages with >= 100 longitudinal kids ===\n")
print(results |> filter(long_kids >= 100), n = 40)

cat("\n=== All languages, by longitudinal coverage ===\n")
print(results, n = 50)

write_csv(results, OUT)
cat(sprintf("\nWrote %s\n", OUT))
