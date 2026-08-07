## 01b_extract_marchman_clean.R -- CLEANED Marchman extraction for the glmer ladder.
##
## Rebuilds Marchman's longitudinal frame applying the three targeted fixes that
## studies/bayes_long/00_prepare_bundles.R identified as defects of the raw
## glmer extraction (01_extract_one.R), so the *standard* glmer ladder can be
## refit on cleaned data and compared to fit_marchman_* (raw):
##
##   1. CHILD KEY = dataset::study_internal_id  (NOT wordbank child_id, which
##      fails to link a child's WG and WS admins -- the raw extraction lost
##      Marchman's entire WG arm, leaving only the 16-30mo WS range).
##   2. ITEM HARMONIZATION by uni_lemma (a WG item and WS item sharing an
##      unambiguous uni_lemma are the same latent item).
##   3. QC crater DROP: remove children whose final vocabulary collapses far
##      below their running peak (Marchman WG comprehension mis-keyed as
##      production: spikes then craters to ~0).
##
## Everything else matches the raw pipeline: monolingual-TD filters, standard
## forms, >=2 administrations, per-unit item filter >=100 obs. Only the three
## fixes above differ, so a raw-vs-clean comparison isolates the cleaning.
##
## Output: fits/glmer_ladder/data_marchman_clean.rds  (same schema as
##         data_marchman.rds -> drop-in for 02_fit_one.R: use slug
##         "marchman_clean").
##
## Usage:  Rscript studies/glmer_ladder/01b_extract_marchman_clean.R

source("model/R/config.R")
suppressPackageStartupMessages({
  library(wordbankr); library(dplyr); library(tidyr)
})

LANGUAGE   <- "English (American)"
DATASET    <- "Marchman"
FORM_KEEP  <- c("WG","WGProd","WGProdShort","WGShort","WS","WSShort")
EXCLUDE    <- c("Edgin","Byers")             # bilingual dropped separately
MIN_ADMINS <- 2                              # match raw glmer longitudinal filter
MIN_ITEM_OBS <- 100                          # match raw per-unit item filter
QC_REL_TOL <- 0.25; QC_PEAK <- 0.10; QC_DROP <- 0.05   # crater rule (bayes_long)

out_dir <- file.path(PATHS$fits_dir, "glmer_ladder")
out_rds <- file.path(out_dir, "data_marchman_clean.rds")

## ---- pull admin + production items, keyed by study_internal_id (fix 1) ----
cat(sprintf("=== pulling %s (for %s) ===\n", LANGUAGE, DATASET))
ad <- get_administration_data(LANGUAGE, include_study_internal_id = TRUE) |>
  filter(!grepl("Bilingual", dataset_origin_name, ignore.case = TRUE),
         !(dataset_name %in% EXCLUDE),
         form %in% FORM_KEEP)
forms <- intersect(FORM_KEEP, unique(ad$form))
cat(sprintf("Forms present: %s\n", paste(forms, collapse = ", ")))

it <- bind_rows(lapply(forms, function(f) {
  cat(sprintf("  items: %s ...\n", f))
  get_instrument_data(language = LANGUAGE, form = f,
                      administration_info = FALSE, item_info = TRUE) |>
    filter(item_kind == "word") |>
    transmute(data_id, item_definition, uni_lemma, form_type,
              produces = as.integer(produces %in% c(TRUE, 1L, "produces")))  # NA -> 0
}))
it <- it |> inner_join(
  ad |> transmute(data_id, age, dataset_name, study_internal_id,
                  ckey = paste(dataset_name, study_internal_id, sep = "::")),
  by = "data_id")

## restrict to the Marchman dataset
it <- it |> filter(dataset_name == DATASET)
cat(sprintf("Marchman rows (pre-collapse): %d over %d data_ids\n",
            nrow(it), dplyr::n_distinct(it$data_id)))

## ---- item harmonization by uni_lemma (fix 2) ----
amb <- it |> distinct(form_type, item_definition, uni_lemma) |>
  filter(!is.na(uni_lemma)) |>
  count(form_type, uni_lemma, name = "n_defs") |>
  filter(n_defs > 1) |> pull(uni_lemma) |> unique()
it <- it |> mutate(item = if_else(!is.na(uni_lemma) & !(uni_lemma %in% amb),
                                  paste0("ul:", uni_lemma),
                                  paste0("id:", item_definition)))

## ---- collapse to one obs per (child, age, item): produced in any admin ----
df <- it |> group_by(ckey, age, item) |>
  summarise(produces = max(produces), .groups = "drop")

## longitudinal: >=2 administrations (distinct child x age)
keep <- df |> distinct(ckey, age) |> count(ckey) |> filter(n >= MIN_ADMINS) |> pull(ckey)
df <- df |> filter(ckey %in% keep)

## per-unit item filter
it_keep <- df |> count(item) |> filter(n >= MIN_ITEM_OBS) |> pull(item)
df <- df |> filter(item %in% it_keep)

## ---- QC crater drop (fix 3) ----
prop <- df |> group_by(ckey, age) |> summarise(v = mean(produces), .groups = "drop")
qc_bad <- prop |> arrange(ckey, age) |> group_by(ckey) |>
  summarise(peak = max(v), last = v[which.max(age)],
            bad = peak >= QC_PEAK & (peak - last) > QC_REL_TOL * peak & (peak - last) >= QC_DROP,
            .groups = "drop") |>
  filter(bad) |> pull(ckey)
df <- df |> filter(!(ckey %in% qc_bad))
cat(sprintf("QC crater-dropped children: %d\n", length(qc_bad)))

## ---- emit in 02_fit_one.R schema (child_id, age, item, produces) ----
df <- df |> transmute(child_id = ckey, age, item, produces) |> arrange(child_id, age, item)

out <- list(
  df = df,
  language = LANGUAGE,
  dataset = DATASET,
  forms_kept = forms,
  n_kids  = dplyr::n_distinct(df$child_id),
  n_admins = nrow(distinct(df, child_id, age)),
  n_items = dplyr::n_distinct(df$item),
  n_obs   = nrow(df)
)
saveRDS(out, out_rds)
cat(sprintf("\nWrote %s\n", out_rds))
cat(sprintf("  kids=%d  items=%d  admins=%d  obs=%d  age=[%d,%d]\n",
            out$n_kids, out$n_items, out$n_admins, out$n_obs,
            min(df$age), max(df$age)))
