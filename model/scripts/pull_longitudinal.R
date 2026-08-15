## Pull English + Norwegian + Japanese longitudinal item-level data from Wordbank,
## across the WG- and WS-family forms. Combining forms extends coverage
## (WG starts at 8 months) and gives every child the maximal item set
## seen across their admin history.
##
## Usage:   Rscript model/scripts/pull_longitudinal.R
## Outputs: fits/long_items.rds
##
## 2026-08-15: aligned with the acceleration repo's clean extraction
## (studies/bayes_long/00_prepare_bundles.R there); see
## journal/paper_models_provenance.md (08-15 currency flags). Three fixes:
##   1. CHILD KEY: `ckey` = dataset_name::study_internal_id. Wordbank
##      `child_id` fails to link a child's WG and WS administrations in
##      some datasets (silently splitting cross-form kids into fake
##      single-form kids; cost Marchman's entire WG arm and ~178
##      Norwegian kids). The legacy `child_id` column is retained for
##      other consumers of long_items.rds; the longitudinal prepare_*
##      scripts key on `ckey`.
##   2. LONGITUDINAL RULE: keep a child with >= 2 admins at distinct ages
##      ACROSS forms (was: >= 2 on the same form, which excluded exactly
##      the cross-form kids the key fix rescues).
##   3. Monolingual-TD exclusions (Bilingual origins, Edgin, Byers) now
##      applied here, once, instead of downstream; `uni_lemma` +
##      `form_type` carried through for cross-form item harmonization in
##      the prepare_* scripts; item responses with NA value count as
##      not-produced (was: dropped).

suppressPackageStartupMessages({
  library(wordbankr)
  library(dplyr)
  library(tidyr)
})

source("model/R/config.R")
source("model/R/helpers.R")

LANGUAGES        <- c("English (American)", "Norwegian", "Japanese")
FORM_KEEP        <- c("WG", "WGProd", "WGProdShort", "WGShort", "WS", "WSShort")
EXCLUDE_DATASETS <- c("Edgin", "Byers")   # clinical samples

all_rows <- list()
for (lang in LANGUAGES) {
  ## Admin table brings study_internal_id (the within-dataset child key that
  ## actually links WG<->WS admins) + dataset fields for the TD exclusions.
  message(sprintf("Fetching administrations: %s ...", lang))
  ad <- get_administration_data(language = lang, include_study_internal_id = TRUE) %>%
    filter(!grepl("Bilingual", dataset_origin_name, ignore.case = TRUE),
           !(dataset_name %in% EXCLUDE_DATASETS),
           form %in% FORM_KEEP) %>%
    transmute(data_id, form, age, dataset_name,
              child_id,
              ## fall back to the wordbank child_id where a dataset carries no
              ## study_internal_id (a missing internal id must not collapse a
              ## whole dataset onto one fake child)
              ckey = if_else(is.na(study_internal_id) | study_internal_id == "",
                             paste(dataset_name, "wb", child_id, sep = "::"),
                             paste(dataset_name, study_internal_id, sep = "::")))
  message(sprintf("    %d admins across %d datasets", nrow(ad),
                  n_distinct(ad$dataset_name)))

  for (form in intersect(FORM_KEEP, unique(ad$form))) {
    message(sprintf("Fetching items: %s / %s ...", lang, form))
    d <- tryCatch(
      get_instrument_data(language = lang, form = form,
                          administration_info = FALSE,
                          item_info = TRUE),
      error = function(e) {
        message(sprintf("    skipped: %s", conditionMessage(e)))
        NULL
      })
    if (is.null(d) || nrow(d) == 0) next
    message(sprintf("    got %d item rows", nrow(d)))
    all_rows[[length(all_rows) + 1]] <- d %>%
      filter(item_kind == "word") %>%
      ## Use the logical `produces` column (wordbankr 2.0 / dataset v2.0:
      ## `value` now holds raw responses like "yes"/"never", not "produces").
      ## NA -> not produced (acceleration-repo convention; was previously
      ## left NA and dropped downstream).
      mutate(language = lang,
             produces = as.integer(produces %in% TRUE)) %>%
      select(language, data_id, item = item_definition, uni_lemma,
             form_type, lexical_category, produces) %>%
      inner_join(ad %>% select(data_id, form, age, dataset_name, child_id, ckey),
                 by = "data_id") %>%
      select(-data_id)
  }
}

d_long <- bind_rows(all_rows)
message(sprintf("\nCombined %d item rows across %d children (ckey)",
                nrow(d_long), n_distinct(paste(d_long$language, d_long$ckey))))

# Longitudinal children: >= 2 admins at distinct ages, ACROSS forms (a
# child with WG at 12mo and WS at 24mo is longitudinal; the downstream
# prepare_* scripts merge same-age admins and re-apply their own
# MIN_ADMINS filter after QC).
long_child_keys <- d_long %>%
  distinct(language, ckey, age) %>%
  count(language, ckey, name = "n_ages") %>%
  filter(n_ages >= 2) %>%
  distinct(language, ckey)
d_long <- d_long %>% inner_join(long_child_keys, by = c("language", "ckey"))

# Word-frequency probabilities are NOT attached here; the per-language
# prepare_* scripts join fresh CHILDES frequencies (normalize-and-match).
d_long <- d_long %>% mutate(prob = NA_real_)

# ---- Reporting ----
cat("\n=== Per-language summary ===\n")
print(d_long %>%
        group_by(language) %>%
        summarise(rows = n(),
                  children = n_distinct(ckey),
                  admins = n_distinct(paste(ckey, age)),
                  items = n_distinct(item),
                  forms = paste(sort(unique(form)), collapse = "+"),
                  age_range = sprintf("%.0f-%.0f", min(age), max(age)),
                  .groups = "drop"))

cat("\n=== Per-language x dataset children ===\n")
print(d_long %>%
        group_by(language, dataset_name) %>%
        summarise(children = n_distinct(ckey),
                  admins = n_distinct(paste(ckey, age)), .groups = "drop"),
      n = 50)

cat("\n=== Cross-form linkage (children with admins on both WG- and WS-family forms) ===\n")
print(d_long %>%
        mutate(fam = if_else(grepl("^WG", form), "WG", "WS")) %>%
        group_by(language, dataset_name) %>%
        summarise(cross_form = length(intersect(unique(ckey[fam == "WG"]),
                                                unique(ckey[fam == "WS"]))),
                  .groups = "drop"),
      n = 50)

# ---- Save ----
out <- file.path(PATHS$fits_dir, "long_items.rds")
saveRDS(d_long, out)
cat(sprintf("\nSaved %s (%.1f MB)\n", out, file.info(out)$size / 1e6))
