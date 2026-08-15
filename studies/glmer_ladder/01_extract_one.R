## Build tidy per-(kid, age, item) production data frames for glmer fitting.
##
## 2026-08-15 REBASE: this script no longer pulls Wordbank itself — it reads
## `fits/long_items.rds`, the single clean extraction produced by
## model/scripts/pull_longitudinal.R, so the glmer ladder and the Bayesian
## bundles share one data standard (see journal/paper_models_provenance.md,
## 08-15 flags). That standard, ported from the acceleration repo:
##   1. child key = dataset::study_internal_id (`ckey`) — links WG<->WS admins
##      that Wordbank child_id fails to (Marchman's WG arm, ~178 NO kids);
##   2. uni_lemma cross-form item harmonization (option a);
##   3. crater/jump local-outlier QC on /J production proportions;
##   4. monolingual-TD exclusions + wordbankr 2.0 (Redivis) compatibility,
##      all handled upstream in the pull.
## The old in-script Wordbank pull keyed on child_id and tested
## `value == "produces"`, both now wrong (see the provenance flags). The
## one-off 01b_extract_marchman_clean.R is superseded by this rebase.
##
## Usage:
##   Rscript studies/glmer_ladder/01_extract_one.R "English (American)" --by-dataset
##   Rscript studies/glmer_ladder/01_extract_one.R Norwegian
##   Rscript studies/glmer_ladder/01_extract_one.R Japanese
##
## (The language must be in pull_longitudinal.R's LANGUAGES and the pull
##  re-run first if long_items.rds predates it.)
##
## Output (contract unchanged, consumed by 02_fit_one.R):
##   fits/glmer_ladder/data_<slug>.rds — list(
##     df = tibble(child_id, age, item, produces), language, dataset,
##     forms_kept, n_kids, n_admins, n_items, n_obs
##   ) plus new fields: qc_admins_removed, ckey_map (ckey <-> wordbank
##   child_ids, for the demographics join in paper/build_cache.R).

source("model/R/config.R")
source("model/R/helpers.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

args <- commandArgs(trailingOnly = TRUE)
BY_DATASET <- "--by-dataset" %in% args
args <- args[args != "--by-dataset"]
LANG <- if (length(args) >= 1) args[1] else "Japanese"

MIN_ADMINS   <- 2
MIN_ITEM_OBS <- 100

slug <- gsub("[^A-Za-z0-9]+", "_", tolower(LANG))
slug <- gsub("^_+|_+$", "", slug)
out_dir <- file.path(PATHS$fits_dir, "glmer_ladder")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== Extracting longitudinal data for %s (from long_items.rds) ===\n", LANG))

long <- readRDS(file.path(PATHS$fits_dir, "long_items.rds"))
stopifnot(LANG %in% unique(long$language))
d <- long %>% filter(language == LANG)
forms_in_lang <- sort(unique(d$form))
cat(sprintf("Forms present: %s\n", paste(forms_in_lang, collapse = ", ")))

## ckey <-> wordbank child_id map (a clean child can span several wb ids);
## saved per unit for the demographics join downstream.
ckey_map_all <- d %>% distinct(ckey, wb_child_id = child_id, dataset_name)

## Clean child key + cross-form item harmonization (shared helpers).
d <- d %>% mutate(child_id = ckey)
n_ids0 <- n_distinct(d$item)
d <- harmonize_cdi_items(d, form_col = "form_type")
cat(sprintf("Harmonization: %d item_definitions -> %d item ids\n",
            n_ids0, n_distinct(d$item_h)))

## Collapse to one administration per (child, age): WG+WS at the same age and
## same-form retests merge; produced if produced in any admin that month.
df_all <- d %>%
  mutate(item = item_h) %>%
  group_by(child_id, age, item) %>%
  summarise(produces = max(produces), .groups = "drop") %>%
  left_join(d %>% distinct(child_id, dataset_name), by = "child_id")

## Per-unit finishing, mirroring the acceleration prep's order:
## >=MIN_ADMINS -> item filter -> QC local-outlier cleaner -> re->=MIN_ADMINS.
write_unit <- function(du, unit_slug, unit_label) {
  keep <- du %>% distinct(child_id, age) %>% count(child_id) %>%
    filter(n >= MIN_ADMINS) %>% pull(child_id)
  du <- du %>% filter(child_id %in% keep)

  item_keep <- du %>% count(item) %>% filter(n >= MIN_ITEM_OBS) %>% pull(item)
  du <- du %>% filter(item %in% item_keep)

  J_qc <- n_distinct(du$item)
  adm_prop <- du %>% group_by(child_id, age) %>%
    summarise(v = sum(produces) / J_qc, .groups = "drop") %>%
    arrange(child_id, age)
  keep_adm <- adm_prop %>% group_by(child_id) %>%
    group_modify(~ mutate(.x, keep = qc_clean_child(.x$age, .x$v,
                                                    min_admins = MIN_ADMINS))) %>%
    ungroup()
  n_out <- sum(!keep_adm$keep)
  du <- du %>% semi_join(keep_adm %>% filter(keep) %>% select(child_id, age),
                         by = c("child_id", "age"))
  keep2 <- du %>% distinct(child_id, age) %>% count(child_id) %>%
    filter(n >= MIN_ADMINS) %>% pull(child_id)
  du <- du %>% filter(child_id %in% keep2)

  n_kids   <- length(unique(du$child_id))
  n_admins <- du %>% distinct(child_id, age) %>% nrow()
  n_items  <- length(unique(du$item))
  n_obs    <- nrow(du)
  cat(sprintf("\n=== %s ===\n  kids=%d  admins=%d  items=%d  obs=%d  (QC removed %d admins)\n",
              unit_label, n_kids, n_admins, n_items, n_obs, n_out))

  out_rds <- file.path(out_dir, sprintf("data_%s.rds", unit_slug))
  saveRDS(list(
    df         = du %>% select(child_id, age, item, produces),
    language   = LANG,
    dataset    = unit_label,
    forms_kept = forms_in_lang,
    n_kids     = n_kids,
    n_admins   = n_admins,
    n_items    = n_items,
    n_obs      = n_obs,
    qc_admins_removed = n_out,
    ckey_map   = ckey_map_all %>% filter(ckey %in% unique(du$child_id))
  ), out_rds)
  cat(sprintf("Wrote %s\n", out_rds))
}

if (BY_DATASET) {
  for (ds in sort(unique(df_all$dataset_name))) {
    ds_slug <- gsub("[^A-Za-z0-9]+", "_", tolower(ds))
    ds_slug <- gsub("^_+|_+$", "", ds_slug)
    write_unit(df_all %>% filter(dataset_name == ds), ds_slug, ds)
  }
} else {
  write_unit(df_all, slug, LANG)
}
