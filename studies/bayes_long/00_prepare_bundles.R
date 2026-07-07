## 00_prepare_bundles.R -- by-dataset Stan bundles for the bayes_long ladder.
##
## Corrects two defects in the glmer extraction (studies/glmer_ladder/01_extract_one.R):
##   1. CHILD KEY = study_internal_id within dataset (NOT wordbank child_id, which
##      fails to link a child's WG and WS administrations -- silently splitting
##      cross-form kids into fake single-form kids; cost ~178 NO kids and
##      Marchman's entire WG arm).
##   2. ITEM HARMONIZATION by uni_lemma cross-link (option a): a WG item and WS
##      item sharing an *unambiguous* uni_lemma are the same latent item; else
##      keep form-specific by item_definition. (Handles in/inside -> inside/in.)
##
## Also: monolingual-TD filters, anchor a0 = 18, full Norwegian, production items.
## Output: fits/bayes_long/bundle_<slug>.rds  (stan_data arrays + index maps + meta)
##
## Usage:  Rscript studies/bayes_long/00_prepare_bundles.R [slug ...]   (default: all)

suppressPackageStartupMessages({
  library(wordbankr); library(dplyr); library(tidyr); library(here)
})

UNITS <- tibble::tribble(
  ~slug,       ~language,             ~dataset,
  "thal",      "English (American)",  "Thal",
  "smith",     "English (American)",  "Smith",
  "marchman",  "English (American)",  "Marchman",
  "norwegian", "Norwegian",           NA,
  "japanese",  "Japanese",            NA
)

FORM_KEEP        <- c("WG","WGProd","WGProdShort","WGShort","WS","WSShort")
EXCLUDE_DATASETS <- c("Edgin","Byers")           # clinical
A0        <- 18                                   # anchor age (explosion milestone)
LOG_H     <- log(365)                             # ~waking hours/month; interpretive offset
MIN_ITEM_OBS <- 100                               # per-unit item filter

OUT_DIR <- here("fits","bayes_long"); dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

args <- commandArgs(trailingOnly=TRUE)
want <- if (length(args)) args else UNITS$slug

## ---- pull one language's admin + production item data (data_id-keyed) ----
pull_language <- function(language) {
  cat(sprintf("\n=== pulling %s ===\n", language))
  ad <- get_administration_data(language, include_study_internal_id=TRUE) |>
    filter(!grepl("Bilingual", dataset_origin_name, ignore.case=TRUE),
           !(dataset_name %in% EXCLUDE_DATASETS),
           form %in% FORM_KEEP)
  forms <- intersect(FORM_KEEP, unique(ad$form))
  it <- bind_rows(lapply(forms, function(f) {
    cat(sprintf("  items: %s ...\n", f))
    get_instrument_data(language=language, form=f,
                        administration_info=FALSE, item_info=TRUE) |>
      filter(item_kind=="word") |>
      transmute(data_id, item_definition, uni_lemma, form_type,
                produces = as.integer(produces %in% TRUE | produces=="produces" | produces==1))
  }))
  # attach admin-level fields (incl. study_internal_id) by data_id
  it <- it |> inner_join(
    ad |> transmute(data_id, age, dataset_name, study_internal_id,
                    ckey = paste(dataset_name, study_internal_id, sep="::")),
    by="data_id")
  list(items=it)
}

## ---- option-a harmonized item id: uni_lemma if unambiguous, else item_definition ----
harmonize_items <- function(it) {
  # a uni_lemma is cross-linkable iff it maps to <=1 item_definition WITHIN EACH
  # form (so in[WG] + inside/in[WS] merge; chicken-animal/chicken-food[WS] do not)
  amb <- it |> distinct(form_type, item_definition, uni_lemma) |>
    filter(!is.na(uni_lemma)) |>
    count(form_type, uni_lemma, name="n_defs") |>
    filter(n_defs > 1) |> pull(uni_lemma) |> unique()
  it |> mutate(item = if_else(!is.na(uni_lemma) & !(uni_lemma %in% amb),
                              paste0("ul:", uni_lemma),
                              paste0("id:", item_definition)))
}

build_bundle <- function(it_unit, slug, label) {
  ## collapse to one obs per (child, age, item): produces if produced in any
  ## admin that month (merges WG+WS at the same age, and same-form retests)
  df <- it_unit |>
    group_by(ckey, age, item) |>
    summarise(produces = max(produces), .groups="drop")

  ## longitudinal: >=2 administrations (distinct child x age)
  keep <- df |> distinct(ckey, age) |> count(ckey) |> filter(n>=2) |> pull(ckey)
  df <- df |> filter(ckey %in% keep)

  ## per-unit item filter
  it_keep <- df |> count(item) |> filter(n>=MIN_ITEM_OBS) |> pull(item)
  df <- df |> filter(item %in% it_keep)

  ## integer indices for Stan
  df <- df |> mutate(admin = paste(ckey, age, sep="@@"))
  child_ix <- tibble(ckey = unique(df$ckey)) |> mutate(ii = row_number())
  admin_ix <- df |> distinct(admin, ckey, age) |>
    left_join(child_ix, by="ckey") |> mutate(aa = row_number())
  item_ix  <- tibble(item = unique(df$item)) |> mutate(jj = row_number())
  obs <- df |> left_join(admin_ix |> select(admin, aa), by="admin") |>
    left_join(item_ix, by="item")

  stan_data <- list(
    N = nrow(obs), A = nrow(admin_ix), I = nrow(child_ix), J = nrow(item_ix),
    aa = obs$aa, jj = obs$jj, y = obs$produces,
    admin_to_child = admin_ix$ii, admin_age = admin_ix$age,
    log_H = LOG_H, a0 = A0)

  ## reporting
  ad_per_kid <- admin_ix |> count(ii) |> pull(n)
  meta <- list(slug=slug, label=label,
               n_kids=nrow(child_ix), n_admins=nrow(admin_ix), n_items=nrow(item_ix),
               n_obs=nrow(obs), age_range=range(admin_ix$age),
               med_admins_per_kid=median(ad_per_kid),
               max_admins_per_kid=max(ad_per_kid))
  cat(sprintf("  [%s] kids=%d admins=%d items=%d obs=%d | age %d-%d | admins/kid med=%d max=%d\n",
              slug, meta$n_kids, meta$n_admins, meta$n_items, meta$n_obs,
              meta$age_range[1], meta$age_range[2], meta$med_admins_per_kid, meta$max_admins_per_kid))

  saveRDS(list(stan_data=stan_data, child_ix=child_ix, item_ix=item_ix,
               admin_ix=admin_ix, meta=meta),
          file.path(OUT_DIR, sprintf("bundle_%s.rds", slug)))
}

## ---- run: pull each needed language once, split into its units ----
units <- UNITS |> filter(slug %in% want)
for (lang in unique(units$language)) {
  pl <- pull_language(lang)
  it <- harmonize_items(pl$items)
  for (i in which(units$language==lang)) {
    u <- units[i,]
    it_u <- if (is.na(u$dataset)) it else filter(it, dataset_name==u$dataset)
    build_bundle(it_u, u$slug, if (is.na(u$dataset)) lang else paste0(lang," / ",u$dataset))
  }
}
cat("\ndone.\n")
