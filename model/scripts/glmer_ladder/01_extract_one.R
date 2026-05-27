## Pull longitudinal CDI data for a given language and save a tidy
## per-(kid, age, item) production data frame for glmer fitting.
##
## We combine WG + WS (and minor form variants) at the item level: the
## same word (e.g., "dog") gets the same item identity regardless of
## which form it was scored on. This means a kid measured at WG (14 mo)
## and WS (24 mo) contributes to both ages on items present on both
## form's checklists.
##
## Usage:
##   Rscript model/scripts/glmer_ladder/01_extract_one.R "English (American)"
##   Rscript model/scripts/glmer_ladder/01_extract_one.R Norwegian
##   Rscript model/scripts/glmer_ladder/01_extract_one.R Japanese
##
## Output: fits/glmer_ladder/data_<lang_slug>.rds
##   contains: df (kid, age, item, produces), child_meta, item_meta,
##              language, forms_kept, n_kids, n_admins, n_items, n_obs

source("model/R/config.R")
suppressPackageStartupMessages({
  library(wordbankr); library(dplyr); library(tidyr); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
LANG <- if (length(args) >= 1) args[1] else "Japanese"

# Form variants we treat as "WG-equivalent" or "WS-equivalent". A few
# languages have minor variants (WGProd, WGShort, WSShort etc.).
# Anything else we drop with a warning — adding new forms here is a
# one-line change.
FORM_KEEP <- c("WG", "WGProd", "WGProdShort", "WGShort",
                "WS", "WSShort")

slug <- gsub("[^A-Za-z0-9]+", "_", tolower(LANG))
out_dir <- file.path(PATHS$fits_dir, "glmer_ladder")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(out_dir, sprintf("data_%s.rds", slug))

cat(sprintf("=== Extracting longitudinal data for %s ===\n", LANG))

## ---- Step 1: admin table; identify longitudinal kids ----
ad <- get_administration_data(language = LANG)
ad_kept <- ad |> filter(form %in% FORM_KEEP)
cat(sprintf("Forms in language: %s\n  Kept: %s\n",
            paste(sort(unique(ad$form)), collapse = ", "),
            paste(sort(unique(ad_kept$form)), collapse = ", ")))
dropped <- setdiff(unique(ad$form), FORM_KEEP)
if (length(dropped) > 0) {
  cat(sprintf("  Dropped (non-standard forms): %s\n",
              paste(dropped, collapse = ", ")))
}

long_kids <- ad_kept |>
  count(child_id) |>
  filter(n >= 2) |>
  pull(child_id)
cat(sprintf("Kids with ≥2 admins (on kept forms): %d\n", length(long_kids)))

ad_long <- ad_kept |> filter(child_id %in% long_kids)
cat(sprintf("Admins from these kids: %d\n", nrow(ad_long)))
cat(sprintf("Age range: [%d, %d] mo\n",
            min(ad_long$age, na.rm = TRUE),
            max(ad_long$age, na.rm = TRUE)))

## ---- Step 2: pull item-level production data ----
## get_instrument_data with administration_info=TRUE attaches age etc.
## We loop over the kept forms because get_instrument_data takes one
## form at a time.
forms_in_lang <- intersect(FORM_KEEP, unique(ad_long$form))
cat(sprintf("Pulling item-level data for forms: %s\n",
            paste(forms_in_lang, collapse = ", ")))

pull_form <- function(f) {
  cat(sprintf("  %s ... ", f))
  t0 <- Sys.time()
  d <- tryCatch(
    get_instrument_data(language = LANG, form = f,
                         administration_info = TRUE,
                         item_info = TRUE),
    error = function(e) { cat(sprintf("ERROR: %s\n", e$message)); NULL }
  )
  if (is.null(d)) return(NULL)
  cat(sprintf("%d rows (%.1f s)\n", nrow(d),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  d
}
items_by_form <- lapply(forms_in_lang, pull_form)
items_by_form <- items_by_form[!sapply(items_by_form, is.null)]

## standardize column set: child_id, age, form, item_definition, item_kind, value (raw)
combine_one <- function(d, form_label) {
  cols <- names(d)
  # column names vary slightly across wordbankr versions: handle both
  child_col <- if ("child_id" %in% cols) "child_id" else "subject_id"
  age_col   <- if ("age" %in% cols) "age" else stop("age column missing")
  item_col  <- if ("item_definition" %in% cols) "item_definition" else "definition"
  uni_col   <- if ("uni_lemma" %in% cols) "uni_lemma" else NA_character_
  kind_col  <- if ("item_kind" %in% cols) "item_kind"
              else if ("type" %in% cols) "type" else NA_character_
  val_col   <- if ("value" %in% cols) "value"
              else if ("response" %in% cols) "response" else NA_character_
  out <- tibble(
    child_id        = d[[child_col]],
    age             = d[[age_col]],
    form            = form_label,
    item_definition = d[[item_col]],
    uni_lemma       = if (!is.na(uni_col)) d[[uni_col]] else NA_character_,
    item_kind       = if (!is.na(kind_col)) d[[kind_col]] else NA_character_,
    value           = if (!is.na(val_col)) d[[val_col]] else NA
  )
  out
}

raw_long <- bind_rows(
  Map(function(d, f) combine_one(d, f),
      items_by_form, forms_in_lang)
)
cat(sprintf("Raw item-level rows: %d\n", nrow(raw_long)))

## ---- Step 3: restrict to production items, score as binary ----
## Production: item_kind == "word" and value indicates "produces".
## value coding varies: in WS, "produces" / NA; in WG, "produces" / "understands" / NA.
## We treat anything matching "produces" (case-insensitive) as 1, else 0.
prod <- raw_long |>
  filter(child_id %in% long_kids) |>
  filter(item_kind == "word" | is.na(item_kind)) |>
  mutate(produces = as.integer(!is.na(value) &
                                tolower(value) == "produces"))

## Some kids may have only WG admins after the merge (since WG comp items dropped) —
## re-check ≥2 admins on production-eligible data.
admins_per_kid_prod <- prod |>
  distinct(child_id, age, form) |>
  count(child_id)
keep_kids <- admins_per_kid_prod |> filter(n >= 2) |> pull(child_id)
prod <- prod |> filter(child_id %in% keep_kids)
cat(sprintf("After production filter, kids retained: %d / %d\n",
            length(keep_kids), length(long_kids)))

## ---- Step 4: build clean df ----
df <- prod |>
  transmute(child_id, age, form,
             item = item_definition,
             produces) |>
  filter(!is.na(produces), !is.na(item), !is.na(age))

## remove items not appearing at least 100 times (avoids tiny-item RE that just adds df)
item_keep <- df |> count(item) |> filter(n >= 100) |> pull(item)
df <- df |> filter(item %in% item_keep)
cat(sprintf("Items kept (≥100 obs): %d\n", length(item_keep)))

## summary
n_kids <- length(unique(df$child_id))
n_admins <- df |> distinct(child_id, age, form) |> nrow()
n_items <- length(unique(df$item))
n_obs <- nrow(df)
cat(sprintf("\n=== Final %s data ===\n  kids=%d  admins=%d  items=%d  obs=%d\n",
            LANG, n_kids, n_admins, n_items, n_obs))

## save
saveRDS(list(
  df       = df,
  language = LANG,
  forms_kept = forms_in_lang,
  n_kids   = n_kids,
  n_admins = n_admins,
  n_items  = n_items,
  n_obs    = n_obs
), out_rds)
cat(sprintf("Wrote %s\n", out_rds))
