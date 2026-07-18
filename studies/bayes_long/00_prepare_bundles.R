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
MIN_ADMINS   <- as.integer(Sys.getenv("MIN_ADMINS", "2"))  # min administrations per child
# Unified local-outlier QC (see clean_child below). A true production trajectory is
# monotone non-decreasing (can't un-produce words) and rises at a plausible rate.
QC_REL_TOL   <- 0.25    # CRATER: an admin >25% below the running peak is a low outlier
QC_PEAK      <- 0.10    #   ...counted only when the running peak reached at least this
QC_DROP      <- 0.05    #   ...and the absolute drop is at least this (ignore report noise)
QC_RATE_MAX  <- 0.40    # JUMP: a rise faster than this (proportion of items / month)...
QC_JUMP_BASE <- 0.10    #   ...launched from a base below this is a high outlier (impossible;
                        #   real fast risers climb from a non-trivial base, so they're spared)

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
                produces = as.integer(produces %in% c(TRUE, 1L, "produces")))  # NA -> 0
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

## ---- unified local-outlier cleaner (per child) ----
## A true production trajectory is monotone non-decreasing (a child cannot un-produce
## words) and its rises are rate-bounded (a child cannot learn ~half the CDI in a
## month). A single administration that violates either is a local outlier -- the same
## mis-keyed WG-comprehension record shows up as an end-crater, a mid spike-then-crash,
## or an impossible end-spike. Greedily remove the single most-severe outlier admin and
## re-check, until the trajectory is clean. Returns a logical KEEP mask over the input
## (age, prop) rows (assumed sorted by age); the >=MIN_ADMINS filter then drops any child
## left with too few waves.
##   - CRATER admin: value >QC_REL_TOL below the running peak (floors QC_PEAK, QC_DROP)
##                   -> the LOW point is the artifact -> remove it.
##   - JUMP admin:   reached via a rise >QC_RATE_MAX/mo from a base <QC_JUMP_BASE
##                   -> the HIGH point is the artifact (impossible speed) -> remove it.
clean_child <- function(age, prop) {
  keep <- rep(TRUE, length(prop))
  repeat {
    ix <- which(keep); if (length(ix) < 2) break
    a <- age[ix]; p <- prop[ix]; peak <- cummax(p)
    jp <- integer(0); js <- numeric(0); cp <- integer(0); cs <- numeric(0)
    for (j in 2:length(p)) {
      drop_amt <- peak[j-1] - p[j]
      rate     <- (p[j] - p[j-1]) / (a[j] - a[j-1])
      if (rate > QC_RATE_MAX && p[j-1] < QC_JUMP_BASE) {                  # jump: high admin j
        jp <- c(jp, ix[j]); js <- c(js, rate)
      } else if (peak[j-1] >= QC_PEAK && drop_amt > QC_REL_TOL*peak[j-1] && drop_amt >= QC_DROP) {
        cp <- c(cp, ix[j]); cs <- c(cs, drop_amt/peak[j-1])              # crater: low admin j
      }
    }
    ## remove JUMPS first: a spike inflates the running peak, making every real point
    ## after it look like a crater -- kill the spike and those induced craters vanish.
    if (length(jp))      keep[jp[which.max(js)]] <- FALSE
    else if (length(cp)) keep[cp[which.max(cs)]] <- FALSE
    else break
    if (sum(keep) < MIN_ADMINS) break                    # will be dropped by admin-count filter
  }
  keep
}

build_bundle <- function(it_unit, slug, label) {
  ## collapse to one obs per (child, age, item): produces if produced in any
  ## admin that month (merges WG+WS at the same age, and same-form retests)
  df <- it_unit |>
    group_by(ckey, age, item) |>
    summarise(produces = max(produces), .groups="drop")

  ## longitudinal: >=2 administrations (distinct child x age)
  keep <- df |> distinct(ckey, age) |> count(ckey) |> filter(n>=MIN_ADMINS) |> pull(ckey)
  df <- df |> filter(ckey %in% keep)

  ## per-unit item filter
  it_keep <- df |> count(item) |> filter(n>=MIN_ITEM_OBS) |> pull(item)
  df <- df |> filter(item %in% it_keep)

  ## QC: unified local-outlier cleaner (see clean_child). Per child, greedily remove
  ## administrations that violate monotone + rate-bounded growth (impossible craters
  ## OR impossible jumps -- both are the mis-keyed WG-comprehension artifact). Removes
  ## the outlier admin, not the whole child; children left with <MIN_ADMINS waves are
  ## then dropped by the admin-count re-filter below.
  ## proportion over the FULL checklist J (not per-administered items): WS is a superset
  ## of WG, so a WG admin scored over its ~396 easy items inflates vs a WS admin over ~680.
  ## Using sum(produces)/J puts both forms on one scale, so a monotonicity violation is a
  ## real decline, not a cross-form item-difficulty artifact. (Rescues ~56 Marchman kids
  ## the per-admin proportion over-excluded; keeps the genuinely-bad craters/tents.)
  J_qc <- n_distinct(df$item)
  prop <- df |> group_by(ckey, age) |> summarise(v = sum(produces)/J_qc, .groups="drop") |> arrange(ckey, age)
  n_admin_before <- nrow(prop)
  keep_adm <- prop |> group_by(ckey) |>
    group_modify(~ mutate(.x, keep = clean_child(.x$age, .x$v))) |> ungroup()
  n_out <- sum(!keep_adm$keep)
  df <- df |> semi_join(filter(keep_adm, keep) |> select(ckey, age), by=c("ckey","age"))
  ## re-apply the >=MIN_ADMINS filter (children may have lost waves to the cleaner)
  keep2 <- df |> distinct(ckey, age) |> count(ckey) |> filter(n >= MIN_ADMINS) |> pull(ckey)
  n_kid_dropped <- length(setdiff(unique(keep_adm$ckey), keep2))
  df <- df |> filter(ckey %in% keep2)

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

  ## reporting (>=3-admin bundles get an _a<N> suffix so they don't clobber the base)
  out_slug <- if (MIN_ADMINS > 2) sprintf("%s_a%d", slug, MIN_ADMINS) else slug
  ad_per_kid <- admin_ix |> count(ii) |> pull(n)
  meta <- list(slug=out_slug, label=label,
               n_kids=nrow(child_ix), n_admins=nrow(admin_ix), n_items=nrow(item_ix),
               n_obs=nrow(obs), age_range=range(admin_ix$age),
               med_admins_per_kid=median(ad_per_kid),
               max_admins_per_kid=max(ad_per_kid),
               qc_admins_removed=n_out, qc_kids_dropped=n_kid_dropped,
               qc_rule=sprintf("local-outlier: crater(>%.2f below peak, floors %.2f/%.2f) | jump(>%.2f/mo from base<%.2f)",
                               QC_REL_TOL, QC_PEAK, QC_DROP, QC_RATE_MAX, QC_JUMP_BASE))
  cat(sprintf("  [%s] kids=%d admins=%d items=%d obs=%d | age %d-%d | admins/kid med=%d max=%d | QC: %d admins removed, %d kids dropped\n",
              out_slug, meta$n_kids, meta$n_admins, meta$n_items, meta$n_obs,
              meta$age_range[1], meta$age_range[2], meta$med_admins_per_kid, meta$max_admins_per_kid,
              n_out, n_kid_dropped))

  saveRDS(list(stan_data=stan_data, child_ix=child_ix, item_ix=item_ix,
               admin_ix=admin_ix, meta=meta),
          file.path(OUT_DIR, sprintf("bundle_%s.rds", out_slug)))
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
