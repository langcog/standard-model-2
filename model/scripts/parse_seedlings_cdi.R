## Parse the SEEDLingS WG CDI export (cdi_ht_raw_temp.csv from
## Dong & Bergelson 2026) into the long item-level format expected by
## prepare_seedlings.R.
##
## Schema in: wide format. 398 `Talk_<item>` columns whose values are
##   "Understands and Says" / "Understands" / NA.
##   "Understands and Says" -> produces = 1; everything else -> 0.
## 31 `Understand_<phrase>` columns and 60 `Gestures_<x>` columns are
## ignored here -- we only model production.
##
## subj alignment caveat: the source `subj` column uses bare integer
## strings ("1".."44") while the public Egan-Dailey `lena_data.csv`
## uses zero-padded "01".."44". We coerce both sides to integer and
## back to padded string so they line up. (Mike has flagged this as
## pending verification; if alignment is wrong we just rerun.)
##
## Inputs:
##   data/raw_data/seedlings/cdi_ht_raw_temp.csv
##
## Output:
##   data/raw_data/seedlings/cdi_items_long.csv
##   data/raw_data/seedlings/cdi_seedlings_short_code_map.csv

source("model/R/config.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(wordbankr)
})

OUT_DIR <- file.path(PROJECT_ROOT, "data/raw_data/seedlings")

# -------------------------------------------------------------------- #
# Manual overrides for sense-disambiguated items (WG form).            #
# Format: <Talk_short> -> <Wordbank item_definition>                   #
# -------------------------------------------------------------------- #
manual_overrides <- c(
  Talk_chicken_animal = "chicken (animal)",
  Talk_chicken_food   = "chicken (food)",
  Talk_fish_animal    = "fish (animal)",
  Talk_fish_food      = "fish (food)",
  Talk_water_food     = "water (beverage)",
  Talk_water_outside  = "water (not beverage)",
  Talk_swing_outside  = "swing (object)",
  Talk_swing_action   = "swing (action)",
  Talk_slide_outside  = "slide (object)",
  Talk_slide_action   = "slide (action)",
  Talk_watch_object   = "watch (object)",
  Talk_watch_action   = "watch (action)",
  Talk_clean_action   = "clean (action)",
  Talk_clean_descriptor = "clean (description)",
  Talk_drink_beverage = "drink (beverage)",
  Talk_drink_action   = "drink (action)",
  Talk_dry            = "dry (description)",  # WG has only "dry (description)"
  Talk_orange         = "orange (food)",      # WG only has orange (food)
  Talk_work           = "work (place)",       # WG only has work (place)
  Talk_toy            = "toy (object)",
  Talk_dress          = "dress (object)",
  Talk_little         = "little (description)",

  # Compounds / multi-word items
  Talk_cockadoodledoo = "cock-a-doodle-doo",
  Talk_bellybutton    = "belly button",
  Talk_booboo         = "owie/boo boo",
  Talk_babysitter     = "babysitter",
  Talk_babysittername = "babysitter's name",
  Talk_childname      = "child's own name",
  Talk_pattycake      = "patty cake",
  Talk_peekaboo       = "peekaboo",
  Talk_nightnight     = "night night",
  Talk_thankyou       = "thank you",
  Talk_dont           = "don't",
  Talk_wanna          = "wanna/want to",
  Talk_shh            = "shh/shush/hush",

  # Sound effects (collapsed in cdi_ht, spaced in WB)
  Talk_baabaa         = "baa baa",
  Talk_choochoo       = "choo choo",
  Talk_quackquack     = "quack quack",
  Talk_uhoh           = "uh oh",
  Talk_woofwoof       = "woof woof",
  Talk_yumyum         = "yum yum",

  # Names with WB asterisks
  Talk_daddy          = "daddy*",
  Talk_grandma        = "grandma*",
  Talk_grandpa        = "grandpa*",
  Talk_mommy          = "mommy*",
  Talk_church         = "church*",

  # Compound multi-word common items
  Talk_firetruck       = "fire truck",
  Talk_icecream        = "ice cream",
  Talk_teddybear       = "teddy bear",
  Talk_highchair       = "high chair",
  Talk_livingroom      = "living room",
  Talk_rockingchair    = "rocking chair",

  # Closed class
  Talk_under          = "under",
  Talk_up             = "up",
  Talk_in             = "in",
  Talk_inside         = "inside",
  Talk_not            = "not",
  Talk_all            = "all",
  Talk_slide          = "slide (object)",   # WG only has slide (object)
  Talk_grr            = "grrr"
)

# Short codes that are NOT vocabulary items (parent-meta questions
# like "does your child parrot?"). Drop these silently.
non_vocab_drops <- c("Talk_parroting", "Talk_labeling")

# -------------------------------------------------------------------- #
# Helper: fingerprint normalization                                    #
# -------------------------------------------------------------------- #
fingerprint <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# -------------------------------------------------------------------- #
# 1. Read raw                                                          #
# -------------------------------------------------------------------- #
d <- read.csv(file.path(OUT_DIR, "cdi_ht_raw_temp.csv"),
              check.names = FALSE, stringsAsFactors = FALSE)
# First column is an unnamed row index from R's write.csv; drop it
if (names(d)[1] == "") d <- d[, -1, drop = FALSE]
cat(sprintf("Raw rows: %d, cols: %d\n", nrow(d), ncol(d)))

talk_cols <- setdiff(grep("^Talk_", names(d), value = TRUE), non_vocab_drops)
cat(sprintf("Talk_ vocab columns: %d (dropped %d non-vocab)\n",
            length(talk_cols), length(non_vocab_drops)))

# -------------------------------------------------------------------- #
# 2. Build mapping Talk_<short> -> WG item_definition                  #
# -------------------------------------------------------------------- #
wb_wg <- get_item_data(language = "English (American)", form = "WG") %>%
  filter(item_kind == "word") %>%
  select(item_definition, category)
fp_wb <- fingerprint(wb_wg$item_definition)

build_mapping <- function(short_codes) {
  out <- tibble(short = short_codes,
                item_definition = NA_character_,
                status = NA_character_)
  for (i in seq_along(short_codes)) {
    s <- short_codes[i]
    if (s %in% names(manual_overrides)) {
      out$item_definition[i] <- manual_overrides[[s]]
      out$status[i] <- "manual_disambig"
      next
    }
    # strip Talk_ prefix and try to match
    bare <- sub("^Talk_", "", s)
    cand <- wb_wg$item_definition[fp_wb == fingerprint(bare)]
    if (length(cand) == 1) {
      out$item_definition[i] <- cand
      out$status[i] <- "auto_exact"
      next
    }
    if (length(cand) > 1) {
      out$item_definition[i] <- cand[1]
      out$status[i] <- sprintf("auto_ambig:%s", paste(cand, collapse="|"))
      next
    }
    # try fuzzy
    d_ed <- as.integer(adist(fingerprint(bare), fp_wb)[1, ])
    best <- which.min(d_ed)
    if (d_ed[best] <= 2) {
      out$item_definition[i] <- wb_wg$item_definition[best]
      out$status[i] <- sprintf("auto_fuzzy_dist%d", d_ed[best])
    } else {
      out$status[i] <- "unmapped"
    }
  }
  out
}

map_df <- build_mapping(talk_cols)
cat("\nMapping status counts:\n")
print(table(sub(":.*", ":<...>", map_df$status)))
cat(sprintf("Unmapped Talk_ short codes: %d\n",
            sum(map_df$status == "unmapped")))
if (any(map_df$status == "unmapped")) {
  cat("Unmapped short codes:\n")
  print(map_df %>% filter(status == "unmapped") %>% pull(short))
}
if (any(grepl("^auto_fuzzy", map_df$status))) {
  cat("\nFuzzy matches (review):\n")
  print(map_df %>% filter(grepl("^auto_fuzzy", status)) %>% as.data.frame())
}

write_csv(map_df, file.path(OUT_DIR, "cdi_seedlings_short_code_map.csv"))

# -------------------------------------------------------------------- #
# 3. Pivot to long & code production                                   #
# -------------------------------------------------------------------- #
# Sequential subject-id mapping: bare integer "1" .. "44" -> "01" .. "44"
# (drops alphanumeric subj like "08semns1" — those are duplicate/alt
# admins per the upstream README convention; keeping only canonical kids)
d_meta <- d %>%
  mutate(subject_id_int = suppressWarnings(as.integer(subj)),
         subject_id = sprintf("%02d", subject_id_int)) %>%
  # Keep only Egan-Dailey final-sample subjects (drops subj 24, who was
  # enrolled but excluded from the published sample)
  filter(!is.na(subject_id_int),
         subject_id_int >= 1, subject_id_int <= 99,
         SeedlingsFinalSample == "Y") %>%
  select(subject_id, age = month, raw_subj = subj, all_of(talk_cols))

cat(sprintf("\nAfter subject-id filter: %d admins (%d unique subjects)\n",
            nrow(d_meta), n_distinct(d_meta$subject_id)))

d_long <- d_meta %>%
  pivot_longer(-c(subject_id, age, raw_subj),
               names_to = "short", values_to = "raw") %>%
  inner_join(map_df %>% filter(!is.na(item_definition)) %>%
               select(short, item_definition, mapping_status = status),
             by = "short") %>%
  mutate(produces = as.integer(!is.na(raw) & raw == "Understands and Says"),
         form = "WG") %>%
  rename(item = item_definition) %>%
  select(subject_id, age, form, item, produces, raw, short, mapping_status)

cat(sprintf("Long-format rows: %d (%d admins x %d items)\n",
            nrow(d_long),
            n_distinct(paste(d_long$subject_id, d_long$age)),
            n_distinct(d_long$item)))

# Drop empty-form admins (parent did not return CDI)
admin_filled <- d_long %>%
  group_by(subject_id, age) %>%
  summarise(any_filled = any(!is.na(raw) & nzchar(raw)), .groups = "drop")
keep <- admin_filled %>% filter(any_filled) %>%
  select(subject_id, age)
d_long <- d_long %>% inner_join(keep, by = c("subject_id","age"))
cat(sprintf("After empty-admin filter: %d rows (%d admins, %d subjects)\n",
            nrow(d_long),
            n_distinct(paste(d_long$subject_id, d_long$age)),
            n_distinct(d_long$subject_id)))

# Cross-check produces totals
agg <- d_long %>% group_by(subject_id, age) %>%
  summarise(n_produces = sum(produces), .groups = "drop")
cat("\nProduces by age (median across kids):\n")
print(agg %>% group_by(age) %>%
      summarise(n_admins = n(),
                median_n = median(n_produces),
                p25 = quantile(n_produces, 0.25),
                p75 = quantile(n_produces, 0.75),
                .groups = "drop"))

write_csv(d_long %>% select(subject_id, age, form, item, produces),
          file.path(OUT_DIR, "cdi_items_long.csv"))
cat(sprintf("\nWrote %s\n", file.path(OUT_DIR, "cdi_items_long.csv")))
