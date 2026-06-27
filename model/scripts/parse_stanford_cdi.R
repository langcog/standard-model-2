## Parse Stanford CDI files (totlot2 / totlot3 / tlo) into long
## item-level format and emit a reviewable short-code-to-Wordbank-item
## mapping. Combines WG + WS for every cohort so the longitudinal age
## window isn't artificially truncated.
##
## RUN LOCALLY ONLY. Requires wordbankr (Sherlock can't reach it).
##
## Study → public dataset crosswalk:
##   totlot3 = TL3 = AM2018   (IDs 11xxx; WG 13-18 + WS 20-27)
##   tlo     = TLO = FMW2013   (IDs 20xxx; WG@18 + WS@24,30)
##   totlot2 = TL2 = FM2012    (WG 15-19 + WS 14-32; processing-only, no LENA)
##
## Inputs (data/peekbank/<dataset_name>/ — see README label correspondence):
##   adams_marchman_2018/TL3_compiled_WS.csv, TL3_compiled_WG.xlsx       (AM2018)
##   fmw_2013/TLO_18m_WG.xlsx, TLO_24_WS.xlsx, TLO_30m_WS.xlsx           (FMW2013;
##                                          age in filename, 'misc' rows dropped)
##   fernald_marchman_2012/TL2_WG_compiled.xlsx, TL2_WS_compiled.xlsx    (FM2012)
##
## Outputs:
##   data/peekbank/cdi_short_code_map_{ws,wg}.csv
##       Mapping from internal short codes to Wordbank item_definitions.
##       Columns: short, item_definition, status. LOOSE END is hand
##       review of rows where status != "auto_exact"/"manual_disambig".
##   data/intermediates/stanford_cdi_items_long.csv
##       Long format, one row per (lab_subject_id, study, age, form,
##       item, produces). Items are Wordbank item_definitions.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(wordbankr)
})

RAW_DIR <- file.path(PROJECT_ROOT, "data/raw")
INT_DIR <- file.path(PROJECT_ROOT, "data/intermediates")

# -------------------------------------------------------------------- #
# 1.  Manual disambiguator overrides (Marchman-lab CDI conventions).   #
#     These are deterministic, not guesses.                            #
# -------------------------------------------------------------------- #
manual_overrides <- c(
  # Sense-disambiguated items (1=noun/object, 2=second sense)
  chicken1 = "chicken (animal)",   chicken2 = "chicken (food)",
  fish1    = "fish (animal)",      fish2    = "fish (food)",
  can1     = "can (object)",       can2     = "can (auxiliary)",
  water1   = "water (beverage)",   water2   = "water (not beverage)",
  watch1   = "watch (object)",     watch2   = "watch (action)",
  swing1   = "swing (object)",     swing2   = "swing (action)",
  slide1   = "slide (object)",     slide2   = "slide (action)",
  drink1   = "drink (beverage)",   drink2   = "drink (action)",
  clean1   = "clean (action)",     clean2   = "clean (description)",
  dry1     = "dry (action)",       dry2     = "dry (description)",
  orange1  = "orange (food)",      orange2  = "orange (description)",
  work1    = "work (place)",       work2    = "work (action)",
  toy      = "toy (object)",
  dress    = "dress (object)",

  # Closed-class suffixed codes
  ifconn = "if", soconn = "so", andconn = "and", thenconn = "then",
  atprep = "at", byprep = "by", toprep = "to", upprep = "up",
  withprep = "with", underprep = "under",
  notquant = "not", allquant = "all",
  ontopof  = "on top of", nextto = "next to", inside = "inside/in",

  # Compounds / abbreviated multi-word items
  cockddld   = "cock-a-doodle-doo",
  bllybttn   = "belly button",
  playdogh   = "play dough",
  rockingchair = "rocking chair",
  refrigerator = "refrigerator",
  livingrm   = "living room",
  highchar   = "high chair",
  washingmachine = "washing machine",
  applesac   = "apple sauce",
  pentbttr   = "peanut butter",
  greenbns   = "green beans",
  frnchfrs   = "french fries",
  potatchp   = "potato chips",
  hamburgr   = "hamburger",
  strwbrry   = "strawberry",
  babysttr   = "babysitter",
  babysittername = "babysitter's name",
  petname    = "pet's name",
  childname  = "child's own name",
  givemefv   = "gimme five",
  gonnagty   = "gonna get you",
  gopotty    = "go potty",
  thslttlp   = "this little piggy",
  patycak    = "patty cake",
  turnarnd   = "turn around",
  callnphn   = "call (on phone)",
  nghtnght   = "night night",
  shopping   = "shopping",
  thankyou   = "thank you",
  sobig      = "so big",
  peekaboo   = "peekaboo",

  # Sound effects + onomatopoeia (collapsed in TL3, spaced in WB)
  baabaa     = "baa baa",
  choochoo   = "choo choo",
  quackqck   = "quack quack",
  uhoh       = "uh oh",
  woofwoof   = "woof woof",
  yumyum     = "yum yum",

  # Body parts / specials
  buttocks   = "buttocks/bottom*",
  owie       = "owie/boo boo",
  penis      = "penis*",
  vagina     = "vagina*",

  # Places / people with WB *
  church     = "church*",
  daddy      = "daddy*",
  grandma    = "grandma*",
  grandpa    = "grandpa*",
  mommy      = "mommy*",

  # Helping verbs (TL3 omits the slashed forms)
  did        = "did/did ya",
  gonna      = "gonna/going to",
  gotta      = "gotta/got to",
  hafta      = "hafta/have to",
  lemme      = "lemme/let me",
  need       = "need/need to",
  try        = "try/try to",
  wanna      = "wanna/want to",

  # Misc
  shh        = "shh/shush/hush",
  soda       = "soda/pop",
  tissue     = "tissue/kleenex",
  little     = "little (description)",

  # Animals/objects with abbreviated spellings
  alligatr   = "alligator",
  buttrfly   = "butterfly",
  helicptr   = "helicopter",
  motrcycl   = "motorcycle",
  firetrck   = "fire truck",
  bicycle    = "bicycle",
  squirrel   = "squirrel",
  teddyber   = "teddy bear",
  undrpnts   = "underpants",
  yesterdy   = "yesterday",
  sandwich   = "sandwich",
  raisin     = "raisin",
  popsicle   = "popsicle",
  popcorn    = "popcorn",
  necklace   = "necklace",
  pumpkin    = "pumpkin",
  pretzel    = "pretzel",
  pancake    = "pancakes",
  vitamins   = "vitamins",
  sprinklr   = "sprinkler",
  sidewalk   = "sidewalk",
  lawnmowr   = "lawn mower",
  snowsuit   = "snowsuit",
  sneaker   = "sneaker",
  slipper   = "slipper",
  pajamas   = "pajamas",
  mittens   = "mittens",
  napkin    = "napkin",
  garbage   = "garbage",
  hammer    = "hammer",
  scissors  = "scissors",
  glasses   = "glasses",
  picture   = "picture",
  pillow    = "pillow",
  comb      = "comb",
  walker    = "walker",
  gasstation = "gas station",
  outside   = "outside",
  downtown  = "downtown",
  camping   = "camping",
  picnic    = "picnic",
  playgrnd  = "playground",
  morning   = "morning",
  tomorrow  = "tomorrow",
  tonight   = "tonight",
  myself    = "myself",
  yourself  = "yourself",
  another   = "another",
  hate      = "hate",
  finish    = "finish",
  kiss      = "kiss",
  pretend   = "pretend",
  scared    = "scared",
  thirsty   = "thirsty",
  windy     = "windy",
  yellow    = "yellow",
  yucky     = "yucky",

  # Resolved fuzzy matches (auto-mapper got these right by edit distance,
  # promoted to manual_disambig after hand verification)
  grr        = "grrr",
  chocolat   = "chocolate",
  spaghett   = "spaghetti",
  telephon   = "telephone",
  tothbrsh   = "toothbrush",
  breakfst   = "breakfast",
  feet       = "foot"   # auto-mapper picked feed; correct WB target is foot
)

# -------------------------------------------------------------------- #
# 2.  Auto-map remaining short codes to Wordbank item_definitions.     #
# -------------------------------------------------------------------- #
fingerprint <- function(x) gsub("[^a-z0-9]", "", tolower(x))

build_mapping <- function(short_codes, wb_items) {
  fp_wb <- fingerprint(wb_items$item_definition)
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
    fp_s <- fingerprint(s)
    cand <- wb_items$item_definition[fp_wb == fp_s]
    if (length(cand) == 1) {
      out$item_definition[i] <- cand
      out$status[i] <- "auto_exact"
      next
    }
    if (length(cand) > 1) {
      # multiple WB items collapse to the same fingerprint; flag for review
      out$item_definition[i] <- cand[1]
      out$status[i] <- sprintf("auto_ambig:%s", paste(cand, collapse = "|"))
      next
    }
    # nearest fingerprint by edit distance
    d <- as.integer(adist(fp_s, fp_wb)[1, ])
    best <- which.min(d)
    if (d[best] <= 2) {
      out$item_definition[i] <- wb_items$item_definition[best]
      out$status[i] <- sprintf("auto_fuzzy_dist%d", d[best])
    } else {
      out$status[i] <- "unmapped"
    }
  }
  out
}

# -------------------------------------------------------------------- #
# 3.  Identify the vocab-section column ranges per file.               #
#     The vocab section runs from "baabaa" (or first sound effect) up  #
#     to and including the connectives ("andconn".."thenconn"). After  #
#     that come sentence-completion ("usepastB1") and morphology.      #
# -------------------------------------------------------------------- #
vocab_columns <- function(nm, form) {
  start <- which(nm == "baabaa")[1]
  stopifnot(!is.na(start))
  end <- if (form == "WS") which(nm == "thenconn")[1]
         else if (form == "WG") which(nm == "some")[1]
         else NA_integer_
  stopifnot(!is.na(end))
  start:end
}

# -------------------------------------------------------------------- #
# 4.  Read each source file, harmonize meta columns, pivot to long.    #
# -------------------------------------------------------------------- #
read_one <- function(path, form_label, age_override = NA_real_,
                      drop_studies = NULL, force_study = NULL) {
  if (grepl("\\.xlsx$", path)) {
    d <- read_excel(path, sheet = 1, .name_repair = "minimal")
  } else {
    d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  # Make column names unique (TL2 has duplicate `id`, `feet`, `combine`)
  names(d) <- make.unique(names(d), sep = ".")
  # Locate the subject-id column robustly: the column named exactly "id"
  # (TL3_WS/TL2 have it first; TLO files have it at position 4). Fall back
  # to the first column if no exact match.
  id_hits <- which(tolower(names(d)) == "id")
  id_col  <- if (length(id_hits)) names(d)[id_hits[1]] else names(d)[1]
  # Standardize the age column name; TLO files have none → use override.
  age_col <- intersect(c("age_mos", "age_cdi", "agecdi", "age"), names(d))[1]
  if (is.na(age_col)) {
    if (is.na(age_override))
      stop(sprintf("No age col in %s and no age_override given", path))
    d[[".age_override"]] <- age_override
    age_col <- ".age_override"
  }
  # Keep core meta + vocab section
  vocab_idx <- vocab_columns(names(d), form_label)
  meta_cols <- c("id" = id_col,
                 "study" = "study",
                 "form"  = "form",
                 "age"   = age_col,
                 "sex"   = intersect(c("sex/gender", "sex"), names(d))[1])
  if (any(is.na(meta_cols))) stop(sprintf("Missing meta cols in %s: %s",
                                          path, paste(names(meta_cols)[is.na(meta_cols)],
                                                      collapse = ",")))
  d_meta  <- d[, meta_cols, drop = FALSE]
  names(d_meta) <- names(meta_cols)
  d_meta$id <- as.character(d_meta$id)   # ids are numeric in some files, char in others
  # force_study: the TLO sheets label rows "tlo" OR "misc", but ALL of
  # them are FMW2013 children (every "misc" subject id is in the 20xxx
  # FMW2013 range and appears in the tlo set at other ages — verified).
  # So we relabel the whole sheet to the study (not drop "misc", which
  # was a bug that discarded 30/34 of the 30-month admins).
  if (!is.null(force_study)) {
    d_meta$study <- force_study
  } else if (length(drop_studies)) {
    keep_study <- !(tolower(trimws(as.character(d_meta$study))) %in%
                    tolower(drop_studies))
    d        <- d[keep_study, , drop = FALSE]
    d_meta   <- d_meta[keep_study, , drop = FALSE]
  }
  d_vocab <- d[, vocab_idx, drop = FALSE]
  # Coerce every vocab column to character so pivot_longer doesn't choke
  # on mixed dbl/chr (TL2 WG has both "1" strings and numeric 1s).
  d_vocab[] <- lapply(d_vocab, as.character)
  # Drop placeholder rows where the parent never returned the form: any
  # row with zero non-empty vocab cells. (~44% of TL3 rows at older ages.)
  n_filled <- rowSums(!is.na(d_vocab) &
                      vapply(d_vocab, function(c) nzchar(trimws(c)), logical(nrow(d_vocab))))
  keep <- n_filled > 0
  d_meta  <- d_meta [keep, , drop = FALSE]
  d_vocab <- d_vocab[keep, , drop = FALSE]

  # Pivot to long, then code PRODUCTION. Critical form difference:
  #   WS: cells are 1 = produces, NA = not. produces = (raw == "1").
  #   WG: cells are 1 = understands, 2 = understands AND says (produces),
  #       NA = neither. So production on WG is raw == "2" ONLY; coding
  #       "1" as produces (the old bug) counted comprehension as
  #       production and inflated WG vocab ~4x.
  prod_codes <- if (form_label == "WG") c("2") else
                 c("1", "produces", "yes", "y", "x")
  d_long <- bind_cols(d_meta, d_vocab) %>%
    pivot_longer(-c(id, study, form, age, sex),
                 names_to = "short", values_to = "raw") %>%
    mutate(rawn = tolower(trimws(as.character(raw))),
           produces = if_else(!is.na(raw) & rawn %in% prod_codes, 1L, 0L),
           # WG also records comprehension (1 = understands, 2 = understands+says);
           # WS is production-only, so comprehension is undefined (NA).
           comprehends = if (form_label == "WG")
             if_else(!is.na(raw) & rawn %in% c("1", "2"), 1L, 0L) else NA_integer_) %>%
    select(-raw, -rawn)
  attr(d_long, "form_label") <- form_label
  d_long
}

cat("Reading all Stanford CDI source files...\n")
report <- function(tag, d) cat(sprintf("  %-12s %5d rows (subjects: %d, admins: %d, ages: %s)\n",
            tag, nrow(d), n_distinct(d$id),
            n_distinct(paste(d$id, d$age)),
            paste(range(suppressWarnings(as.integer(d$age)), na.rm = TRUE),
                  collapse = "-")))
# FM2012 (TL2 / totlot2): WG + WS  [raw dir renamed -> fernald_marchman_2012/]
tl2_wg <- read_one(file.path(RAW_DIR, "FM2012/TL2_WG_compiled.xlsx"), "WG"); report("TL2 WG", tl2_wg)
tl2_ws <- read_one(file.path(RAW_DIR, "FM2012/TL2_WS_compiled.xlsx"), "WS"); report("TL2 WS", tl2_ws)
# AM2018 (TL3 / totlot3): WS (existing csv) + WG  [raw dir -> adams_marchman_2018/]
tl3_ws <- read_one(file.path(RAW_DIR, "AM2018/TL3_compiled_WS.csv"),  "WS"); report("TL3 WS", tl3_ws)
tl3_wg <- read_one(file.path(RAW_DIR, "AM2018/TL3_compiled_WG.xlsx"), "WG"); report("TL3 WG", tl3_wg)
# FMW2013 (TLO): WG@18, WS@24, WS@30 — age comes from the filename.  [raw dir -> fmw_2013/]
# force_study="tlo": every row in these sheets is an FMW2013 child
# (the "misc" label is a within-study annotation, not a different
# study), so keep them all.
tlo_wg18 <- read_one(file.path(RAW_DIR, "FMW2013/TLO_18m_WG.xlsx"), "WG", age_override = 18, force_study = "tlo"); report("TLO WG18", tlo_wg18)
tlo_ws24 <- read_one(file.path(RAW_DIR, "FMW2013/TLO_24_WS.xlsx"),  "WS", age_override = 24, force_study = "tlo"); report("TLO WS24", tlo_ws24)
tlo_ws30 <- read_one(file.path(RAW_DIR, "FMW2013/TLO_30m_WS.xlsx"), "WS", age_override = 30, force_study = "tlo"); report("TLO WS30", tlo_ws30)

# -------------------------------------------------------------------- #
# 5.  Build mappings (one per form), apply, and emit outputs.          #
# -------------------------------------------------------------------- #
wb_ws <- get_item_data(language = "English (American)", form = "WS") %>%
  filter(item_kind == "word") %>% select(item_id, item_definition, category)
wb_wg <- get_item_data(language = "English (American)", form = "WG") %>%
  filter(item_kind == "word") %>% select(item_id, item_definition, category)

ws_short <- unique(c(tl3_ws$short, tl2_ws$short, tlo_ws24$short, tlo_ws30$short))
wg_short <- unique(c(tl2_wg$short, tl3_wg$short, tlo_wg18$short))

map_ws <- build_mapping(ws_short, wb_ws)
map_wg <- build_mapping(wg_short, wb_wg)

cat("\nWS mapping status:\n"); print(table(sub(":.*", ":<...>", map_ws$status)))
cat("WG mapping status:\n");   print(table(sub(":.*", ":<...>", map_wg$status)))

write_csv(map_ws, file.path(INT_DIR, "cdi_short_code_map_ws.csv"))
write_csv(map_wg, file.path(INT_DIR, "cdi_short_code_map_wg.csv"))
cat(sprintf("\nWrote cdi_short_code_map_ws.csv (%d rows) + _wg.csv (%d rows)\n",
            nrow(map_ws), nrow(map_wg)))

# Apply mapping & filter to mapped items only (drop unmapped)
apply_map <- function(d_long, map_df) {
  d_long %>%
    inner_join(map_df %>% filter(!is.na(item_definition)) %>%
                 select(short, item_definition, status),
               by = "short")
}

cdi_long <- bind_rows(
  apply_map(tl2_wg,   map_wg) %>% mutate(form = "WG", source_file = "TL2_WG_compiled.xlsx"),
  apply_map(tl2_ws,   map_ws) %>% mutate(form = "WS", source_file = "TL2_WS_compiled.xlsx"),
  apply_map(tl3_ws,   map_ws) %>% mutate(form = "WS", source_file = "TL3_compiled_WS.csv"),
  apply_map(tl3_wg,   map_wg) %>% mutate(form = "WG", source_file = "TL3_compiled_WG.xlsx"),
  apply_map(tlo_wg18, map_wg) %>% mutate(form = "WG", source_file = "TLO_18m_WG.xlsx"),
  apply_map(tlo_ws24, map_ws) %>% mutate(form = "WS", source_file = "TLO_24_WS.xlsx"),
  apply_map(tlo_ws30, map_ws) %>% mutate(form = "WS", source_file = "TLO_30m_WS.xlsx")
) %>%
  rename(lab_subject_id = id, item = item_definition) %>%
  mutate(age = suppressWarnings(as.integer(age))) %>%
  filter(!is.na(age), !is.na(lab_subject_id), nzchar(lab_subject_id)) %>%
  select(lab_subject_id, study, age, form, sex, item, produces, comprehends,
         short, mapping_status = status, source_file)

# -------------------------------------------------------------------- #
# 5b. Drop corrupted FM2012 (totlot2) 24mo WS admins.                  #
#     In TL2_WS_compiled.xlsx these 11 admins have biologically        #
#     impossible production profiles: ~0% of the 100 easiest words     #
#     marked, yet 27-36 of the 40 hardest items (function words, rare  #
#     nouns, late grammar/complexity) marked — and most sit between a  #
#     500+ word 21mo and a 600+ word 30mo on the same form. The source #
#     cells genuinely contain only these sparse hard-item marks (no    #
#     column shift — a shift would preserve the ~550 mark count, not   #
#     collapse it to 17-276 — and there is no alternate source file),  #
#     so the true 24mo vocab is unrecoverable. We drop these admins    #
#     rather than emit fake mid-trajectory dips. Localized to 24mo WS; #
#     other totlot2 ages/kids parse correctly.                         #
fm2012_bad_ws24 <- c("10006", "10008", "10013", "10016", "10024", "10025",
                     "10031", "10032", "10040", "10043", "10053")
n_before <- n_distinct(paste(cdi_long$lab_subject_id, cdi_long$age, cdi_long$form))
cdi_long <- cdi_long %>%
  filter(!(study == "totlot2" & age == 24L & form == "WS" &
           lab_subject_id %in% fm2012_bad_ws24))
n_after <- n_distinct(paste(cdi_long$lab_subject_id, cdi_long$age, cdi_long$form))
cat(sprintf("Dropped %d corrupted FM2012 24mo WS admins (misaligned source data).\n",
            n_before - n_after))

cat(sprintf("\nLong-format CDI: %d rows (%d subjects x %d admins x %d items)\n",
            nrow(cdi_long), n_distinct(cdi_long$lab_subject_id),
            n_distinct(paste(cdi_long$lab_subject_id, cdi_long$age, cdi_long$form)),
            n_distinct(cdi_long$item)))

write_csv(cdi_long, file.path(INT_DIR, "stanford_cdi_items_long.csv"))
cat(sprintf("Wrote stanford_cdi_items_long.csv\n"))
