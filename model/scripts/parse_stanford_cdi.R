## Parse Stanford CDI files (totlot2 / totlot3) into long item-level
## format and emit a reviewable short-code-to-Wordbank-item mapping.
##
## RUN LOCALLY ONLY. Requires wordbankr (Sherlock can't reach it).
##
## Inputs:
##   data/raw_data/peekbank/TL3_compiled_WS.csv   (csv; ~22 + 25 mo, WS form)
##   data/raw_data/peekbank/TL2_WS_compiled.xlsx  (~18 mo, WS form)
##   data/raw_data/peekbank/TL2_WG_compiled.xlsx  (~16 mo, WG form)
##
## Outputs:
##   data/raw_data/peekbank/cdi_short_code_map_ws.csv
##   data/raw_data/peekbank/cdi_short_code_map_wg.csv
##       Mapping from internal short codes to Wordbank item_definitions.
##       Columns: short, item_definition, status. All entries are
##       used in production at the moment; LOOSE END is hand review of
##       the rows where status != "auto_exact" or "manual_disambig".
##   data/raw_data/peekbank/stanford_cdi_items_long.csv
##       Long format, one row per (lab_subject_id, study, age, form,
##       item, produces). Items are Wordbank item_definitions.

source("model/R/config.R")
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(wordbankr)
})

OUT_DIR <- file.path(PROJECT_ROOT, "data/raw_data/peekbank")

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
read_one <- function(path, form_label) {
  if (grepl("\\.xlsx$", path)) {
    d <- read_excel(path, sheet = 1, .name_repair = "minimal")
  } else {
    d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  # Make column names unique (TL2 has duplicate `id`, `feet`, `combine`)
  names(d) <- make.unique(names(d), sep = ".")
  # Standardize the age column name
  age_col <- intersect(c("age_mos", "age_cdi", "agecdi", "age"), names(d))[1]
  if (is.na(age_col)) stop(sprintf("No age col in %s", path))
  # Keep core meta + vocab section
  vocab_idx <- vocab_columns(names(d), form_label)
  meta_cols <- c("id" = names(d)[1],   # first id column
                 "study" = "study",
                 "form"  = "form",
                 "age"   = age_col,
                 "sex"   = intersect(c("sex/gender", "sex"), names(d))[1])
  if (any(is.na(meta_cols))) stop(sprintf("Missing meta cols in %s: %s",
                                          path, paste(names(meta_cols)[is.na(meta_cols)],
                                                      collapse = ",")))
  d_meta  <- d[, meta_cols, drop = FALSE]
  names(d_meta) <- names(meta_cols)
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

  # Pivot to long: produces = 1 if cell is "1" (or "produces"); else 0
  d_long <- bind_cols(d_meta, d_vocab) %>%
    pivot_longer(-c(id, study, form, age, sex),
                 names_to = "short", values_to = "raw") %>%
    mutate(produces = case_when(
      is.na(raw)                    ~ 0L,
      tolower(trimws(as.character(raw))) %in% c("1", "produces", "yes", "y", "x") ~ 1L,
      TRUE                          ~ 0L
    )) %>%
    select(-raw)
  attr(d_long, "form_label") <- form_label
  d_long
}

cat("Reading TL3 (WS), TL2 WS, TL2 WG...\n")
tl3_ws <- read_one(file.path(OUT_DIR, "TL3_compiled_WS.csv"),     "WS")
tl2_ws <- read_one(file.path(OUT_DIR, "TL2_WS_compiled.xlsx"),    "WS")
tl2_wg <- read_one(file.path(OUT_DIR, "TL2_WG_compiled.xlsx"),    "WG")
cat(sprintf("  TL3 WS: %d rows (subjects: %d, admins: %d)\n",
            nrow(tl3_ws), n_distinct(tl3_ws$id),
            n_distinct(paste(tl3_ws$id, tl3_ws$age))))
cat(sprintf("  TL2 WS: %d rows (subjects: %d, admins: %d)\n",
            nrow(tl2_ws), n_distinct(tl2_ws$id),
            n_distinct(paste(tl2_ws$id, tl2_ws$age))))
cat(sprintf("  TL2 WG: %d rows (subjects: %d, admins: %d)\n",
            nrow(tl2_wg), n_distinct(tl2_wg$id),
            n_distinct(paste(tl2_wg$id, tl2_wg$age))))

# -------------------------------------------------------------------- #
# 5.  Build mappings (one per form), apply, and emit outputs.          #
# -------------------------------------------------------------------- #
wb_ws <- get_item_data(language = "English (American)", form = "WS") %>%
  filter(item_kind == "word") %>% select(item_id, item_definition, category)
wb_wg <- get_item_data(language = "English (American)", form = "WG") %>%
  filter(item_kind == "word") %>% select(item_id, item_definition, category)

ws_short <- unique(c(tl3_ws$short, tl2_ws$short))
wg_short <- unique(tl2_wg$short)

map_ws <- build_mapping(ws_short, wb_ws)
map_wg <- build_mapping(wg_short, wb_wg)

cat("\nWS mapping status:\n"); print(table(sub(":.*", ":<...>", map_ws$status)))
cat("WG mapping status:\n");   print(table(sub(":.*", ":<...>", map_wg$status)))

write_csv(map_ws, file.path(OUT_DIR, "cdi_short_code_map_ws.csv"))
write_csv(map_wg, file.path(OUT_DIR, "cdi_short_code_map_wg.csv"))
cat(sprintf("\nWrote cdi_short_code_map_ws.csv (%d rows) + _wg.csv (%d rows)\n",
            nrow(map_ws), nrow(map_wg)))

# Apply mapping & filter to mapped items only (drop unmapped)
apply_map <- function(d_long, map_df) {
  d_long %>%
    inner_join(map_df %>% filter(!is.na(item_definition)) %>%
                 select(short, item_definition, status),
               by = "short")
}

tl3_ws_m <- apply_map(tl3_ws, map_ws)
tl2_ws_m <- apply_map(tl2_ws, map_ws)
tl2_wg_m <- apply_map(tl2_wg, map_wg)

cdi_long <- bind_rows(
  tl3_ws_m %>% mutate(form = "WS", source_file = "TL3_compiled_WS.csv"),
  tl2_ws_m %>% mutate(form = "WS", source_file = "TL2_WS_compiled.xlsx"),
  tl2_wg_m %>% mutate(form = "WG", source_file = "TL2_WG_compiled.xlsx")
) %>%
  rename(lab_subject_id = id, item = item_definition) %>%
  mutate(age = suppressWarnings(as.integer(age))) %>%
  filter(!is.na(age), !is.na(lab_subject_id), nzchar(lab_subject_id)) %>%
  select(lab_subject_id, study, age, form, sex, item, produces,
         short, mapping_status = status, source_file)

cat(sprintf("\nLong-format CDI: %d rows (%d subjects x %d admins x %d items)\n",
            nrow(cdi_long), n_distinct(cdi_long$lab_subject_id),
            n_distinct(paste(cdi_long$lab_subject_id, cdi_long$age, cdi_long$form)),
            n_distinct(cdi_long$item)))

write_csv(cdi_long, file.path(OUT_DIR, "stanford_cdi_items_long.csv"))
cat(sprintf("Wrote stanford_cdi_items_long.csv\n"))
