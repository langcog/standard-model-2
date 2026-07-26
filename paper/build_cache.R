## build_cache.R — populates outputs/paper/cache/ with the small RDS
## summaries that standard_model.qmd loads in its setup chunk.
##
## Run from project root:
##   Rscript outputs/paper/build_cache.R
##
## Re-run when:
##   - underlying fits change
##   - dataset survey changes
##   - we add languages or studies

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(here)
  library(lme4); library(broom.mixed); library(posterior); library(wordbankr)
})

## By-study units: English fans into its three contributed datasets
## (Thal/Smith/Marchman); Norwegian (Kristoffersen) and Japanese
## (Tsuji+Hagihara) are single bundles kept whole. The Bayesian D/D'
## fits stay pooled-per-language; only the glmer ladder + demographics
## run by-study.
PAPER_LANGS <- c("thal", "smith", "marchman", "norwegian", "japanese")
LANG_LABELS <- c(thal = "English (Thal)", smith = "English (Smith)",
                 marchman = "English (Marchman)",
                 norwegian = "Norwegian", japanese = "Japanese")
# Each unit -> Wordbank language + (optional) dataset_name filter, used to
# pull the right per-child demographics. The English units share one
# English pull, split by dataset_name; NO/JA take the whole language.
UNIT_LANG    <- c(thal = "English (American)", smith = "English (American)",
                  marchman = "English (American)",
                  norwegian = "Norwegian", japanese = "Japanese")
UNIT_DATASET <- c(thal = "Thal", smith = "Smith", marchman = "Marchman",
                  norwegian = NA, japanese = NA)

CACHE <- here("paper", "cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

## (Table 1 is built inline in the tbl-datasets chunk of the manuscript,
##  not cached, so its numbers/names can be edited directly there.)

## ---- 1. glmer-ladder predictions + empirical points -------------
sc <- readRDS(here("fits", "glmer_ladder", "sim_cache.rds"))

# keep all 7 model variants for the four paper languages: the main-text
# figure filters to the log-age models, the supplement to the linear ones.
qt <- sc$qtiles |>
  filter(lang_slug %in% PAPER_LANGS)
ep <- sc$emp |>
  filter(lang_slug %in% PAPER_LANGS)

saveRDS(list(qtiles = qt, emp = ep, langs = PAPER_LANGS,
             lang_labels = LANG_LABELS[PAPER_LANGS]),
        file.path(CACHE, "fig2_glmer_ladder.rds"))
cat(sprintf("Wrote %s (qtiles=%d, emp=%d rows)\n",
            file.path(CACHE, "fig2_glmer_ladder.rds"),
            nrow(qt), nrow(ep)))

## ---- 1b. AIC/BIC model-comparison summary -----------------------
# dAIC / dBIC of every model variant relative to the best model within
# each language (main-text comparison table + inline log-vs-linear range).
aic_summary <- as.data.frame(sc$summ) |>
  filter(lang_slug %in% PAPER_LANGS) |>
  group_by(lang_slug) |>
  mutate(dAIC = AIC - min(AIC), dBIC = BIC - min(BIC)) |>
  ungroup() |>
  select(lang_slug, language, model, AIC, BIC, dAIC, dBIC, df, n_obs, n_kids)
saveRDS(aic_summary, file.path(CACHE, "aic_summary.rds"))
cat(sprintf("Wrote %s (%d rows)\n",
            file.path(CACHE, "aic_summary.rds"), nrow(aic_summary)))

## ---- 3. BLUPs from D_log fits + Wordbank demographics ----------
extract_blups <- function(lang_slug) {
  fp <- here("fits", "glmer_ladder",
             sprintf("fit_%s_D_log.rds", lang_slug))
  if (!file.exists(fp)) return(NULL)
  m  <- readRDS(fp)
  re <- ranef(m)$child
  if (is.null(re)) re <- ranef(m)[[grep("child|kid", names(ranef(m)),
                                         ignore.case = TRUE)[1]]]
  slope_col <- intersect(c("log_age", "la"), colnames(re))[1]
  tibble::tibble(
    child_id = as.integer(rownames(re)),
    xi       = re[, "(Intercept)"],
    zeta     = re[, slope_col],
    lang_slug = lang_slug
  )
}

cat("Extracting BLUPs from D_log fits ...\n")
blups <- bind_rows(lapply(PAPER_LANGS, function(l) {
  b <- extract_blups(l); if (is.null(b)) cat(" - skipping", l, "\n")
  b
}))
cat(sprintf("  total BLUPs: %d kids across %d languages\n",
            nrow(blups), length(unique(blups$lang_slug))))

cat("Pulling Wordbank demographics ...\n")
# Pull each unique language once (keeping dataset_name), then assign to
# units. English's three datasets share one pull, split by dataset_name.
get_demo_lang <- function(lang_label) {
  forms <- c("WS", "WG")
  ad <- bind_rows(lapply(forms, function(f) {
    tryCatch(
      get_administration_data(language = lang_label, form = f,
                              include_demographic_info = TRUE,
                              include_birth_info = TRUE),
      error = function(e) NULL)
  }))
  if (is.null(ad) || nrow(ad) == 0) return(NULL)
  ad |>
    group_by(child_id, dataset_name) |>
    summarise(
      sex          = sex[!is.na(sex)][1],
      birth_order  = birth_order[!is.na(birth_order)][1],
      caregiver_ed = caregiver_education[!is.na(caregiver_education)][1],
      .groups = "drop"
    )
}
demo_by_lang <- setNames(lapply(unique(UNIT_LANG), get_demo_lang),
                         unique(UNIT_LANG))

demos <- bind_rows(lapply(PAPER_LANGS, function(u) {
  d <- demo_by_lang[[ UNIT_LANG[[u]] ]]
  if (is.null(d)) return(NULL)
  if (!is.na(UNIT_DATASET[[u]])) d <- filter(d, dataset_name == UNIT_DATASET[[u]])
  d |> select(-dataset_name) |> mutate(lang_slug = u)
}))
cat(sprintf("  demographics: %d rows\n", nrow(demos)))

# join + recode
recode_demo <- function(df) {
  df |>
    mutate(
      sex = factor(sex, levels = c("Female", "Male")),
      birth_order_2 = case_when(
        birth_order == "First" ~ "First",
        birth_order %in% c("Second", "Third", "Fourth", "Fifth",
                            "Sixth", "Seventh", "Eighth") ~ "Later",
        TRUE ~ NA_character_
      ) |> factor(levels = c("First", "Later")),
      caregiver_ed_3 = case_when(
        caregiver_ed %in% c("None", "Primary", "Some Secondary",
                             "Secondary", "Some College") ~ "<=Some College",
        caregiver_ed == "College" ~ "College",
        caregiver_ed %in% c("Some Graduate", "Graduate") ~ "Graduate+",
        TRUE ~ NA_character_
      ) |> factor(levels = c("<=Some College", "College", "Graduate+"))
    )
}
joined <- blups |>
  left_join(demos, by = c("lang_slug", "child_id")) |>
  recode_demo()
saveRDS(joined, file.path(CACHE, "blups_demographics.rds"))
cat(sprintf("Wrote %s (%d rows)\n",
            file.path(CACHE, "blups_demographics.rds"), nrow(joined)))

## ---- 4. IO summary (wide-delta) for figure 5 -------------------
b <- readRDS(here("fits", "io_pooled_subset_data.rds"))
fit <- readRDS(here("fits", "io_pooled_widedelta.rds"))
fitg <- readRDS(here("fits", "io_pooled_gamma_widedelta_add.rds"))
d <- fit$draws(format = "df"); dg <- fitg$draws(format = "df")
io_summary <- list(
  obs = b$df |>
    group_by(study, ii, aa, age) |>
    summarise(prop = mean(produces), .groups = "drop"),
  params = list(
    sigma_r            = median(d$sigma_r),
    sigma_alpha        = median(d$sigma_alpha),
    sigma_zeta_base    = median(d$sigma_zeta),
    sigma_zeta_add     = median(dg$sigma_zeta),
    gamma_add          = median(dg$gamma),
    kappa_pop          = median(d$delta) + 1
  )
)
saveRDS(io_summary, file.path(CACHE, "fig5_io_summary.rds"))
cat(sprintf("Wrote %s\n", file.path(CACHE, "fig5_io_summary.rds")))

## ---- 5. Figure 4: per-word "exposures to 50% production" ---------
# For a typical English child, solve the age at which each word reaches
# 50% production from the fitted linear predictor, then count that word's
# cumulative input exposures by that age. Uses the Bayesian longitudinal
# fit `long_no_freq_slopes`, which estimates per-child variation in BOTH
# efficiency (sigma_alpha) and acceleration (sigma_zeta) — the Bayesian
# counterpart of the glmer M4 (D_log). The typical-child curve uses the
# population values (log_alpha = 0; kappa = 1 + median(delta)). CHILDES
# word frequencies via childes-db. (Mirrors model/scripts/exposure_to_learn.R.)
# Pair the recent post-dedup EN fit with its OWN bundle. long_no_freq_slopes's
# psi has 682 items (re-pulled from the GCP node, 2026-06-08 extraction), which
# matches long_subset_data.rds (the dedup'd June bundle). The earlier 671-item
# psi that triggered the count-mismatch was a STALE May extraction — refreshed.
lb   <- readRDS(here("fits", "long_subset_data.rds"))
lsd  <- lb$stan_data
# log_H: log waking hours/month (~365); a0: anchor age (mo); mu_r: typical
# log input rate (tokens/hr).
log_H <- lsd$log_H; a0 <- lsd$a0; mu_r <- lsd$mu_r
ldraws    <- as.data.frame(readRDS(here("fits", "summaries",
                                        "long_no_freq_slopes.draws.rds")))
kappa_typ <- 1 + median(ldraws$delta)           # population scaling exponent
xi_typ    <- mu_r                               # typical intercept: log_r = mu_r, log_alpha = 0
psi <- read_csv(here("fits", "summaries", "long_no_freq_slopes_psi.csv"),
                show_col_types = FALSE) |> rename(delta_j = delta_j_median)
FREQ_MIN <- 1e-5                                # drop CHILDES no-match floor items

# Lexical-class fix: the upstream scheme lumps CDI *people* (mailman, aunt,
# doctor, ...) and *places/events* (woods, camping, park, zoo, ...) into
# "other" alongside genuine non-nouns. Reclassify those concrete nouns as
# "nouns"; keep only sound effects, social routines, and time words as
# "other" (these 42 items stay):
stay_other <- c(
  # animal sounds / sound effects
  "baa baa","choo choo","cockadoodledoo","grrr","meow","moo","ouch",
  "quack quack","uh oh","vroom","woof woof","yum yum",
  # games / social routines
  "bye","call (on phone)","give me five!","go potty","gonna get you!",
  "hello","hi","night night","no","pattycake","peekaboo","please",
  "shh/shush/hush","so big!","thank you","this little piggy","turn around","yes",
  # time words
  "after","before","day","later","morning","night","now","time",
  "today","tomorrow","tonight","yesterday")

# Join word difficulties (delta_j) to words. Prefer joining by item NAME so
# the mapping survives item-set changes between the fit and the data bundle.
# Fall back to the positional jj/row-number join ONLY when the two are the
# same length; otherwise the rows misalign and delta_j gets bolted onto the
# wrong words, silently breaking Figure 4 (this happened when a May fit with
# 671 items was joined to a re-cleaned June bundle with 682).
wi <- lb$word_info |> mutate(jj = row_number())
psi_item_col <- intersect(c("item", "item_definition", "word"), names(psi))
if (length(psi_item_col) > 0) {
  exposure_items <- wi |>
    inner_join(psi |> transmute(item = .data[[psi_item_col[1]]], delta_j),
               by = "item")
} else {
  if (nrow(psi) != nrow(wi)) {
    stop(sprintf(paste0(
      "build_cache section 5: psi has %d items but word_info has %d, and psi ",
      "has no item-name column. A positional (row-number) join would MISALIGN ",
      "delta_j with words and silently break Figure 4. Re-fit the longitudinal ",
      "model on the current data bundle (so its psi matches the bundle read above), ",
      "or export psi with an item column."), nrow(psi), nrow(wi)))
  }
  exposure_items <- wi |> left_join(psi |> select(jj, delta_j), by = "jj")
}

exposure_items <- exposure_items |>
  filter(!is.na(delta_j), prob >= FREQ_MIN) |>
  mutate(lexical_class = as.character(lb$class_levels[cc]),
         lexical_class = if_else(lexical_class == "other" & !(item %in% stay_other),
                                 "nouns", lexical_class),
         lexical_class = factor(lexical_class, levels = lb$class_levels),
         t_50   = a0 * exp((delta_j - log_H - xi_typ) / kappa_typ),
         N_word = exp(xi_typ) * exp(log_H) * t_50 * prob) |>
  # Drop the residual "other" class (social routines etc.): it's a mishmash of
  # heterogeneous categories, so it's excluded from Figure 4. The reclassify
  # step above first rescues concrete nouns mislabeled "other" into "nouns".
  filter(lexical_class != "other") |>
  mutate(lexical_class = droplevels(lexical_class)) |>
  select(item, lexical_class, delta_j, prob, t_50, N_word)

# Sanity guard: in a valid fit, word difficulty (delta_j) is strongly NEGATIVELY
# correlated with the data's per-item PRODUCTION rate (easy words are produced by
# more children). It is NOT correlated with frequency: for CDI words frequency is
# ~orthogonal to difficulty (high-frequency function words like "in"/"if" are
# late-learned, cor ~ 0), so the old frequency-based guard false-alarmed on good
# fits -- the actual misalignment test is against production (see journal entry 34).
.prod <- lb$df |> dplyr::group_by(jj) |>
  dplyr::summarise(p = mean(produces), .groups = "drop") |>
  dplyr::left_join(lb$word_info |> dplyr::select(jj, item), by = "jj")
.align <- cor(exposure_items$delta_j,
              .prod$p[match(exposure_items$item, .prod$item)], use = "complete.obs")
if (is.na(.align) || .align > -0.5) {
  stop(sprintf(paste0(
    "build_cache section 5: cor(delta_j, production rate) = %.2f (expected strongly ",
    "negative). delta_j looks misaligned with words -- Figure 4 would be wrong. ",
    "Check that psi and the data bundle come from the same item set."),
    .align))
}
saveRDS(list(items = exposure_items,
             meta  = list(mu_r = mu_r, kappa_typ = kappa_typ,
                          n_items = nrow(exposure_items))),
        file.path(CACHE, "fig4_exposure.rds"))
cat(sprintf("Wrote %s (%d items, kappa=%.2f)\n",
            file.path(CACHE, "fig4_exposure.rds"),
            nrow(exposure_items), kappa_typ))

## ---- 6. Figure 6A: children-vs-LLM scaling slopes ---------------
# Port of model/scripts/chang_bergen_comparison.R. Places two populations
# on a common "logit per natural-log experience" axis:
#   LLMs  : per-(model, word) sigmoid slope = (1/ParamScale)/ln(10)
#           from Chang & Bergen (2022) fits (data/chang_bergen_2022/).
#   Kids  : per-child kappa_i = 1 + delta + zeta_i, sampled from the
#           Bayesian D posterior (long_no_freq_slopes; EN + Norwegian).
# Cached as a tidy long frame + per-population summary so the manuscript
# chunk only plots and reports numbers.
CB_DIR <- here("data", "chang_bergen_2022")
LMS <- c("bert", "bilstm", "gpt2", "lstm")
lm_slopes <- bind_rows(lapply(LMS, function(lm) {
  d <- read.delim(file.path(CB_DIR, sprintf("%s_sigmoids.txt", lm)),
                  stringsAsFactors = FALSE)
  d$model <- lm; d
})) |>
  mutate(surprisal_range = ParamUpper - ParamLower,
         slope_natural   = (1 / ParamScale) / log(10)) |>
  # drop degenerate sigmoid fits (numerical edges / never-learned words)
  filter(is.finite(slope_natural), ParamScale > 0.01, ParamScale < 10,
         surprisal_range > 1.0) |>
  mutate(label = factor(model, levels = c("bert", "gpt2", "bilstm", "lstm"),
                        labels = c("BERT", "GPT-2", "BiLSTM", "LSTM")),
         category = "LLMs (Chang & Bergen 2022)") |>
  select(label, category, slope_natural)

# Kids: per-child kappa_i BLUPs (posterior medians) from the CURRENT bayes_long
# M3 fits, 3+-administrations (_a3) variant -- the same fits behind Figs 1-2 and
# the SI BLUP analyses. English pools the three longitudinal datasets
# (thal/smith/marchman); Norwegian is its own fit. (Replaces the retired
# long_no_freq_slopes posterior simulation from the pre-bayes_long pipeline;
# real per-child estimates rather than MVN draws from population parameters.)
BL_SUMM <- here("fits", "bayes_long", "summaries")
blup_kappa <- function(slugs, label) bind_rows(lapply(slugs, function(s)
  read.csv(file.path(BL_SUMM, paste0(s, "_a3_m3_child.csv"))))) |>
  transmute(label = label, category = "Children (this work)",
            slope_natural = kappa_median)
kid_slopes <- bind_rows(
  blup_kappa(c("thal", "smith", "marchman"), "Children (English)"),
  blup_kappa("norwegian", "Children (Norwegian)"))

# ---- Our CHILDES-trained GPT-2 (this work): two experience axes ----
# Same per-word C&B slope statistic; category "LLMs (this work)".
# (a) TRAINING-step axis: per-word sigmoid over training checkpoints (3 seeds, L1).
childes_train <- bind_rows(lapply(c(42, 0, 123), function(s) {
  read.delim(here("fits", "llm", "sigmoids",
                  sprintf("gpt2_childes_seed%d_sigmoids.txt", s)),
             stringsAsFactors = FALSE)
})) |>
  mutate(surprisal_range = ParamUpper - ParamLower,
         slope_natural   = (1 / ParamScale) / log(10)) |>
  filter(is.finite(slope_natural), ParamScale > 0.01, ParamScale < 10,
         surprisal_range > 1.0) |>
  transmute(label = "GPT-2 / CHILDES (training)", category = "LLMs (this work)",
            slope_natural)

# (b) DEVELOPMENT axis: per-word sigmoid over DISTINCT-INPUT budgets, refit across
# the converged-model ladder. Update fits/llm/ladder_bestval_finer.csv when the
# final seeds land and re-run this script; the figure regenerates from it.
four_pl_scale <- function(x, y) tryCatch({
  f <- nls(y ~ lo + (up - lo) / (1 + exp((x - mid) / sc)),
           start = list(up = max(y), lo = min(y), mid = mean(x), sc = 0.5),
           lower = c(up = min(y)-5, lo = min(y)-5, mid = min(x)-3, sc = 1e-3),
           upper = c(up = max(y)+5, lo = max(y)+5, mid = max(x)+3, sc = 50),
           algorithm = "port", control = nls.control(maxiter = 200, warnOnly = TRUE))
  c(sc = unname(coef(f)["sc"]), rng = unname(coef(f)["up"] - coef(f)["lo"]))
}, error = function(e) c(sc = NA_real_, rng = NA_real_))
ladder <- readr::read_csv(here("fits", "llm", "ladder_bestval_finer.csv"),
                          show_col_types = FALSE) |>
  mutate(l10w = log10(words))
n_budgets <- max(dplyr::count(ladder, seed, word)$n)   # full-ladder rung count
childes_dev <- ladder |>
  group_by(seed, word) |> filter(dplyr::n() == n_budgets) |>
  group_modify(~{ p <- four_pl_scale(.x$l10w, .x$surprisal)
                  tibble::tibble(sc = p["sc"], rng = p["rng"]) }) |>
  ungroup() |>
  filter(is.finite(sc), sc > 0.01, sc < 10, rng > 1) |>
  transmute(label = "GPT-2 / CHILDES (development)", category = "LLMs (this work)",
            slope_natural = 0.434 / sc)

# group column for the single-panel figure: kids by language; BookCorpus pooled;
# our two CHILDES axes kept separate (training vs development).
slopes <- bind_rows(
  kid_slopes    |> mutate(label = as.character(label), group = label),
  lm_slopes     |> mutate(label = as.character(label), group = "LMs: C&B 2022 (4 architectures)"),
  childes_train |> mutate(group = "LMs: CHILDES (training)"),
  childes_dev   |> mutate(group = "LMs: CHILDES (development)")) |>
  mutate(group = factor(group, levels = c(
           "Children (English)", "Children (Norwegian)", "LMs: C&B 2022 (4 architectures)",
           "LMs: CHILDES (training)", "LMs: CHILDES (development)")),
         category = factor(category, levels = c("Children (this work)",
                            "LLMs (Chang & Bergen 2022)", "LLMs (this work)")))
slope_summary <- slopes |>
  group_by(group, category) |>
  summarise(median = median(slope_natural),
            q025   = quantile(slope_natural, 0.025),
            q975   = quantile(slope_natural, 0.975),
            n      = n(), .groups = "drop")
# Per-architecture C&B slopes, kept separately for SI: Acceleration in Other
# Architectures (the pooled C&B group is no longer shown in the main figure).
slopes_cb_arch <- slopes |>
  filter(group == "LMs: C&B 2022 (4 architectures)") |>
  transmute(arch = as.character(label), slope_natural)
saveRDS(list(slopes = slopes, summary = slope_summary,
             slopes_cb_arch = slopes_cb_arch),
        file.path(CACHE, "fig6_llm_slopes.rds"))
cat(sprintf("Wrote %s (%d slope rows)\n",
            file.path(CACHE, "fig6_llm_slopes.rds"), nrow(slopes)))

## ---- 7. Input-rate variation table (sigma_r supplement) ----------
# Per-source input rate in tokens/hr and tokens/month, plus the implied
# sigma_r (= SD of log tokens/hr across children) and the 1-SD
# multiplicative factor exp(sigma_r). Two provenance tiers:
#   * computed  : recomputed from the committed per-dyad CSV
#                 data/sperry/hourly_tokens_Sperry_HartRisley.csv, using the
#                 best channel per study (Sperry adult-to-child; S&W all
#                 speech; HR/W&F mother CDS), both whole-sample pooled and
#                 within-SES-stratum (means subtracted per band first).
#   * literature: canonical numbers from input_estimation/validation_set.csv
#                 (curated from the source PDFs).
# Tokens/month = tokens/hr x H, H = 365 waking hr/mo (12 hr/day x 30.44
# days/mo; matches MODEL_CONSTANTS$log_H). Mirrors
# model/scripts/input_rate_table.R but pulls the cross-cultural rows from
# the validation set rather than its hand-entered estimates.
H_PER_MONTH <- 365

ir <- read_csv(here("data", "sperry", "hourly_tokens_Sperry_HartRisley.csv"),
               show_col_types = FALSE) |>
  rename(mother = mother_child_tokens_hr, all_speech = all_tokens_hr,
         adult_child = adult_child_tokens_hr) |>
  mutate(best = case_when(dataset == "Sperry" ~ adult_child,
                          dataset == "Soderstrom & Wittebolle" ~ all_speech,
                          TRUE ~ mother),
         channel = case_when(dataset == "Sperry" ~ "adult-to-child",
                             dataset == "Soderstrom & Wittebolle" ~ "all speech",
                             TRUE ~ "mother CDS"))
ir_summ <- function(x) {
  x <- x[!is.na(x) & x > 0]; lx <- log(x); m <- mean(lx); s <- sd(lx)
  tibble(n = length(x), mean_hr = exp(m), lo_hr = exp(m - s),
         hi_hr = exp(m + s), sigma_r = s)
}
ir_pooled <- ir |> group_by(dataset, channel) |> reframe(ir_summ(best)) |>
  mutate(sample = "all kids pooled")
ir_within <- ir |> filter(dataset %in% c("Hart & Risley", "Sperry"),
                          !is.na(best), best > 0) |>
  mutate(lx = log(best)) |>
  group_by(dataset, channel, sample) |> mutate(dev = lx - mean(lx)) |>
  group_by(dataset, channel) |>
  reframe(n = n(), mean_hr = exp(mean(lx)), sigma_r = sd(dev),
          lo_hr = exp(mean(lx) - sd(dev)), hi_hr = exp(mean(lx) + sd(dev))) |>
  mutate(sample = "within-stratum pooled")
ir_comp <- bind_rows(ir_pooled, ir_within) |>
  transmute(source = dataset, sample, channel, n, mean_hr, lo_hr, hi_hr,
            sigma_r, kind = "computed")

# Cross-cultural / LENA rows: canonical numbers from the validation set.
vs <- read_csv(here("studies", "input_estimation", "validation_set.csv"),
               show_col_types = FALSE)
vrow <- function(src, samp = NULL) {
  r <- vs[vs$source == src, ]
  if (!is.null(samp)) r <- r[r$sample_label == samp, ]
  r[1, ]
}
seed <- vrow("Bergelson et al. 2018 (Day by day)", "SEEDLingS subset")
tsel <- vrow("Casillas Brown & Levinson 2020")
clo  <- vrow("Coffey et al. 2024 (lower bound)")
chi  <- vrow("Coffey et al. 2024 (upper bound)")
ir_lit <- tibble(
  source  = c("Bergelson 2018 (SEEDLingS, LENA)",
              "Casillas 2020 (Tseltal Mayan)",
              "Coffey et al. 2024 (plausible envelope)"),
  sample  = c("US daylong, 6--7 mo", "Chiapas daylong", "cross-cultural bounds"),
  channel = c("all adult (AWC)", "child-directed", "CDS to all-input"),
  n       = c(seed$n_kids, tsel$n_kids, NA_integer_),
  mean_hr = c(seed$tokens_per_hour_mean, tsel$tokens_per_hour_mean, NA_real_),
  lo_hr   = c(exp(seed$log_r_mean - seed$log_r_sd), NA_real_,
              clo$tokens_per_hour_mean),
  hi_hr   = c(exp(seed$log_r_mean + seed$log_r_sd), NA_real_,
              chi$tokens_per_hour_mean),
  sigma_r = c(seed$log_r_sd, NA_real_, NA_real_),
  kind    = "literature")

row_order <- c("Hart & Risley", "Soderstrom & Wittebolle", "Sperry",
               "Weisleder & Fernald", "Bergelson 2018 (SEEDLingS, LENA)",
               "Casillas 2020 (Tseltal Mayan)",
               "Coffey et al. 2024 (plausible envelope)")
input_rate_table <- bind_rows(ir_comp, ir_lit) |>
  mutate(mean_mo = mean_hr * H_PER_MONTH,
         source = factor(source, levels = row_order)) |>
  arrange(source, sample) |>
  mutate(source = as.character(source))
saveRDS(input_rate_table, file.path(CACHE, "input_rate_table.rds"))
cat(sprintf("Wrote %s (%d rows)\n",
            file.path(CACHE, "input_rate_table.rds"), nrow(input_rate_table)))

## ---- 8. io-partition scalars (population pi_alpha, EN + NO) --------
# Population input-efficiency partition reported inline in the paper section
# "Population input-related variation": pi_alpha = share of efficiency (xi)
# variance attributable to input, from the imputed-input D fits
# (long_no_freq_slopes[_norwegian], sigma_r pinned). Combined N from the fit
# bundles. NO summary is the 2026-06-11 post-dedup refit (experiments.md #35).
io_pa <- function(tag) {
  s <- as.data.frame(readRDS(here("fits", "summaries", paste0(tag, ".summary.rds"))))
  r <- s[s$variable == "pi_alpha", ]
  c(pi = r$mean, lo = r$q025, hi = r$q975)
}
io_nk <- function(bundle) {
  b <- readRDS(here("fits", bundle))
  if (!is.null(b$n_kids)) b$n_kids else length(unique(b$df$child_id))
}
en_pa <- io_pa("long_no_freq_slopes")
no_pa <- io_pa("long_no_freq_slopes_norwegian")
io_partition <- list(
  en_pi = unname(en_pa["pi"]), en_lo = unname(en_pa["lo"]), en_hi = unname(en_pa["hi"]),
  no_pi = unname(no_pa["pi"]), no_lo = unname(no_pa["lo"]), no_hi = unname(no_pa["hi"]),
  n_en  = io_nk("long_subset_data.rds"),
  n_no  = io_nk("long_subset_data_nor.rds"))
io_partition$n_total <- io_partition$n_en + io_partition$n_no
saveRDS(io_partition, file.path(CACHE, "io_partition.rds"))
cat(sprintf("Wrote %s (EN pi_alpha %.3f, NO pi_alpha %.3f, N %d)\n",
            file.path(CACHE, "io_partition.rds"),
            io_partition$en_pi, io_partition$no_pi, io_partition$n_total))

cat("\nAll caches built.\n")
