## Combine literature_estimates.csv + local_summary.csv into the
## canonical validation set.
##
## Outputs
##   input_estimation/validation_set.csv  — the deliverable
##
## Adds two derived columns relative to the input files:
##   tokens_per_month  = tokens_per_hour_mean * 365     (model H = 365)
##   reported_or_computed: "computed-from-data" for local rows,
##                        "reported" for literature rows
##
## Usage:  Rscript input_estimation/build_validation_set.R

suppressPackageStartupMessages({
  library(readr); library(dplyr)
})

ROOT <- Sys.getenv("STANDARD_MODEL_ROOT", unset = "/Users/mcfrank/Projects/standard_model_2")
OUT  <- file.path(ROOT, "input_estimation")

# ------------------------------------------------------------
# Run compute_local_rates if its outputs are missing
# ------------------------------------------------------------
local_path <- file.path(OUT, "local_summary.csv")
if (!file.exists(local_path)) {
  source(file.path(OUT, "compute_local_rates.R"), local = new.env())
}
local <- read_csv(local_path, show_col_types = FALSE, progress = FALSE)
lit   <- read_csv(file.path(OUT, "literature_estimates.csv"),
                  show_col_types = FALSE, progress = FALSE)

# ------------------------------------------------------------
# Harmonize columns
# ------------------------------------------------------------
# local_summary has: source, sample_label, language, measure_type, method,
#   citation_path, n_kids, log_r_mean, log_r_sd, tokens_per_hour_mean,
#   tokens_per_hour_median, age_range_mo
# literature_estimates has: source, sample_label, language, n_kids,
#   age_range_mo, measure_type, method, tokens_per_hour_mean,
#   tokens_per_hour_sd, log_r_mean, log_r_sd, citation_path, notes

local_norm <- local %>%
  transmute(
    source, sample_label, language, n_kids, age_range_mo,
    measure_type, method,
    tokens_per_hour_mean,
    tokens_per_hour_sd = NA_real_,    # we report log-scale sd, not raw sd
    log_r_mean, log_r_sd,
    citation_path,
    reported_or_computed = "computed-from-data",
    notes = sprintf("Per-child means computed in input_estimation/compute_local_rates.R; n=%d distinct kids.",
                    n_kids))

# Literature estimates: set reported_or_computed = "reported" by default.
# Some Sperry/HR/WF rows in the literature CSV are actually computed-from-data
# (since the data live in our pooled CSV); flag those individually.
lit_norm <- lit %>%
  mutate(
    reported_or_computed = case_when(
      grepl("^Sperry", source)             ~ "computed-from-paper-data",
      grepl("^Weisleder", source)          ~ "computed-from-paper-data",
      grepl("^Soderstrom", source)         ~ "computed-from-paper-data",
      grepl("^Bergelson .* 2018", source)  ~ "computed-from-paper-data",
      TRUE                                  ~ "reported"))

# Combine
val <- bind_rows(local_norm, lit_norm) %>%
  # Derived: tokens / month using model's H = 365 hr/month
  mutate(
    tokens_per_month = round(tokens_per_hour_mean * 365, 0),
    log_r_mean = round(log_r_mean, 4),
    log_r_sd   = round(log_r_sd, 4),
    tokens_per_hour_mean = round(tokens_per_hour_mean, 1),
    tokens_per_hour_sd   = round(tokens_per_hour_sd, 1)) %>%
  select(source, sample_label, language, n_kids, age_range_mo,
         measure_type, method,
         tokens_per_hour_mean, tokens_per_hour_sd,
         tokens_per_month,
         log_r_mean, log_r_sd,
         reported_or_computed,
         citation_path, notes) %>%
  arrange(measure_type, source, sample_label)

write_csv(val, file.path(OUT, "validation_set.csv"))
message(sprintf("[validation_set] %d rows -> input_estimation/validation_set.csv",
                nrow(val)))

# ------------------------------------------------------------
# σ_r implication report
# ------------------------------------------------------------
# What does the validation set tell us about the σ_r quantity that
# anchors RQ3 / π_α decomposition? Print four scenarios.

cat("\n=== Implications for σ_r (the model's prior on log r_i across children) ===\n")

scenarios <- list(
  "Within-Western-CDS (Sperry pooled, n=42)" =
    val %>% filter(grepl("^Sperry|POOLED", source),
                   measure_type == "CDS-any-adult",
                   reported_or_computed != "reported"),
  "Cross-site Western CDS (Sperry sites + HR groups, computed)" =
    val %>% filter(reported_or_computed == "computed-from-data" |
                   reported_or_computed == "computed-from-paper-data",
                   grepl("CDS-(any-adult|mother)", measure_type),
                   grepl("^Sperry|^Hart|^Weisleder", source)),
  "All-adult tokens (CDS+ODS) within Western samples" =
    val %>% filter(reported_or_computed != "reported",
                   measure_type == "all-adult",
                   !grepl("Spanish|Manitoba|Tseltal|Yélî|Pirahã", language)),
  "Local data only (BabyView + SEEDLingS per-child)" =
    val %>% filter(grepl("^this paper", source))
)

for (nm in names(scenarios)) {
  s <- scenarios[[nm]]
  if (nrow(s) == 0) next
  cat(sprintf("\n%s\n  rows: %d (sources: %s)\n",
              nm, nrow(s),
              paste(unique(s$source), collapse = "; ")))
  cat(sprintf("  log_r_mean range: [%.2f, %.2f]   tokens/hr: [%.0f, %.0f]\n",
              min(s$log_r_mean, na.rm = TRUE),
              max(s$log_r_mean, na.rm = TRUE),
              min(s$tokens_per_hour_mean, na.rm = TRUE),
              max(s$tokens_per_hour_mean, na.rm = TRUE)))
  finite_sds <- s$log_r_sd[is.finite(s$log_r_sd)]
  if (length(finite_sds) > 0) {
    cat(sprintf("  within-group log_r_sd: median %.3f, range [%.3f, %.3f] across %d groups\n",
                median(finite_sds),
                min(finite_sds), max(finite_sds), length(finite_sds)))
  }
}

# Bracket on σ_r via the spread of group-mean log_r values
cat("\n--- Between-group spread (σ_r implied if treating each estimate as one 'population') ---\n")

between_western <- val %>%
  filter(reported_or_computed != "reported",
         grepl("CDS-(any-adult|mother)", measure_type),
         grepl("English \\(US|English \\(CA", language)) %>%
  pull(log_r_mean)
if (length(between_western) >= 2) {
  cat(sprintf("Western CDS group means (n=%d): SD across group means = %.3f log-units\n",
              length(between_western), sd(between_western)))
}

between_cross <- val %>%
  filter(grepl("CDS", measure_type) | measure_type == "all-input") %>%
  pull(log_r_mean) %>% .[is.finite(.)]
if (length(between_cross) >= 2) {
  cat(sprintf("All CDS+all-input group means (n=%d): SD across group means = %.3f log-units\n",
              length(between_cross), sd(between_cross)))
}

cat("\nFor reference, the model currently uses σ_r = 0.534 (pooled Sperry+HR+WF n=42)\n")
cat("and the sensitivity sweep tested σ_r ∈ {0.30, 0.53, 0.80, 1.20}.\n")
cat("π_α (efficiency share of variance) ranged 0.70 (σ_r=1.20) to 0.98 (σ_r=0.30).\n")
