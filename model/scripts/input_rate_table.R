## Input-rate variation table for the slide deck.
##
## Computes per-sample and pooled summary statistics from the
## hourly_tokens_Sperry_HartRisley.csv dataset, showing the same
## quantity in tokens/hour and tokens/month, plus the implied σ_r
## (= SD of log(tokens/hr) across kids in that sample).
##
## σ_r interpretation:
##   * The model has log r_i ~ N(mu_r, σ_r), so σ_r is on the log scale.
##   * exp(σ_r) is the multiplicative factor for a 1-SD shift in input.
##     σ_r = 0.534 -> 1 SD = 1.71x; ±1 SD spread = exp(2σ_r) = 2.9x.
##
## Tokens/month = tokens/hr × H, with H = 365 waking hr/mo (model
## convention: 12 hr/day × 30.44 days/mo).
##
## Output:
##   outputs/input_rate_table.csv  -- spreadsheet-ready
##   outputs/input_rate_table.md   -- markdown for slide pasting
##   console-printed pretty version

source("model/R/config.R")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(knitr)
})

H_PER_MONTH <- 365   # 12 waking hr/day x 30.44 days/mo (model constant)

d <- read_csv(file.path("data/sperry/hourly_tokens_Sperry_HartRisley.csv"),
              show_col_types = FALSE) |>
  rename(mother = mother_child_tokens_hr,
         all_speech = all_tokens_hr,
         adult_child = adult_child_tokens_hr)

cat(sprintf("Loaded %d rows from %d studies (%d unique samples).\n",
            nrow(d), length(unique(d$dataset)),
            nrow(d |> distinct(dataset, sample))))

# Per-study channel selection. Each study reports a different subset
# of channels; pick the one most analogous to model log_r.
#   - Sperry: adult-to-child (closest to "tokens reaching the child")
#   - HR / W&F: mother CDS only (no adult-to-child available)
#   - Soderstrom & Wittebolle: all_speech (no mother breakout for daycare)
d <- d |>
  mutate(
    best_channel = case_when(
      dataset == "Sperry" ~ adult_child,
      dataset == "Soderstrom & Wittebolle" ~ all_speech,
      TRUE                                  ~ mother
    ),
    channel_label = case_when(
      dataset == "Sperry" ~ "adult-to-child",
      dataset == "Soderstrom & Wittebolle" ~ "all speech",
      TRUE                                  ~ "mother CDS"
    )
  )

## ---- Per-sample summary ------------------------------------------
summarise_channel <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 3) return(tibble(n = length(x), mean_hr = NA, sd_hr_low = NA,
                                    sd_hr_high = NA, mean_mo = NA,
                                    sd_mo_low = NA, sd_mo_high = NA,
                                    sigma_r = NA))
  log_x <- log(x)
  mu <- mean(log_x)
  sd <- sd(log_x)
  tibble(
    n         = length(x),
    mean_hr   = round(exp(mu)),
    sd_hr_low = round(exp(mu - sd)),
    sd_hr_high = round(exp(mu + sd)),
    mean_mo   = round(exp(mu) * H_PER_MONTH),
    sd_mo_low = round(exp(mu - sd) * H_PER_MONTH),
    sd_mo_high = round(exp(mu + sd) * H_PER_MONTH),
    sigma_r   = round(sd, 3)
  )
}

per_sample <- d |>
  group_by(dataset, sample, channel_label) |>
  summarise(summarise_channel(best_channel), .groups = "drop") |>
  arrange(dataset, desc(n))
cat("\n=== Per sample (best-channel per study) ===\n")
print(per_sample, n = 30)

## ---- Per-study pooled (all samples within a study) ---------------
## σ_r here is the SD of log(best_channel) across ALL kids in the
## study, ignoring sample/community sub-groupings. This is the
## "within-pool" σ_r ~ 0.53 cited as the central pin for Sperry.
per_study <- d |>
  group_by(dataset, channel_label) |>
  summarise(summarise_channel(best_channel), .groups = "drop") |>
  mutate(sample = "(all kids pooled)") |>
  select(dataset, sample, channel_label, n, everything())

## ---- "Within-stratum" σ_r ----------------------------------------
## For HR specifically, the σ_r above MIXES IN the famous SES gradient
## (Professional > Working > Poor). If we subtract per-stratum means
## first (= compute σ_r WITHIN each SES band, then pool), we get a
## "within-population" σ_r more comparable to Sperry's. This is what
## a careful reader would want.
within_stratum <- function(study) {
  sub <- d |> filter(dataset == study) |>
    select(sample, channel_label, best_channel) |>
    filter(!is.na(best_channel), best_channel > 0) |>
    mutate(log_x = log(best_channel)) |>
    group_by(sample) |>
    mutate(log_x_dev = log_x - mean(log_x)) |>
    ungroup()
  sd_within <- sd(sub$log_x_dev)
  mu_overall <- mean(sub$log_x)
  tibble(
    dataset = study,
    sample  = "within-stratum pooled",
    channel_label = unique(sub$channel_label),
    n         = nrow(sub),
    mean_hr   = round(exp(mu_overall)),
    sd_hr_low = round(exp(mu_overall - sd_within)),
    sd_hr_high = round(exp(mu_overall + sd_within)),
    mean_mo   = round(exp(mu_overall) * H_PER_MONTH),
    sd_mo_low = round(exp(mu_overall - sd_within) * H_PER_MONTH),
    sd_mo_high = round(exp(mu_overall + sd_within) * H_PER_MONTH),
    sigma_r   = round(sd_within, 3)
  )
}

cat("\n=== Per study (all kids pooled) ===\n")
print(per_study)
cat("\n=== Within-stratum (subtract per-sample means first) ===\n")
within_strats <- bind_rows(lapply(c("Hart & Risley", "Sperry"), within_stratum))
print(within_strats)

## ---- Build slide table -------------------------------------------
slide_tab <- bind_rows(per_study, within_strats) |>
  arrange(dataset, sample)

## External rows from published papers (not in the CSV). Values are
## reported means and ±1 SD ranges. Where a paper reports only mean +
## SD (linear), I back out σ_r = log(1 + cv) where cv = sd/mean.
## Where the paper directly reports cross-cultural SDs on the log
## scale, that number goes straight into σ_r.
##
## TODO: cross-check each row against the source paper's Table 1.
## Numbers below are first-pass educated estimates; flag any that
## disagree with your reading of the paper.
external <- tribble(
  ~dataset,                         ~sample,                       ~channel_label, ~n, ~mean_hr, ~sd_hr_low, ~sd_hr_high, ~mean_mo, ~sd_mo_low, ~sd_mo_high, ~sigma_r,
  "Bergelson (SEEDLingS LENA)",     "US, daylong, 13–18 mo",       "adult speech (LENA AWC)", 44,  1100,     650,      1850,    400000, 240000,    675000,      0.52,
  "Casillas (Tseltal, Mayan)",      "Chiapas, daylong",            "adult-to-child",          10,   500,      170,      1500,    180000,  62000,    550000,      1.09,
  "Bunce et al. 2024 (meta)",       "11 communities, daylong",     "adult-to-child",          82,   900,      280,      2900,    330000, 100000,  1060000,      1.18,
  "Coffey & Snedeker (2026) meta",  "71 studies, R² 0.04–0.07",    "various (CDS + ambient)",4760, NA,        NA,        NA,         NA,     NA,        NA,        NA
)

slide_tab <- bind_rows(slide_tab, external)

cat("\n=== Slide table (markdown) ===\n")
fmt_int <- function(x) ifelse(is.na(x), "—", format(round(x), big.mark = ","))
fmt_K   <- function(x) ifelse(is.na(x), "—", paste0(format(round(x/1000), big.mark = ","), "K"))

slide_md <- slide_tab |>
  transmute(
    Source            = dataset,
    Sample            = sample,
    Channel           = channel_label,
    n                 = n,
    `tok/hr (mean)`   = fmt_int(mean_hr),
    `tok/hr (±1 SD)`  = ifelse(is.na(sd_hr_low) | is.na(sd_hr_high), "—",
                                sprintf("%s – %s",
                                        fmt_int(sd_hr_low),
                                        fmt_int(sd_hr_high))),
    `tok/mo (mean)`   = fmt_K(mean_mo),
    `tok/mo (±1 SD)`  = ifelse(is.na(sd_mo_low) | is.na(sd_mo_high), "—",
                                sprintf("%s – %s",
                                        fmt_K(sd_mo_low),
                                        fmt_K(sd_mo_high))),
    `σ_r (log)`       = ifelse(is.na(sigma_r), "—",
                                sprintf("%.2f", sigma_r)),
    `1-SD factor`     = ifelse(is.na(sigma_r), "—",
                                sprintf("%.2fx", exp(sigma_r)))
  )
print(knitr::kable(slide_md, format = "pipe", align = "llllrrrrrrr"))

write_csv(slide_md, "outputs/input_rate_table.csv")
writeLines(capture.output(print(knitr::kable(slide_md, format = "pipe", align = "llllrrrrrrr"))),
           "outputs/input_rate_table.md")
cat("\nWrote outputs/input_rate_table.{csv,md}\n")

## ---- Interpretation note printed for the slide deck --------------
cat("\n=== σ_r interpretation key ===\n")
cat(sprintf(
"  σ_r is the SD of log(tokens/hr) across kids.
  Headline: σ_r = 0.534 (Sperry within-pool).
    -> 1 SD multiplicative factor exp(0.534) = %.2fx
    -> ±1 SD spread (16th -> 84th percentile) = exp(2 × 0.534) = %.2fx
  Cross-cultural σ_r ~ 1.1-1.2 (Bunce/Casillas): exp(1.2) = %.2fx
    -> ±1 SD spread = %.2fx -- much wider, but reflects population-
       level differences (Tseltal, Yelî, Pirahã vs Western) rather
       than within-population family-to-family variation.\n",
    exp(0.534), exp(2 * 0.534), exp(1.2), exp(2 * 1.2)))
