# input_estimation

Validation set and supporting analyses for the **per-child input-rate
prior** that drives RQ3 / π_α in the standard-model project.

The full model (`notes/model_explainer.pdf`) has a single identifiability
lever for the variance-decomposition question (RQ3): the population SD
of `log r_i` across children, called **σ_r**. We pin σ_r externally from
data and check the conclusion under a sensitivity sweep. This directory
collects the data and citations that justify the chosen pin and bracket
the sensitivity range.

## What's here

```
input_estimation/
├── README.md                      ← this file
├── compute_local_rates.R          ← reads BabyView + SEEDLingS + Sperry CSV,
│                                    writes the three local_*.csv files below
├── literature_estimates.csv       ← hand-curated rows from PDFs in
│                                    papers/input_estimation/
├── build_validation_set.R         ← combines local + literature and prints
│                                    a σ_r implication report
├── local_per_recording.csv        ← (generated) one row per recording/dyad
├── local_per_child.csv            ← (generated) one row per child after
│                                    aggregating recordings
├── local_summary.csv              ← (generated) one row per (source, sample,
│                                    measure_type) group
└── validation_set.csv             ← (generated) the canonical output:
                                     local + literature merged, with units
                                     standardised and σ_r-relevant columns
```

## Units

Every estimate is reported in **adult tokens per hour of waking
observation**, matching the model's `r_i` parameter. Two parallel
columns:

| column | what | scale |
|---|---|---|
| `tokens_per_hour_mean` | geometric mean tokens/hr | natural |
| `log_r_mean` | natural log of `tokens_per_hour_mean` | log |
| `tokens_per_hour_sd` | (rare) raw SD across children, when reported | natural |
| `log_r_sd` | SD of `log r_i` across children — **this is σ_r** | log |
| `tokens_per_month` | `tokens_per_hour_mean × 365` (model H = 365 hr/mo) | natural |

`H = 365 hr/mo` is the model's convention (12 waking hours/day × 30.44
days/month). It matches `MODEL_CONSTANTS$log_H` in `model/R/config.R`.

## Measure types (controlled vocabulary)

Distinguishes child-directed from overheard speech where the source
reports the breakdown.

| `measure_type` | meaning |
|---|---|
| `CDS-mother` | tokens spoken to the child by mother only (Hart & Risley 1995, some Sperry analyses) |
| `CDS-any-adult` | tokens spoken to the child by any adult (Sperry "adult_child_tokens_hr"; Bergelson 2019) |
| `CDS-target-child-directed` | TCDS in the cross-cultural literature (Bunce et al. 2024) |
| `CDS-child-directed-vocalizations` | non-token-counted vocalizations directed to child (Cristia 2023 review) |
| `CDS-optimal` | "ideal" CDS under maximum filtering (Coffey 2024 lower-bound row) |
| `all-adult` | every adult token in the recording, CDS + ODS lumped (LENA AWC, Sperry "all_tokens_hr") |
| `ODS-only` | tokens not directed to the child (Kachergis 2022 standard-model row) |
| `ADS` | adult-directed speech overheard by the child (Bergelson 2019) |
| `all-input` | total speech input under the most liberal Coffey 2024 bound (CDS+ODS+ADS+overlap) |
| `all-input-with-overlap` | above, with overlapping speakers counted multiply (Yélî Dnye max) |

When a paper reports min/hr of speech rather than tokens/hr (e.g.,
Casillas, Bunce), we apply Coffey 2024's 0.3-s-per-word conversion
(200 words/min, from the Brent corpus CDS rate) and flag the
derivation in `notes`.

## reported_or_computed

| value | meaning |
|---|---|
| `computed-from-data` | We have per-recording rows in `local_per_recording.csv` |
| `computed-from-paper-data` | The paper's data is in `data/raw_data/sperry/hourly_tokens_Sperry_HartRisley.csv` (per-dyad) and our summary aggregates it |
| `reported` | The number is taken from the paper's text or table directly; we did not re-derive |

## Citations and PDFs

PDFs are in `papers/input_estimation/`:

| paper | file | role |
|---|---|---|
| Coffey, Räsänen, Scaff & Cristia 2024 (Interspeech) | `coffey_etal_2024.pdf` | The headline integration: derives 1–3,300 hr/year envelope; cites every other key source. |
| Sperry, Sperry & Miller 2019 (Child Dev) | `sperry_etal_2019.pdf` | The five-site re-examination of HR; per-dyad data feeds the model's pooled prior. |
| Weisleder & Fernald 2013 (Psych Sci) | `weisleder_fernald_2013.pdf` | Spanish-speaking low-SES, 29 children at 19 mo; CDS vs ODS distinction. |
| Bergelson, Soderstrom et al. 2023 (PNAS) | `bergelson_etal_2023.pdf` | 1,001 children, 6 continents — largest cross-cultural LENA dataset. SES not significant predictor of child speech. |
| Bergelson, Casillas et al. 2019 (Devel Sci) | `bergelson_etal_2019.pdf` | 61 NA children, 4 cities — establishes ~10.8 min/hr CDS as North American baseline. |
| Bergelson, Amatuni et al. 2018 (Devel Sci) | `bergelson_etal_2018.pdf` | The SEEDLingS 6/7 mo paper. Shows hour-long video vs. day-long audio give different pictures. |
| Casillas, Brown & Levinson 2020 (Child Dev) | `casillas_etal_2020_tseltal.pdf` | Tseltal Mayan, 10 kids: 3.6 min/hr CDS — about ⅓ of NA estimate. |
| Casillas, Brown & Levinson 2021 (J Child Lang) | `casillas_etal_2021_yeli.pdf` | Yélî Dnye, 19 kids in PNG: similar to Tseltal, low CDS rates. |
| Cychosz et al. 2021 (Devel Sci) | `cychosz_etal_2021.pdf` | 49 children, 5 language/cultural contexts; canonical babble across cultures (vocal output, not adult talk). |
| Dailey & Bergelson 2022 (Devel Sci) | `dailey_bergelson_2022.pdf` | Quantitative meta-analysis of SES × input quantity; 19 studies, 1,991 kids. |
| Räsänen et al. 2019 (Speech Comm) | `rasanen_etal_2019.pdf` | LENA AWC validation across 6 language environments. |
| Bunce et al. 2024 (preprint) | `bunce_etal_2024.pdf` | 5-culture TCDS+ADS comparison: NA/UK English, Argentinian Spanish, Tseltal, Yélî Dnye. |
| Soderstrom et al. 2021 (Collabra) | `soderstrom_etal_2021.pdf` | ACLEW cross-cultural annotation methodology. |
| Rowe 2012 (Child Dev) | `rowe_2012.pdf` | Quantity vs. quality of CDS in 50 dyads at 18/30/42 mo. |
| Dupoux 2018 (Cognition) | `dupoux_2018.pdf` | The reverse-engineering roadmap that motivates needing input-quantity bounds. |

**Cited but PDF unavailable** (paywalled when we tried, or only abstract):

- Cristia 2023 (Devel Sci e13265): the systematic review of infant-directed
  vocalization across 29 reports. Numbers extracted from Coffey 2024's
  reading and Bunce 2024's citations.
- Loukatou, Scaff, Demuth, Cristia, Havron 2022 (J Child Lang):
  French vs. Sesotho input split.
- Hart & Risley 1995 (book): the original 30-million-word-gap
  monograph. Numbers extracted from Coffey 2024 and Sperry 2019 reading.
- Gilkerson et al. 2017 (AJSLP): LENA Foundation normative sample.
  Adult-words-per-12-hr-recording extracted from secondary citations.

## What this validates

The model's external prior is currently:

```
ξ_i = log r_i + log α_i ~ N(μ_r, σ_r² + σ_α²)
μ_r = log 1198 ≈ 7.09           # model_explainer.pdf §4.3
σ_r = 0.534                     # = sd(log adult_child_tokens_hr) over n=42 Sperry+ rows
```

(Note: `μ_r = log 1198` in `model_explainer.pdf` rounds the Sperry-pool
geometric mean to `exp(7.09) ≈ 1198`, but the actual pooled mean of
`log r` is 7.34 → `exp(7.34) ≈ 1537` — a ~28% discrepancy worth
documenting in the explainer in passing.)

`build_validation_set.R` prints four σ_r scenarios. Headline numbers:

| scenario | n groups | within-group σ_r (median) | between-group σ_r |
|---|---:|---:|---:|
| Within-Western-CDS (Sperry pooled) | 7 | **0.534** | — |
| Cross-site Western CDS (Sperry sites + HR + WF) | 17 | 0.571 | 0.349 |
| All-adult tokens, Western only | 9 | 0.425 | — |
| All CDS+all-input groups (incl. Tseltal/Yélî/Pirahã) | 27 | — | **1.237** |

**Implication for the model's σ_r sweep** (which tested
σ_r ∈ {0.30, 0.53, 0.80, 1.20}):

- **σ_r = 0.30** is below any defensible single-site estimate (every
  Sperry/HR/WF subgroup has within-group SD ≥ 0.35). Treat as a lower
  *robustness* check, not a plausible value.
- **σ_r = 0.534** matches Sperry within-pool SD exactly. Defensible if
  the Wordbank sample is well-modeled as Sperry-like.
- **σ_r ≈ 0.80** is in the right range if Wordbank includes families
  spanning the welfare-to-professional Hart & Risley range (HR
  Welfare-group within-SD = 0.52; cross-SES SD = 0.52 per H&R, larger if
  also including Sperry's racial/site variation, around 0.6–0.9).
- **σ_r = 1.20** requires invoking non-WEIRD variation (Tseltal/Yélî
  Dnye CDS rates push the cross-cultural log-r SD above 1.0). The
  Wordbank CDI:WS sample is overwhelmingly WEIRD, so this is the
  *upper bracket*, not a plausible point estimate.

So the model's headline conclusion ("input explains ≤ 30% of
between-child variance even at σ_r = 1.2; ≤ 6% at σ_r = 0.3") is
**robust to σ_r over the range supported by Western-input data**, and
even the cross-cultural upper bracket falls within the existing
sensitivity sweep.

## Caveats

1. **Tokens vs. types.** Everything here is tokens. Type-level
   estimates (which would be much smaller) are not in this validation
   set.
2. **Speech rate conversion.** Min/hr → tokens/hr uses the Brent-corpus
   CDS rate of 0.3 s/word (≈ 200 words/min). Coffey 2024 cites this
   from Brent & Siskind 2001, but adult-directed speech rates can be
   different; tokens/hr from min/hr conversions in this CSV should be
   treated as ±20% accurate.
3. **WEIRD over-representation.** Most well-quantified data is from
   English (and a bit of Spanish) in North America. The Tseltal,
   Yélî Dnye, Pirahã, Tsimane', and Sesotho rows are individually
   small-n and methodologically heterogeneous.
4. **Maternal vs. all-adult vs. overheard.** Be careful comparing
   estimates across `measure_type` — the same population can yield a
   3× different rate depending on which speakers are counted.
5. **Reactivity / Hawthorne effects.** Some studies (Hart & Risley,
   Sperry, Casillas) involve experimenter visits, others (LENA) use
   passive recordings. Bergelson 2018 ("Day by day") shows hour-long
   video and daylong audio give different (2-4× different) pictures
   of the same children. The BabyView prior `β_react ~ N(0.4, 0.4)`
   in `prepare_babyview.R` is meant to absorb this.
6. **Age range.** Most useful estimates are for 6–24 mo. Older-child
   data (Sperry 18-48 mo, Hart & Risley 7-36 mo) is sparser.

## How to use

To regenerate everything (after editing `literature_estimates.csv` or
data):

```bash
Rscript input_estimation/compute_local_rates.R
Rscript input_estimation/build_validation_set.R
```

To use the validation set in a model fit / sensitivity check, read
`validation_set.csv` and filter to the rows whose `measure_type` and
`language` match your target population. The pooled
`POOLED (model's external prior on log r_i)` row is what
`load_input_rate_prior()` in `model/R/helpers.R` currently uses.
