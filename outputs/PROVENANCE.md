# Slide deck provenance map

Maps each load-bearing asset in `outputs/slides/standard_model_pptx.pdf`
to the script that produces it, the fits it depends on, and the bundle
data it consumes.

**Anchor**: revised PowerPoint deck rendered to PDF at
`outputs/slides/standard_model_pptx.pdf` (May 23, 2026; 46 slides). Mike
edits the `.pptx` directly; the `.qmd` source is retired.

**Conventions:**

- *Script*: file that, when run, produces the asset.
- *Fit*: posterior summaries in `fits/summaries/<tag>.{summary.rds,
  draws.rds}` and `<tag>_psi.csv` (per-item delta_j medians).
- *Bundle*: `fits/<dataset>_subset_data*.rds` containing the Stan input
  and the raw item-level data.
- *Static schematic*: figure with no model dependency (illustrations).

---

## Slides 1–9 (intro)

| Slide | Asset | Source | Status |
|---|---|---|---|
| 1 | Title text | (no figure) | — |
| 2 | AI disclosure | (no figure) | — |
| 3 | "Vocabulary growth" | (no figure) | — |
| 4 | `precursor_acceleration.png` | Pre-existing Wordbank empirical fan; **byte-match confirmed in PPTX media (image2.png)** | current |
| 5 | `precursor_acceleration.png` (same figure reused) | — | current |
| 6 | `kachergis_buckets.png` | Static schematic from Kachergis et al. 2021 | static |
| 7 | `irt.png` | Static schematic (converted from `irt.pdf` by `outputs/slides/build_pptx.py` for PPT compatibility) | static |
| 8 | Coffey & Snedeker forest plot | External — copied from C&S 2026 manuscript | static |
| 9 | "This work" outline | (no figure) | — |

## Slides 10–15 (building the model)

The 4-step build: `pure → +α → +κ_pop → +ζ_i`. **No s_i step** (per §23 of `experiments.md`).

| Slide | Asset | Source | Status |
|---|---|---|---|
| 10 | "Building the model" section divider | — | — |
| 11 | `kachergis_buckets.png` | Static | static |
| 12 | "Rasch IRT" — equation only | — | — |
| 13 | "Rasch IRT + accumulator dynamics" — equation | — | — |
| 14 | "Rasch + per-child efficiency" — equation, defines π_α | — | — |
| 15 | "Rasch + efficiency slope" — equation, defines κ_pop | — | — |

## Slides 16–23 (longitudinal results)

| Slide | Asset | Source | Status |
|---|---|---|---|
| 16 | "+ per-child efficiency slope" — equation | — | — |
| 17 | "Results: EN + NO longitudinal" divider | — | — |
| 18 | "Data and compute" — text | — | — |
| 19 | `quantile_demo_4panel_I200.png` | Script: [`model/scripts/quantile_demo_4panel_I200.R`](../model/scripts/quantile_demo_4panel_I200.R). Fits: `long_demo_pure`, `long_demo_alpha`, `long_demo_kappa`, `long_no_freq_slopes_english_I200`. Bundle: `long_subset_data_I200.rds`. | current |
| 20 | `theta_spaghetti_4panel_I200.png` | Script: [`model/scripts/theta_spaghetti_4panel_I200.R`](../model/scripts/theta_spaghetti_4panel_I200.R). Same fits + bundle as slide 19. | current |
| 21 | `m_best_quantile_EN_NO_wordbank.png` | Script: [`model/scripts/quantile_demo_mbest_EN_NO_wordbank.R`](../model/scripts/quantile_demo_mbest_EN_NO_wordbank.R). Fits: `long_no_freq_slopes` (EN), `long_no_freq_slopes_norwegian` (NO). Bundles: `long_subset_data.rds`, `long_subset_data_nor.rds`. Empirical: `long_items.rds` (wordbank-wide) via `build_xsec_empirical_wordbank`. | current |
| 22 | `exposure_to_learn_EN.png` | Script: [`model/scripts/exposure_to_learn.R`](../model/scripts/exposure_to_learn.R). Fit: `long_no_freq_slopes` (EN). Bundle: `long_subset_data.rds` (uses `word_info$prob` for per-word CHILDES freq). | current |
| 23 | "Best-fitting model: parameter posteriors" table (EN, NO M_best columns) | Script: [`model/scripts/param_table.R`](../model/scripts/param_table.R). Outputs: `outputs/param_table.{csv,md,xlsx}`. Pulls from `fits/summaries/long_no_freq_slopes{,_norwegian}.summary.rds`. | current |
| 24 | "Cross-language acceleration" text | Numbers from §25 of `experiments.md`; same fits as slide 23. | current |

## Slides 25–28 (input quantity)

| Slide | Asset | Source | Status |
|---|---|---|---|
| 25 | "Input quantity" divider | — | — |
| 26 | "Where does σ_r come from?" — text + formula | (no figure) | — |
| 27 | "Empirical estimates of input-rate variation" — text/table | Sources: Sperry et al. 2019, Hart & Risley 1995, Weisleder & Fernald 2013, etc. | external |
| 28 | `sigma_r_analytical_sensitivity.png` | Pre-existing; **byte-match confirmed in PPTX media (image37.png)**. Script provenance: TBD (likely a one-off analysis script from §18); content currently consistent. | current |
| 29 | "Caveat: input quantity, not quality" — text | (no figure) | — |

## Slides 30–34 (IO + Peekbank extensions)

| Slide | Asset | Source | Status |
|---|---|---|---|
| 30 | "Extensions: directly observed input and processing" divider | — | — |
| 31 | "Model extensions: testing π_α with richer data" — text | — | — |
| 32 | `babyview_io_input.png` + `seedlings_io_input.png` (side-by-side) | Pre-existing; **byte-match confirmed (image43, image46)**. Numbers in caption (π_α = 0.84 / 0.94) from `io_no_freq_slopes` and `io_no_freq_slopes_seedlings` summaries. | current |
| 33 | `m_best_input_quartile_io.png` ("Minimal input effects") | Script: [`model/scripts/quantile_demo_io_input.R`](../model/scripts/quantile_demo_io_input.R). Fits: `io_no_freq_slopes` + `io_no_freq_slopes_seedlings`. Bundles: `babyview_subset_data.rds`, `seedlings_subset_data.rds`. | current |
| 34 | "Processing speed readout: Peekbank LWL" — text only | Numbers from `long_proc_no_freq_slopes.draws_full.rds` (γ_rt, π_α, ρ(ζ, rtslope)). | current |
| 35 | `m_best_rt_quartile_proc.png` ("Peekbank RT results") | Script: [`model/scripts/quantile_demo_proc_rt.R`](../model/scripts/quantile_demo_proc_rt.R). Fit: `long_proc_no_freq_slopes`. Bundle: `stanford_linked_subset_data.rds`. Uses `lwl_log_rt`, `lwl_log_age` from stan_data. | current |
| 36 | "π_α across five samples" table | Script: [`model/scripts/param_table.R`](../model/scripts/param_table.R) (same as slide 23). 5-fit table. | current |

## Slides 37–41 (LLMs)

| Slide | Asset | Source | Status |
|---|---|---|---|
| 37 | "Connecting to LLMs" divider | — | — |
| 38 | "Why look at LLMs?" — text | — | — |
| 39 | "Per-word sigmoids in LMs" — Chang & Bergen 2022 formula | External | static |
| 40 | "Could the gap be an input-distribution artifact?" — Feng et al. 2026 setup | External | static |
| 41 | `feng_chang_bergen_slope_comparison.png` ("CDS training does not move LM slope") | Pre-existing; **byte-match confirmed (image52.png)**. Generated by LM-side pipeline (separate from this repo's R/Stan scripts). | current |
| 42 | `D1_scaling_disanalogy.png` ("Aggregate view: scaling laws") | Pre-existing; **byte-match confirmed (image53.png)**. Same LM-side pipeline. | current |

## Slides 43–46 (synthesis + appendix)

| Slide | Asset | Source | Status |
|---|---|---|---|
| 43 | "Synthesis" divider | — | — |
| 44 | "Technical summary" — text bullets | Numbers from slide 23 table | current |
| 45 | "Two signatures, twin constraints on theory" — `precursor_acceleration.png` + `precursor_variability.png` | Pre-existing; **both byte-match (image2.png, image55.png)** | current |
| 46 | `m_best_spaghetti.png` (appendix: per-child trajectories) | Pre-existing; **byte-match confirmed (image56.png)**. Generated from M_best posterior (script provenance to confirm — likely [`model/scripts/m_best_spaghetti.R`](../model/scripts/m_best_spaghetti.R) if it exists, else the §17 era script). | current |

---

## Numerical claims cited inline

Each is sourced from the same `fits/summaries/<tag>.summary.rds` files that
back the headline table:

| Slide | Claim | Source |
|---|---|---|
| 23 | EN: δ=10.31, σ_α=1.81, σ_ζ=3.81, π_α=0.92, ρ=−0.09 | `long_no_freq_slopes.summary.rds` |
| 23 | NO: δ=11.47, σ_α=2.05, σ_ζ=4.79, π_α=0.94, ρ=−0.13 | `long_no_freq_slopes_norwegian.summary.rds` |
| 24 | κ_pop ≈ 10–12 in both languages (= 1+δ); EN admins/kid ≈ 3, NO ≈ 8 | Same |
| 28 | π_α range 0.91–0.94 across σ_r sweep | §18 sensitivity analysis |
| 32 | BabyView π_α = 0.84, SEEDLingS π_α = 0.94 | `io_no_freq_slopes{,_seedlings}.summary.rds` |
| 34 | γ_rt = 0.084 [0.044, 0.125]; π_α = 0.90 [0.86, 0.94]; ρ(ζ, rtslope) = +0.05 [−0.26, +0.34] | `long_proc_no_freq_slopes.draws_full.rds` |
| 36 | π_α table for 5 samples | `param_table.csv` |
| 42 | "1 + δ ≈ 10.4" in the children panel | `long_no_freq_slopes.summary.rds`. Note: slide says 10.4; current posterior median is 10.31. Difference is the rounding choice; update if presenting. |

---

## Shared infrastructure

Used by multiple plot scripts:

| File | Purpose |
|---|---|
| [`model/R/config.R`](../model/R/config.R) | `PATHS`, `DEFAULT_PRIORS`, `MODEL_CONSTANTS`. Defines `s_prior_mean=0, s_prior_sd=0.001` (the post-excision default). |
| [`model/R/helpers.R`](../model/R/helpers.R) | `variant_hyperpriors`, `fit_variant_cmdstanr` (with the `output_dir` patch and `save_object` `try()` wrap), `summarize_fit`. |
| [`model/R/datasets.R`](../model/R/datasets.R) | Dataset registry: `english`, `english_I200`, `norwegian`, `babyview`, `seedlings`, `stanford_linked`. |
| [`model/R/empirical_xsec_helper.R`](../model/R/empirical_xsec_helper.R) | `build_xsec_empirical` (bundle-internal) + `build_xsec_empirical_wordbank` (wordbank-wide) + `fit_xsec_quantile_fan` (GAMLSS BEINF). |
| [`sherlock/recover_from_csvs.R`](../sherlock/recover_from_csvs.R) | Scalar + delta_j recovery from persistent cmdstanr CSV output. Used when `save_object` crashes (disk-full). |
| [`model/stan/log_irt_long.stan`](../model/stan/log_irt_long.stan) | Cleaned Stan model: `η = ξ_i + log H + (1+δ+ζ_i)·log(age/a_0) − δ_j`. No β_c, no λ_j, flat δ_j hyperprior. |
| [`model/stan/log_irt_long_io.stan`](../model/stan/log_irt_long_io.stan) | IO extension: adds per-kid observed log_r channel. |
| [`model/stan/log_irt_long_io_comp.stan`](../model/stan/log_irt_long_io_comp.stan) | IO + comprehension channel. |
| [`model/stan/log_irt_long_proc.stan`](../model/stan/log_irt_long_proc.stan) | + Peekbank LWL RT channel (γ_rt link). |

## Reproducing a fit from scratch

```sh
# Bundle prep (deterministic, seeded)
Rscript model/scripts/prepare_longitudinal_data.R "English (American)" 500 671
Rscript model/scripts/prepare_longitudinal_norwegian.R 500 680

# Fit
Rscript model/scripts/fit_longitudinal.R long_no_freq_slopes english
Rscript model/scripts/fit_longitudinal.R long_no_freq_slopes norwegian

# Extract scalars + delta_j medians
Rscript sherlock/extract_summary_table_only.R long_no_freq_slopes
Rscript sherlock/extract_scalar_draws.R long_no_freq_slopes
Rscript sherlock/extract_delta_j_slim.R long_no_freq_slopes

# Render plots
Rscript model/scripts/quantile_demo_mbest_EN_NO_wordbank.R
Rscript model/scripts/exposure_to_learn.R
Rscript model/scripts/param_table.R
```

For IO and Peekbank fits, use [`model/scripts/fit_io.R`](../model/scripts/fit_io.R) and
[`model/scripts/fit_stanford_linked.R`](../model/scripts/fit_stanford_linked.R)
respectively. The GCP batch driver in
[`gcp/run_fit.sh`](../gcp/run_fit.sh) wraps the fit + extract pipeline
into a single invocation.
