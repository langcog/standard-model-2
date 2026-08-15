# Provenance — "Children, but not language models, show accelerating returns in word learning"

**What this is.** Every claim in the paper and supplement, traced to the cache object that
carries it, the script that produced that object, and the fits underneath. It doubles as
the **manifest for the repository split**: what moves to `acceleration` is what appears
here.

**How it was built.** Not by reading code. Every build script and a full render of all
three output targets were run under an R `trace()` on `readRDS`, `readLines`, `file`,
`gzfile`, `read.csv` and `readr::read_*`, logging every path opened. The manifest below is
what the build *actually touched* — 6,511 traced opens, reduced to 319 distinct project
paths. Grep had already proved unreliable: two scripts show no `fits/` reads under a regex
and plainly read fits.

**Last verified:** 2026-08-15, against `origin/master`.

---

## 1. External inputs — the complete set

The entire manuscript build, from raw fits to rendered PDF/docx, needs exactly this:

| input | files | what it is |
|---|--:|---|
| `fits/bayes_long/` | 234 | Stan bundles, summaries, draws, diagnostics, psi/child exports |
| `fits/llm/ladder_bestval_finer.csv` | 1 | CHILDES developmental ladder, 10 seeds × 18 budgets |
| `fits/llm/register_bestval.csv` | 1 | BabyLM / ClimbMix composition control |
| `fits/english_word_freq.rds` | 1 | CHILDES unigram frequencies for the exposure axis |
| `data/harmonized/input_level.csv` | 1 | per-recording input rates (BabyView, SEEDLingS, AM2018, FMW2013) |
| `paper/cache/*.rds` | 21 | committed intermediates; the manuscript renders from these alone |

**Nothing else is touched.** Not `fits/glmer_ladder`, `fits/io_anchored`, `fits/summaries`,
`fits/io_pooled*`, `studies/io_proc_glmer`, `studies/cross_sectional_demographics`,
`studies/input_estimation`, `model/`, or `data/` beyond the one CSV. That is a much smaller
footprint than the repository's 43 GB of fits suggests, and it is measured rather than
assumed.

All of `fits/bayes_long` and `fits/llm` is archived publicly at
**<https://redivis.com/datasets/datapages.acceleration:a1c7>** (v1.0, 314 files, verified
by md5 after upload).

## 2. Claim → display item → cache → script

### Main text

| claim | item | cache | produced by |
|---|---|---|---|
| The accumulator family, and what each rung predicts | Fig 1A | `fig1_fan$conceptual` | `build_cache_short.R` |
| M3 tracks five samples in three languages; M0/M1 cannot | Fig 1B | `fig1_fan` (`fan`, `ablations`, `spag_all`) | `build_cache_short.R` |
| Exposures needed per word fall with age; M0 predicts almost no decline | Fig 2 | `fig2_efficiency`, `si_fig2_m0` | `build_cache_short.R`, `04_fig2_m0_counterfactual.R` |
| κ ≈ 11–13 for children, ≈ 1 for LMs; children vary, LMs do not | Fig 3 | `fig6_llm_slopes` (`slopes`, `vary`, `vary_cv`, `fig3`) | `build_cache_short.R` |
| Sample sizes, ages, N per dataset | inline | `si_inline`, `bayes_long_sample` | `build_cache_short.R`, `qc_exclusion_report.R` |

### Supplement

| claim | item | cache | produced by |
|---|---|---|---|
| Cross-sectional data cannot separate level from slope | fig. S1 | *(none — self-contained simulation)* | the chunk itself |
| Dataset inventory | table S1 | `si_datasets` | `build_cache_short.R` |
| Which children the QC filter removes | fig. S2 | `qc_spaghetti_data` | `qc_exclusion_report.R` |
| Full LOO ladder, both admin thresholds | table S2 | `si_loo` | `build_cache_short.R` |
| Per-child acceleration does not forecast | table S3 | `si_cv_depth` | `05_cv_depth_ladder.R` |
| Acceleration itself does forecast (M2 vs M2₀) | table S4 | `si_cv_headline` | `08_cv_headline_report.R` |
| M3 matches a non-parametric GAMLSS fit | fig. S3 | `si_gamlss` | `build_cache_si_gamlss.R` |
| κ robust across 2+/3+/pooled settings | table S5 | `si_settings` | `build_cache_si_settings.R` |
| κ robust to the exclusion filter | table S6 | `si_qc_sensitivity` | `06_qc_sensitivity_report.R` |
| Log-age beats linear age | table S7 | `si_loglin` | `build_cache_short.R` |
| 2PL fits better; κ and AoA barely move | table S8 | `si_2pl` | `03_compare_1pl_2pl.R` |
| κ is not an artifact of item selection | fig. S4 | `si_item_subset` | `11_item_subset_report.R` |
| Per-child κ distributions | fig. S5, table S9 | `si_blups`, `si_ranef` | `build_cache_short.R` |
| Input rate is flat across development | fig. S6 | `si_input_stability` | `07_input_stability.R` |
| M0's implied AoA range is impossible | fig. S7, table S10 | `si_fig2_m0` | `04_fig2_m0_counterfactual.R` |
| κ is not an artifact of the estimator | fig. S8, table S11 | `si_perword_4pl` | `10_perword_4pl.R` |
| The LM result is not GPT-2-specific | fig. S9 | `fig6_llm_slopes$slopes_cb_arch` | `build_cache_short.R` |

## 3. Traps — things that will bite whoever rebuilds this

**`qc_exclusion_report.R` cannot run offline.** It pulls Wordbank live and the endpoint it
uses (`https://langcog.github.io/wordbank-datapage//db_args`) currently returns 404, so it
fails outright. Its two outputs — `bayes_long_sample.rds` and `qc_spaghetti_data.rds` — are
committed, so the manuscript builds; but the sample sizes quoted in the main text and
fig. S2 cannot be re-derived from scratch until that endpoint is fixed or the script gains
a mirror fallback. `build_cache_short.R` already has one (`WB_MIRROR` →
`paper/cache/wordbank_items_en_ws.rds`) and runs offline cleanly.

**`10_perword_4pl.R` silently truncates when run bare.** Its default argument is
`c("thal_a3", "marchman_a3")` — two of five datasets. Running `Rscript
studies/bayes_long/10_perword_4pl.R` with no arguments overwrites the five-dataset cache
with a two-dataset one, and nothing errors. It must be invoked with all five slugs. (This
was found the hard way while tracing.)

**Three caches have more than one producer.** `fig2_efficiency` is written by both
`build_cache_short.R` and `04_fig2_m0_counterfactual.R`; `fig6_llm_slopes` by both
`build_cache_short.R` and the input paper's `build_cache.R`; `si_cv_headline` and
`bayes_long_sample` are read in `_setup_shared.R`. Run order matters, and the input
paper's `build_cache.R` must not be run after the split expecting the same file.

**`fig. S1` depends on nothing.** It is a self-contained simulation with a fixed seed, so
it needs no fit and no data — worth knowing when checking whether the archive is complete.

## 4. What stays behind

The input paper keeps `paper/input_paper/`, `build_cache.R`, `build_input_cache.R`,
`studies/io_proc_glmer`, `studies/cross_sectional_demographics`, `studies/input_estimation`,
`studies/proc_dp`, `studies/io_pooled`, and the corresponding fits. Its `build_cache.R`
reads `fits/bayes_long` and `fits/llm` today; after the split it must either stop doing so
or vendor what it needs, since those directories move.

`journal/PROVENANCE.md` (data sources) and `journal/paper_models_provenance.md` (last
audited 2026-06-13) describe the **input paper** and stay with it. This file supersedes
them for anything concerning the acceleration paper.
