# Experiments log

A running record of fits, findings, and backlog for the standard-model
project. For the durable model specification, see
[`model_explainer.pdf`](model_explainer.pdf). The shared
**child- and item-sampling strategy** (used identically across Wordbank
longitudinal, BabyView, and future datasets) is documented in §
"Sampling strategy: children and items" of the explainer.

---

## Status key

- 🟢 completed
- 🟡 running / active
- ⚪ queued / backlog

---

## 🟢 1. Parameter recovery on synthetic data

**Setup:** `simulate_data()` → 250 children × 150 words × 3 classes,
*N* = 37 500 obs. True σ_α = 0.5, μ_c = (6.5, 8.0, 9.5), τ_c = (0.5, 0.7,
0.7), s = 4.5, δ = 0.1. Fit with 4 chains × 1500 iter (750 warm-up).
Runtime ~4 min.

**Result:** every scalar parameter recovered within its 95 % credible
interval. Word-level ψ<sub>j</sub> recovered with correlation 0.955 and
exactly nominal CrI coverage. 1 divergence out of 3 000. `π_α` recovered
as 0.624 [0.50, 0.72] against true 0.61.

**Key finding:** the collapsed-ξ parameterization (`log_irt.stan`) is
well-identified and samples cleanly on synthetic data.

**Artifacts:** `model/fits/recovery.rds`, `model/figs/recovery_*.png`.

---

## 🟢 2. Cross-sectional subset fits (English WS, 500 × 200)

**Setup:** Stratified subsample of 500 children × 200 items from
Wordbank English CDI:WS; *N* = 100 000. External prior on input rate
from `adult_child_tokens_hr` in
`hourly_tokens_Sperry_HartRisley.csv` (*n* = 42 samples): μ_r = 7.34,
σ_r = 0.53. All fits use 4 chains × 1500 iter.

### 2a. Baseline (Rasch, s / δ free)

| | posterior | CrI | Rhat |
|---|---:|---|---|
| σ_α | 1.98 | [1.85, 2.12] | 1.06 |
| π_α | 0.93 | [0.92, 0.94] | 1.06 |
| s | 12.2 | [10.3, 13.8] | 1.02 |
| δ | 2.31 | [1.65, 2.98] | 1.01 |

**Key finding:** model lands in an extreme corner — pure linear
accumulation cannot produce the per-word growth-curve steepness the
data show (∼12× too slow per word under Rasch). The model cranks δ, s,
and σ_α to compensate. RQ1 answer: r(ψ<sub>j</sub>, log p<sub>j</sub>) = 0.68 — frequency
positively correlated with threshold, a Braginsky-2019 style result.

### 2b. Diagnostic 2 × 2 (fix_delta / fix_s / both_fixed)

Same data, but with (δ, s) each fixed or free:

| variant | σ_α | s | δ | π_α |
|---|---:|---:|---:|---:|
| baseline | 1.98 | 12.2 | 2.31 | 0.93 |
| fix_delta | 2.20 | **15.0** (prior bound) | 0 | 0.94 |
| fix_s | 2.06 | 2.00 | 4.05 | 0.94 |
| both_fixed | — | — | — | chain stuck, killed |

**Key finding:** σ_α robust at ~2 across every variant. s and δ
trade off enormously (fix one, the other compensates). The extreme
"93% efficiency" number is *not* an artifact of the s/δ trade-off.

### 2c. 2PL (item discrimination λ<sub>j</sub>)

| | posterior | CrI | Rhat |
|---|---:|---|---|
| σ_α | 2.12 | [1.95, 2.31] | 1.02 |
| σ_λ | 0.275 | [0.24, 0.31] | 1.01 |
| π_α | 0.94 | — | 1.02 |

**Key finding:** σ_α did **not** shrink under 2PL (prediction falsified).
σ_λ is real but modest (~1.7× range across words at ±1 SD). **Item
discrimination is not absorbing the Rasch σ_α inflation.**

By class, highest λ in verbs and adjectives, lowest in nouns and
"other" — opposite of the concrete-noun-sharpness prediction.

### 2d. 2PL + per-child slopes ζ<sub>i</sub>

| | posterior | CrI | Rhat |
|---|---:|---|---|
| σ_α | 2.12 | [1.95, 2.31] | 1.02 |
| σ_λ | 0.276 | — | 1.01 |
| **σ_ζ** | **0.16** | **[0.007, 0.57]** | **1.15** |

**Key finding:** σ_ζ is effectively unidentified from cross-sectional
data. Rhat 1.15, n_eff 23, posterior hugs zero with a long right
tail. The heteroskedasticity argument for identifying per-child slopes
doesn't work in practice with a 14-mo age range.

Also: new "within-age ability SD" PPC panel shows observed SD rises
from ~1.0 at 16 mo to ~1.6 at 30 mo, while model predicts a flat ~2.2.
Confirms σ_α is inflated relative to real within-age dispersion.

**Artifacts:** `model/fits/wordbank_*.rds`, `model/figs/ppc_*.png`.

---

## 🟢 3. Longitudinal linear mixed model on admin-level totals

**Setup:** Wordbank `admins.feather` + `items.feather`. Filter to
English WS longitudinal children (≥2 admins): 1 653 kids, 3 786 admins,
median 2 admins per child, median span 7 mo. Outcome: logit of
production proportion per admin. Fit:
`logit_prop ~ log_age + (log_age | child_id)` via `lme4::lmer`,
log-age centered at 24 mo.

**Result:** fixed effects (intercept −0.26, slope 7.14 per log-age
unit). Random: SD(intercept) 1.49, **SD(slope) 2.43**,
**cor(intercept, slope) +0.72**, residual SD 0.92. LRT for adding
random slope: p < 5 × 10⁻⁶¹.

**Key finding:** children *do* genuinely vary in vocabulary growth
rate with a real, statistically massive effect — but ξ and ζ are
strongly positively correlated, not independent as our Stan model
assumes. Consistent with Peekbank's "faster processors have both
higher current vocab AND faster vocab growth" finding.

**Artifacts:** `model/fits/longitudinal_slopes.rds`.

---

## 🟢 4. σ_r sensitivity sweep (2PL)

Two passes: (4a) on the small laptop subset, (4b) replicated on the
full CDI:WS data on Sherlock. Both tell the same story.

### 4a. Laptop subset (500 × 200, 2 chains × 1500 iter)

| σ_r | σ_α | π_α | σ_xi² |
|---:|---:|---:|---:|
| 0.30 | 2.18 | **0.981** | 4.84 |
| 0.53 | 2.13 | **0.941** | 4.82 |
| 0.80 | 2.05 | **0.868** | 4.84 |
| 1.20 | 1.83 | **0.699** | 4.79 |

### 4b. Full CDI:WS replication (Sherlock, 4 chains × 1500 iter)

Quantitatively almost identical — confirms the decomposition is robust
to sample size and well-mixed:

| σ_r | σ_α | π_α |
|---:|---:|---:|
| 0.30 | 2.18 | **0.98** |
| 0.53 | 2.13 | **0.94** |
| 0.80 | 2.05 | **0.87** |
| 1.20 | 1.83 | **0.70** |

CrI ribbons on π_α are very tight at every σ_r (±0.02–0.04); σ_α has
overlapping CrIs across σ_r values, indicating that the data identify
σ_xi² ≈ 4.8 regardless of how the prior splits it.

### Key finding

Total child-level variance σ_xi² = 4.8 is **tightly identified** by
the data, **but the decomposition into input vs. efficiency is
entirely determined by σ_r**. π_α ranges from 0.70 to 0.98 across
plausible σ_r values. **Even at the widest σ_r = 1.2, input explains at
most 30% of child-level variance** — a robust lower bound on the
efficiency share.

The shape of the curve is informative: π_α moves only gently between
σ_r = 0.3–0.8, then steepens. This is the geometry of σ_α²/(σ_α²+σ_r²):
when σ_α ≫ σ_r, doubling σ_r barely shifts the proportion; when σ_r
approaches σ_α, π_α becomes more sensitive.

### Paper-ready phrasing

> Depending on how much input rates vary across the population the
> CDI sample represents (σ_r ∈ [0.3, 1.2]), individual differences in
> learning efficiency account for **70–98% of between-child variance in
> vocabulary size at 16–30 months**, with input rate explaining the
> complement. Total child-level variance on the logit scale is tightly
> identified at σ_xi² ≈ 4.8; its decomposition between input and
> efficiency is prior-bound.

### What this rules in / out

- **Rules out**: "child differences are mostly about input quantity" —
  that requires σ_r ≳ 2, far outside any defensible value from
  Sperry / Hart & Risley / Weisleder & Fernald.
- **Rules in**: "child differences are predominantly individual
  differences in how efficiently exposure becomes vocabulary." Matches
  the Peekbank picture of stable processing-speed traits being the
  dominant axis of individual variation.
- **Open**: the specific number (70 vs 94%) hinges on what σ_r really
  is for a Wordbank-like sample. Tightening this is the BabyView /
  Seedlings within-vs-between-child input decomposition in the backlog.

**Artifacts:**
`model/fits/sensitivity_sigma_r_2pl.rds`,
`model/fits/wordbank_2pl_sigmaR_*.rds` (full-data replication),
`model/figs/sensitivity_*_2pl.png`.

---

## 🟢 5. Longitudinal accumulator fit (English, 2PL + slopes)

**Setup:** Wordbank longitudinal English WS, *N* = 318 400 obs (600
children × 1 554 admins × 200 items, median 2.6 admins/child). Uses
`log_irt_long.stan` — bivariate (ξ, ζ) per child (LKJ on correlation)
and admin-level age indexing. 4 chains × 1500 iter × 750 warm-up,
adapt_delta 0.95. Sampling time 4 hr.

**Posterior summary:**

| param | median | 95% CrI | Rhat | n_eff |
|---|---:|---|---:|---:|
| σ_α | 2.18 | [2.00, 2.36] | 1.03 | 118 |
| σ_xi | 2.25 | [2.07, 2.42] | 1.03 | 118 |
| σ_zeta | 2.58 | [2.34, 2.86] | 1.04 | 94 |
| **ρ_ξζ** | **0.50** | **[0.43, 0.57]** | 1.02 | 193 |
| π_α | 0.944 | [0.933, 0.952] | 1.03 | 119 |
| s | 7.41 | [6.42, 8.31] | 1.06 | 58 |
| δ | 6.35 | [5.83, 6.84] | 1.06 | 68 |
| σ_λ | 0.30 | [0.27, 0.33] | 1.01 | 341 |

**Key findings:**

1. **σ_ζ = 2.58** — now well-identified (was 0.16 unidentified in the
   cross-sectional fit). Per-child growth-rate variation is *real* and
   of comparable magnitude to σ_α. This was the main test — and the
   cross-sectional identifiability failure, not an absence of signal.

2. **ρ_ξζ = 0.50** — intercept and slope are genuinely coupled, not
   independent. Consistent (on a different scale) with the LMM result
   of +0.72 on admin totals. The model must allow correlated random
   effects.

3. **s moved from 12 to 7 mo** and δ rose from 2.3 to 6.4. Once
   genuine per-child slope variance is in the model, the population
   structural parameters occupy more physically plausible values
   (production plausibly starts around 7 mo; single population δ need
   not absorb all child-level slope heterogeneity).

4. **σ_α and π_α are stable** at ~2.18 and 0.94 — consistent with every
   cross-sectional fit. The variance-decomposition conclusion is robust
   to adding a third (growth-rate) component.

**Implications for the theory.** There are now **three** child-level
variance components in the fitted model:

- input rate σ_r² (pinned externally; 0.28 at baseline σ_r = 0.53)
- level-of-efficiency σ_α² (4.75)
- growth-rate deviation σ_ζ² × L(age)² (up to ~6.7 × L² — varies by age)

"Individual differences in how well a child learns language" is not a
single scalar; it has a level component (ξ) and a rate component (ζ)
that are moderately correlated (0.5). This is the Peekbank picture
operationalized inside the accumulator model.

**Caveats.** Rhat 1.03–1.06 on structural parameters; n_eff 58–200.
Chains are directionally reliable but not fully mixed — longer runs
will tighten intervals but not move medians much.

**Artifacts:** `model/fits/long_2pl_slopes.rds`.

**Posterior-predictive checks (marginal).** Two PPCs, both sampling
new hypothetical children from the fitted bivariate-normal population
distribution (not using any real child's inferred posterior):

- *Within-age SD of logit(vocab/J) vs. age*: smooth predicted ribbon
  rises from SD ≈ 1.0 at 15 mo to ≈ 1.8 at 30 mo. Observed bin-level
  SDs scatter around the ribbon; the few outliers (ages 20, 21, 25, 26)
  all have small *n* (< 30 per bin) and are consistent with sampling
  noise. Large-*n* bins (416, 364, 142 admins at 16, 28, 30 mo) all
  sit essentially on the curve.
- *Trajectory spaghetti, matched age schedules*: hypothetical children
  simulated at the same admin ages as real children with ≥ 3 admins.
  Spread, fanning, and low-age floor all match; the model's
  deterministic within-child growth curve plus item-level Bernoulli
  noise misses some of the within-child irregularity in real trajectories,
  but reproduces the population-level shape.

**Artifacts:** `model/figs/long_2pl_slopes_ppc_variance_marginal.png`,
`model/figs/long_2pl_slopes_ppc_spaghetti_marginal.png`.

---

## 🟢 6. Pivot to lean baseline + variant grammar refactor

**What changed.** After §5 we had a working `long_2pl_slopes` fit but
diagnostic clutter: 2PL discrimination averaged out at the population
level, start time `s` was poorly identified, and per-child slopes
needed longitudinal data to identify. We rewrote `DEFAULT_PRIORS` so
the new lean baseline is Rasch + frequency + per-class psi + free δ,
with `s`, `σ_λ`, `σ_ζ` all pinned at zero via tight priors. Variants
opt in to extra components (see explainer §"Lean baseline + opt-in
variants").

**Implementation.** `variant_hyperpriors()` in `helpers.R` is the
single source of truth; the `long_` and `io_` prefixes are stripped
inside so all pipelines share variant names.

**Artifacts:** explainer `Lean baseline + opt-in variants` and
`Datasets and observation channels` sections.

---

## 🟡 7. Cross-language ablation set (Wordbank longitudinal, lean)

**Setup.** Lean baseline + three ablations (drop slopes, pin δ=0,
free s) × two languages (English, Norwegian). Plus the existing
`long_2pl_slopes` fit on disk for the 2PL comparison. 4 chains × 1000
iter × 500 warm-up, adapt_delta 0.95 throughout. Submitted via
`sherlock/submit_ablations.sh`.

| dataset | variant | reference vs. ablation | status |
|---|---|---|---|
| english | `long_slopes` | reference | 🟢 done (3:13) |
| english | `long_baseline` | drops ζ_i | 🟢 done (3:16) |
| english | `long_fix_delta_slopes` | pins δ=0 | 🟡 running |
| english | `long_free_s_slopes` | frees s | 🟡 running |
| english | `long_2pl_slopes` | adds λ_j | 🟢 from §5 |
| norwegian | `long_slopes` | reference | 🟡 running (~7 hr) |
| norwegian | `long_baseline` | drops ζ_i | 🟡 running |
| norwegian | `long_fix_delta_slopes` | pins δ=0 | 🟡 running |
| norwegian | `long_free_s_slopes` | frees s | 🟡 running |
| norwegian | `long_2pl_slopes` | adds λ_j | 🟢 already on disk |

**Bug history.** First submission crashed: cross-sectional bundle
build had a J/cc dimension mismatch (item-class vector longer than J)
because WG and WS forms have item-name collisions. Fixed in
`prepare_longitudinal_data.R` by deduplicating `word_info` to one row
per `jj`. Second issue: `s_prior_mean = 0` was on the parameter
boundary, blowing up Stan's bounded-parameter transform; fixed to
`(0.5, 0.05)`. See commits `329b8e4` and `d6d0fad`.

**Pending analysis.** Once all six are in: cross-language comparison
plot of fitted vs. observed admin-level vocab vs. age, with each
ablation as a separate panel. Confirms whether the lean baseline
captures population structure as well as the heavier `long_2pl_slopes`,
and whether the ablations break in the expected directions.

---

## 🟡 8. Input-observed fits (BabyView + SEEDLingS)

**Setup.** Same lean baseline + slopes (`io_slopes` variant), but
running on the input-observed Stan model `log_irt_io.stan` which
adds a per-recording measurement layer on `log r_obs` for each
child. Per-bundle differences:

- **BabyView** (head-mounted video, on-camera observer):
  β_react ~ N(0.4, 0.4) — active inflation parameter.
- **SEEDLingS** (LENA all-day audio, passive observation):
  β_react pinned at 0 via N(0, 0.001). Critical correction: passive
  recording has no observer effect, so reactivity inflation does
  not apply. Caught after submitting the first version with the
  BabyView prior (cancelled at 5:50 elapsed; resubmitted with the
  fix).

**BabyView (`io_slopes/babyview`, 23140187):** 🟢 done in 41 min.
20 children × 101 admins × 200 items × 5 688 video recordings.
Posterior medians: σ_α = 1.20, π_α = 0.86, β_react = 0.31 (95% CrI
[-0.42, +1.04]; barely identified given small N), σ_within = 0.70.

**SEEDLingS (`io_slopes/seedlings`, 23177145):** 🟡 pending. 44
children × 514 CDI admins × 200 items × 525 LENA recordings. CDI
input-data path: `cdi_ht_raw_temp.csv` from
[BergelsonLab/WordExposure](https://github.com/BergelsonLab/WordExposure)
(Dong & Bergelson 2026); auto-mapper resolves all 396 WG short
codes; SeedlingsFinalSample filter restricts to the canonical 44
Egan-Dailey subjects.

**Artifacts:** `model/fits/io_slopes.rds` (BabyView), pending
`io_slopes_seedlings.rds`.

---

## 🟡 9. Joint vocab + LWL processing channel (Stanford-linked)

**Setup.** New Stan model `log_irt_long_proc.stan` extends
`log_irt_long.stan` with a second observation channel: per-LWL admin
log RT modeled as a function of per-child log_α and a per-child
RT-by-log-age slope. log_α is now a free per-child latent (not a
shrinkage estimator), drawn jointly with ζ and rtslope from a 3-D MVN
with LKJ. log_r_dev is independent ~ N(0, σ_r), so
ξ_i = μ_r + log_r_dev_i + log_α_i. The LWL channel breaks the log_r
vs log_α exchangeability that holds in CDI-only fits.

**Bridge.** Used `lab_subject_id` from peekbankr 2022.1 to join
adams_marchman_2018 LWL admins to the Stanford TotLot 3 item-level
CDIs we received from Marchman. 62 of 69 Adams kids matched (the
seven unmatched are post-2018 enrollees not in the file we have);
none of the fmw_2013 kids match because their lab IDs (20xxx) don't
overlap with TotLot 2/3 (11xxx). After fixing a subject-level
matching bug, sample is 62 kids × 247 LWL admins (ages 13-20 mo) +
102 CDI admins × 200 stratified items.

**Smoke fit (laptop, 2 chains × 200 sampling iter):**

| param | mean | 95% CrI |
|---|---:|---|
| **γ_rt** | **0.073** | **[0.032, 0.117]** |
| μ_rtslope | -0.73 | [-1.21, -0.28] |
| σ_rtslope | 0.74 | [0.49, 1.03] |
| σ_α | 1.72 | [1.42, 2.11] |
| π_α | 0.91 | [0.88, 0.94] |
| σ_lwl | 0.23 | [0.21, 0.26] |
| ρ_α_ζ | -0.17 | [-0.44, +0.11] |
| ρ_α_rtslope | -0.27 | [-0.71, +0.21] |

**γ_rt is positive and bounded away from zero.** Each unit of log_α
maps to ~0.073 of log RT; a 1-SD higher-α child (~1.7 units) has
~12% lower mean RT. Modest but real -- this is the first quantitative
estimate of how much LWL processing speed informs the model's
efficiency latent.

μ_rtslope = -0.73 confirms the typical declining-RT pattern in the
13-20 mo window, and σ_rtslope = 0.74 means there is meaningful
per-child variation in how fast RT improves.

**Caveats / pre-bug-fix history.** A first version of the linkage
script joined LWL admins to the peekbank-development d_sub records
by (dataset, age) only -- not by subject_id -- so RT and accuracy
attached to a kid's lab_subject_id silently came from a different
child of the same age. Smoke fit then returned γ_rt ≈ 0 with CrI
[-0.03, +0.03]; production submission was cancelled. Fix: pull a
fresh per-admin LWL summary directly from peekbank 2022.1 with
lab_subject_id attached (`pull_peekbank_lwl.R`). After fix, γ_rt
recovers as above. Lesson: never trust an age-only fuzzy join.

**Pending.** Production fit on Sherlock (job 23190113, ~30-60 min).
Will give proper Rhat/ESS, tighter intervals, and a properly-mixed
posterior on the 3-D MVN correlations.

**Artifacts:** `model/stan/log_irt_long_proc.stan`,
`model/fits/long_proc_slopes.rds` (smoke; will be overwritten by
production fit), `data/raw_data/peekbank/peekbank_2022_lwl_summary.csv`,
`data/raw_data/peekbank/peekbank_stanford_linked.csv`.

---

## 🟡 10. Difficulty-side ablations (in progress)

**Motivation.** A diagnostic pass over the existing ablation set
(`compare_english_ability.R`, `compare_english_difficulty.R`) revealed
that 4 of 5 ablations live on the **ability** side of the 2PL
factorization
$\eta_{ijt} = \lambda_j (\theta_{it} - \beta_j)$, with
$\theta_{it} = \xi_i + (1{+}\delta{+}\zeta_i)\log\!((t-s)/a_0)$ and
$\beta_j = \psi_j - \log p_j - \log H$. Only the 2PL variant touches
the item-side, and even then through $\lambda_j$ (multiplier on the
gap), not through the structure of $\beta_j$ itself. The diagnostic
$\beta_j$ density across the 5 existing variants showed near-perfect
correlation (r > 0.995) with small parallel translations: free_s
shifts $\beta_j$ down by ~1.3 logits — a side-effect of ability-side
changes, not a probe of difficulty structure.

**Two new variants** to fill the gap:

| variant | change | what it tests |
|---|---|---|
| `no_class_slopes` | data override: cc<-1, C<-1 (single global $\psi \sim N(\mu, \tau)$) | does lexical-class hierarchy add anything beyond per-word $\psi_j$ + frequency? |
| `class_beta_slopes` | $\beta_c \sim N(1, 0.5)$ free per-class slope on $\log p_j$ | does frequency enter with class-specific weight (e.g., weaker for function words)? |

**Implementation.** `class_beta` adds a new parameter `beta_c` to
`log_irt_long.stan` controlled by the `beta_c_prior_sd` hyperprior in
`DEFAULT_PRIORS` (pinned at 0.001 by default → all $\beta_c$ pinned at
1, equivalent to the pre-change behavior). `no_class` is a data-side
override applied via `variant_data_overrides()` in `helpers.R`.

**Smoke fits** on a 30-admin English subset (2 chains × 60 iter,
laptop) confirmed:
- `long_slopes`: $\beta_c$ pinned at [1.00, 1.00, 1.00, 1.00] ✓
- `long_class_beta_slopes`: $\beta_c$ freed to [0.63, 0.42, 0.33, 0.01]
  (subset too small to interpret meaningfully; confirms identifiability)
- `long_no_class_slopes`: data override applied (C=1), runs cleanly

**Submission.** `sherlock/submit_difficulty_ablations.sh`. Will
populate panels alongside the existing ablation comparison plots once
fits complete.

---

## Backlog (⚪)

### Data / robustness
- **Stanford TotLot CDI mapping — hand review.** The auto-mapper in
  [`model/scripts/parse_stanford_cdi.R`](../model/scripts/parse_stanford_cdi.R)
  resolves all 1 076 short codes (680 WS + 396 WG) to Wordbank
  `item_definition`s with status `auto_exact` (≈75 %) or
  `manual_disambig` (≈25 %, mostly the deterministic Marchman
  disambiguator suffixes: `chicken1`/`chicken2`, `ifconn`, `withprep`,
  `notquant`, etc.). All entries are used in production. Loose end: a
  ~20-minute eyeball pass over
  [`data/raw_data/peekbank/cdi_short_code_map_ws.csv`](../data/raw_data/peekbank/cdi_short_code_map_ws.csv)
  and `cdi_short_code_map_wg.csv` to confirm the manual_disambig rows
  (especially for body-parts compounds, helping verbs with slashed
  forms, and place-names with `*` annotations) match what the form
  actually printed. Replace any wrong mapping in
  `manual_overrides` and rerun.
- **Norwegian longitudinal fit.** Data already pulled (1 562 kids,
  4.2M rows). Needs CHILDES-matched word frequencies; after the English
  longitudinal fit lands, adapt.
- **Within-child input variance** from BabyView or Seedlings: use it
  to constrain σ_r more tightly, replacing the pooled Sperry/HR/WF
  estimate. Would pin down the RQ4 answer inside the sensitivity range.
- **Other-language cross-sectional robustness fits** for languages
  with CDI + CHILDES coverage, to test whether π_α and the RQ1 pattern
  are English-specific.

### Model extensions
- **Correlated (ξ, ζ) prior with LKJ.** Implemented in
  `log_irt_long.stan` already — verify it doesn't collapse under the
  data.
- **Comprehension vs. production joint fit** on the WG form. Bivariate
  (ξ<sup>comp</sup>, ξ<sup>prod</sup>) with estimated correlation; tests
  whether "ability" has modality-specific components.

### Instrumentation
- **Observable.js app** for interactive exploration of the fitted
  model once the posteriors stabilize. Sliders for σ_r, σ_α, σ_ζ, s,
  δ → live growth curves and distributions.
