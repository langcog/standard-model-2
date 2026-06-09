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

**Artifacts:** `fits/recovery.rds`, `outputs/figs/recovery_*.png`.

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

**Artifacts:** `fits/wordbank_*.rds`, `outputs/figs/ppc_*.png`.

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

**Artifacts:** `fits/longitudinal_slopes.rds`.

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
`fits/sensitivity_sigma_r_2pl.rds`,
`fits/wordbank_2pl_sigmaR_*.rds` (full-data replication),
`outputs/figs/sensitivity_*_2pl.png`.

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

**Artifacts:** `fits/long_2pl_slopes.rds`.

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

**Artifacts:** `outputs/figs/long_2pl_slopes_ppc_variance_marginal.png`,
`outputs/figs/long_2pl_slopes_ppc_spaghetti_marginal.png`.

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

## 🟡 8. Input-uptake fits (BabyView + SEEDLingS)

**Setup.** Same lean baseline + slopes (`io_slopes` variant; note the
`io_` code prefix is retained as engineering shorthand for the dual
observation channels — input *and* CDI output — but the published
naming for this experiment is *input-uptake*, following the
literature). Runs the input-uptake Stan model `log_irt_io.stan`,
which adds a per-recording measurement layer on `log r_obs` for each
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

**Artifacts:** `fits/io_slopes.rds` (BabyView), pending
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
`fits/long_proc_slopes.rds` (smoke; will be overwritten by
production fit), `data/peekbank/peekbank_2022_lwl_summary.csv`,
`data/peekbank/peekbank_stanford_linked.csv`.

---

## 🟡 10. Difficulty-side ablations (in progress)

**Motivation.** A diagnostic pass over the existing ablation set
(`compare_english_ability.R`, `compare_english_difficulty.R`) revealed
that 4 of 5 ablations live on the **ability** side of the 2PL
factorization
$\eta_{ijt} = \lambda_j (\theta_{it} - \beta_j)$, with
$\theta_{it} = \xi_i + (1{+}\delta{+}\zeta_i)\log\!((t-s)/a_0)$ and
$\beta_j = \delta_j - \log p_j - \log H$. Only the 2PL variant touches
the item-side, and even then through $\lambda_j$ (multiplier on the
gap), not through the structure of $\beta_j$ itself. The diagnostic
$\beta_j$ density across the 5 existing variants showed near-perfect
correlation (r > 0.995) with small parallel translations: free_s
shifts $\beta_j$ down by ~1.3 logits — a side-effect of ability-side
changes, not a probe of difficulty structure.

**Two new variants** to fill the gap:

| variant | change | what it tests |
|---|---|---|
| `no_class_slopes` | data override: cc<-1, C<-1 (single global $\psi \sim N(\mu, \tau)$) | does lexical-class hierarchy add anything beyond per-word $\delta_j$ + frequency? |
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

## 🟢 11. Ability-side tradeoff diagnostics

`model/scripts/ability_tradeoffs.R` produces two figures:

**(A) (s, δ) joint posterior, lean_ref vs free_s, both languages.**
With s pinned, modest negative correlation r ≈ -0.3. With s freed,
the posterior collapses onto a 1-d ridge: r(s, δ) = -0.96 (English),
-0.95 (Norwegian). s and δ are essentially one parameter when both
free; (s=0.5, δ=9.4) and (s=3, δ=8.1) are equivalent points on the
same ridge. **Practical (near-)non-identifiability.** In principle the
pair is identifiable from the *shape* of (1+δ)·log((t-s)/a_0) — taking
the derivative gives (1+δ)/(t-s), which depends differently on s and
δ. Numerically over the data range (16-30 mo), shape difference
between two ridge endpoints is ~0.1 logits — at the noise floor of
Bernoulli observations after ξ_i absorbs the level shift. Pinning s
in the lean baseline is therefore the right default; `free_s` is a
robustness check, not a refinement.

**(B) (ξ_i, ζ_i) per-kid scatter, English vs. Norwegian.** Raw r(ξ, ζ):
+0.39 (English), -0.32 (Norwegian). The flip survives subset filters
(n_admins ≥ 4, no-ceiling), but **vanishes when ξ is re-centered at
each kid's median admin age** (English +0.30, Norwegian +0.08). The
flip is a parameterization artifact: ξ_i is "ability at a_0=20"; when
admin ages skew far from a_0 (Norwegian median = 26 mo, 80% above
a_0), the per-kid posterior ridge in (ξ, ζ) tilts and pulls the
marginal r across kids. Consistent with the Wordbank-book finding
that variance structure (MADM ≈ 1) is universal across languages.

**Action taken.** All `prepare_*.R` scripts now set `a_0` to
`round(median(admin_info$age))` per dataset, so ξ_i is the per-child
posterior at the natural pivot of the data. This is a
reparameterization (likelihood unchanged); existing fits are still
valid but their ξ_i posteriors are interpretable at the old a_0=20.
Re-fitting under the new bundles is queued for after the
difficulty-side jobs land. See `model_explainer.tex` "Why we pin s"
and "Reference age a_0 is dataset-specific" for the durable writeup.

---

## 🟢 12. Input rate is age-invariant (descriptive check)

`model/scripts/input_rate_vs_age.R` regresses per-recording
$\log r_{iv}^{\mathrm{obs}}$ on age within each kid, on both BabyView
(head-cam, FEM-derived adult tokens) and SEEDLingS (LENA AWC).

| dataset | n kids | n recordings | mean within-kid slope | median | SD |
|---|---:|---:|---:|---:|---:|
| BabyView | 20 | 5,688 | −0.006 logits/mo | −0.002 | 0.058 |
| SEEDLingS | 44 | 525 | −0.011 logits/mo | −0.010 | 0.037 |

Pooled (between+within) slope: BabyView +0.011 (p < 1e-7), SEEDLingS
−0.012 (p ≈ 0.06). Both are an order of magnitude smaller than the
$(1+\delta)\log_{age}$ effect over the data range (~5 logits in
$\theta$ vs ~0.25 logits from input-rate growth).

**Verdict.** The model's age-invariance assumption on $\log r_i$ is
empirically defensible. Per-child slope variance $\sigma_\zeta$ is
therefore mostly attributable to age-varying *efficiency*, not
age-varying input rate. Documented in `model_explainer.tex`
§Assumptions.

**Artifacts:** `outputs/figs/io/input_rate_vs_age.png`.

---

## ⚪ 13. Adopt-and-augment plan across datasets (queued)

Once the M0..M5 nested family on English completes, the plan for
extending across datasets is:

1. **English longitudinal (M0..M5 + no_freq + LMM)** — primary
   structural ablation. LOO-CV + posterior shifts identify M_best.
   *In flight on Sherlock as of 2026-05-02.*

2. **Norwegian longitudinal (M_best only)** — cross-language
   replication. Confirms whether the structural choice from English
   transfers to Norwegian, where δ ≈ 11.5 vs English's 9.4 hints at
   stronger acceleration and Norwegian has more longitudinal
   density per kid (8 admins vs 3).

3. **Input-uptake (M_best + uptake channel) on BabyView and
   SEEDLingS.** The `log_irt_io.stan` Stan file gets patched to
   incorporate whatever components M_best adds beyond the current
   spec (likely `time_baseline`, `beta_c`, `log_lik`). Each dataset
   gets one primary fit; per-dataset secondary ablations (β_react
   free vs pinned in BabyView; class-specific σ_within if signal
   suggests it) are local follow-ups.

4. **Processing (M_best + proc channel) on Peekbank-Stanford.**
   `log_irt_long_proc.stan` similarly patched. Estimates γ_rt for
   the LWL-RT-as-α coupling alongside M_best's structure.

This avoids running M0..M5 across all four data settings, which
would be ~24 fits with mostly redundant findings (the core RQs are
properties of the IRT/accumulator structure, not the dataset).
Dataset-specific structural variation can be tested as targeted
secondary analyses if results suggest it.

---

## 🟢 14. M0..M5 nested-family LOO comparison (English longitudinal)

**Setup.** All six spine variants fit on Sherlock with `log_lik` in
generated quantities. `extract_summaries.R` (job 23696322) pulled
small artifacts (`<tag>.summary.rds`, `<tag>.draws.rds`,
`<tag>.loo.rds`) per fit; rsync'd home into
`fits/summaries/`. Ablation analysis in
`model/scripts/nested_family_analysis.R`.

**Variant → tag mapping.**

| label | tag | what it adds |
|---|---|---|
| M0 | `long_m0` | minimal IRT (no time, no freq) |
| M1 | `long_m1` | unit time + unit freq, δ pinned |
| M2 | `long_baseline` | + free δ |
| M3 | `long_slopes` | + per-child slopes ζ_i |
| M4 | `long_class_beta_slopes` | + class-specific β_c on log p_j |
| M5 | `long_m5` | + 2PL discrimination λ_j |
| no_freq | `long_no_freq_slopes` | M3 with log p_j dropped (RQ2 test) |

**Headline scalar posteriors (English).**

| | σ_α | π_α | δ | σ_ζ | σ_λ | s |
|---|---:|---:|---:|---:|---:|---:|
| M0 | 0.88 | 0.73 | pinned | pinned | pinned | 0.51 |
| M1 | 0.88 | 0.73 | pinned | pinned | pinned | **1.47** |
| M2 | 1.76 | 0.92 | 9.13 | pinned | pinned | 0.52 |
| M3 | 1.75 | 0.91 | 9.39 | 3.48 | pinned | 0.52 |
| M4 | 1.72 | 0.91 | 9.36 | 3.45 | pinned | 0.53 |
| M5 | 1.62 | 0.90 | 8.52 | 3.00 | 0.35 | 0.54 |

The s = 1.47 in M1 is a misspecification diagnostic, not a bug: with
both δ and time-baseline fixed, the model has no other knob for
trajectory shape, so s breaks past its tight prior of (0.5, 0.05) by
~20 SD to mimic what δ would otherwise do. In every well-specified
variant (M2–M5) s sits comfortably at ~0.52, exactly where its prior
puts it.

**Step-wise LOO ELPD differences.**

| step | ΔELPD | SE | z | what it adds |
|---|---:|---:|---:|---|
| M1 vs M0 | +8 033 | 36 | 222.6 | unit time + freq |
| **M2 vs M1** | **+23 289** | **178** | **131.0** | **free δ (RQ1)** |
| M3 vs M2 | +2 212 | 69 | 31.9 | per-child slope ζ_i |
| **M4 vs M3** | **+1.5** | **1.25** | **1.2** | **class-specific β_c (RQ2 fail)** |
| M5 vs M4 | +1 079 | 52 | 20.8 | 2PL λ_j |

Off-spine RQ2 robustness check:

| | ΔELPD | SE | z |
|---|---:|---:|---:|
| **no_freq vs M3** | **+0.9** | **1.6** | **0.5** |

Best by LOO: **M5**.

**Findings by RQ.**

- **RQ1 (acceleration). Closed.** Adding free δ is the largest single
  step in the family by a large margin (z = 131). δ in M2 is 9.13
  [9.01, 9.26] — decisively above zero. The McMurray–Mitchell argument
  that the observed growth-curve shape requires acceleration is
  empirically confirmed at industrial significance.
- **RQ2 (frequency reconstructs word difficulty). Sharp null.** Both
  the M4-vs-M3 test (class-specific β_c, z = 1.2) and the no_freq-vs-M3
  test (drop log p_j entirely, z = 0.5) say frequency contributes
  nothing structural over and above the per-word ψ_j parameter. The
  M4 within-fit posterior gives β_c = (0.55, 0.45, 0.36, 0.16) with
  CrIs comfortably below 1 — but cross-validation says ψ_j and
  β_c · log p_j trade off so cleanly that any combination preserving
  their sum predicts equally well out-of-sample. Frequency information
  is preserved through the per-word ψ_j (the
  `r(ψ_j, log p_j) = 0.68` finding from §2a), it just has no separate
  structural channel.
- **RQ3 (input vs efficiency). Closed.** π_α stabilizes at 0.90–0.92
  in every well-specified variant (M2–M5). The σ_r = 0.534 prior from
  the input-estimation work (see §15 / `input_estimation/`) anchors
  this. M0 and M1's lower π_α ≈ 0.73 are an artifact of the model
  being misspecified there.

**Implications for adopt-and-augment.** Strict-LOO-best is M5, but the
gain of M5 over M3 is +1 079 (z = 21), about 20× smaller than the
M2-vs-M1 step. M3 is the cleanest minimal answer to all three RQs;
M5 is a sensible final-table robustness check. M4 should not be
shipped as the production scaffold (it tells the same out-of-sample
story as M3 but with a worse-conditioned interpretation).

**Artifacts.**
- `fits/summaries/long_{m0,m1,baseline,slopes,class_beta_slopes,m5,no_freq_slopes}.{summary,draws,loo}.rds`
- `outputs/figs/longitudinal/nested_family_scalars.png`
- `outputs/figs/longitudinal/nested_family_loo.png`
- `outputs/figs/longitudinal/nested_family_summary.csv`
- `outputs/figs/longitudinal/nested_family_loo_ranking.csv`
- `outputs/figs/longitudinal/nested_family_loo_steps.csv`

---

## 🟢 15. Norwegian nested family LOO (cross-language replication)

**Setup.** Same 7-stage spine as English §14 fit on Norwegian
longitudinal (200 kids × ~1 600 admins, median 8 admins/kid vs.
English's 3). The three structurally heaviest fits
(`long_baseline_norwegian`, `long_slopes_norwegian`,
`long_class_beta_slopes_norwegian`) had to be **resubmitted with
cmdstanr** because the original April rstan fits predated the
`log_lik` patch in `log_irt_long.stan` (commit b65cfcc) and so had no
LOO. Refit jobs 23816977–79 finished May 4; old rstan fits archived
locally as `*_oldnoll.rds`.

**Headline scalar posteriors (Norwegian).**

| tag | δ | σ_α | σ_ζ | π_α | ρ(ξ,ζ) |
|---|---:|---:|---:|---:|---:|
| `long_m0_norwegian` | 0.01 | 1.42 | — | 0.88 | — |
| `long_m1_time_only_norwegian` | 0.01 | 0.03 | 0.03 | 0.00 | — |
| `long_m1_norwegian` | 0.01 | 0.03 | 0.03 | 0.00 | — |
| `long_baseline_norwegian` (+δ) | 10.98 | 1.88 | — | 0.92 | — |
| `long_no_freq_slopes_norwegian` (+ζ, no freq) | 11.43 | 2.10 | 3.74 | 0.94 | −0.32 |
| `long_slopes_norwegian` (+δ+ζ+freq) | 11.47 | 2.10 | 3.72 | 0.94 | −0.31 |
| `long_class_beta_slopes_norwegian` (+β_c) | 11.44 | 2.10 | 3.74 | 0.94 | −0.32 |
| `long_m5_norwegian` (+2PL) | 7.97 | 1.57 | 2.66 | 0.90 | −0.43 |

**Step-wise LOO ELPD differences (Norwegian).**

| step | ΔELPD | SE | z | reading |
|---|---:|---:|---:|---|
| M0 → +time | +59 136 | 233 | 254 | time helps massively (as in English) |
| +time → +freq | −46 | 4.5 | −10 | freq null/slightly hurts (replicates English) |
| **+freq → +δ alone** | **−235** | **86.7** | **−2.7** | **δ alone HURTS (NOT in English)** |
| +time → +ζ alone (no_freq path) | +2 746 | 96 | +29 | ζ alone helps massively |
| +δ → +δ+ζ | +3 028 | 83 | +37 | adding ζ on top of δ is the big win |
| +ζ → +ζ+freq | +1 | 1.5 | +0.7 | freq still null |
| +β_c (Mboth → Mclass) | **−2.7** | **1.1** | **−2.3** | **class-specific β_c null (replicates English RQ2)** |
| → +2PL | +1 425 | 59 | +24 | 2PL real, matches English ~+1 079 |

**Striking finding: δ replicates as exponent, not as parameter.**
In English, freeing δ on top of the unit accumulator is the largest
structural step (+23 289 ELPD, z = 131). In Norwegian, freeing δ
**alone** *reduces* ELPD by 235; freeing per-child slopes ζ_i *alone*
(without δ) gains 2 746 ELPD instead. Once both are free, the
posterior medians are δ = 11.5, σ_ζ = 3.7, almost identical to
English (δ = 9.4, σ_ζ = 3.5). The cross-language disagreement is
about *parameterization*, not about *acceleration*: the population
mean of the per-child scaling exponent (1 + δ + ζ_i) is ~10.4 in
English and ~12.5 in Norwegian — comparable.

**Driver: longitudinal density.** Norwegian has 8 admins/kid (vs.
English's 3), so per-child slopes ζ_i are well-identified.
Population δ becomes redundant when ζ_i can carry the per-kid
acceleration directly. With sparse longitudinal data (English), the
data cannot pin ζ_i and the population term δ pools the action
instead.

**Implication for paper framing.** The "all children show
super-linear scaling" headline is robust across languages. The clean
parameterization for the comparison is `(1 + δ + ζ_i)` — not δ alone.
The disanalogy figure (`outputs/figs/schematic/D1_scaling_disanalogy.png`)
should label the kid scaling exponent as `1 + δ + ζ_i` rather than
`1 + δ`, and the population mean of that quantity is the natural
cross-language summary.

**RQ2 cross-language.** Norwegian replicates English's sharp RQ2 null:
adding class-specific β_c on top of `long_slopes_norwegian` *reduces*
ELPD by 2.7 (z = −2.3, comparable to English's near-zero step).
Frequency contributes nothing structural in either language once
per-word ψ_j is free.

**Artifacts.** `fits/summaries/long_*_norwegian.{summary,draws,loo}.rds`
(8 fits × 3 files each, complete).

---

## 🟢 16. SEEDLingS parent-report-noise correction (comp + std channels)

**Concern.** σ_α in CDI-only fits absorbs parent-report
measurement error. The CDI ↔ later-CELF correlation r = 0.63 implies
CDI reliability ≤ 0.42 (lower bound; developmental drift attenuates
further), suggesting 30–50% of CDI between-child variance may be
parent-report noise rather than true between-child efficiency.
The SEEDLingS extreme `π_α = 0.98` is the most likely candidate to
be inflated by this.

**Approach.** Add second observation channels on the SEEDLingS sample
that pin `log α_i` from non-CDI signals:
- **comp**: SEEDLingS WG comprehension scores enter the io model as
  a second measurement layer parallel to production, sharing the
  per-child latent. Stan extension committed in `f12d2b0` (`comp_*`
  variant family).
- **std**: SEEDLingS preschool-age standardized scores (CELF, QUILS,
  PVT at ~4;6) enter as a non-CDI readout of `log α_i`. Stan
  extension in `42d5a4e` (`std_*` variant family). Same modular
  pattern as `comp`.

**Four-way SEEDLingS comparison.**

| variant | σ_α | σ_ζ | π_α | δ | ESS_bulk(δ) |
|---|---:|---:|---:|---:|---:|
| baseline (`io_no_freq_slopes_seedlings`) | 2.59 [2.13, 3.21] | 3.62 | 0.98 [0.97, 0.99] | 8.05 | ~3 700 |
| **+comp** (`io_comp_no_freq_slopes_seedlings`) | **1.39 [1.13, 1.74]** | **2.73** | **0.94 [0.89, 0.97]** | 7.76 | ~3 600 |
| +std (`io_std_no_freq_slopes_seedlings`) | 0.96 [0.80, 1.23] | 1.92 | 0.88 [0.80, 0.94] | 3.77 | ~3 700 |
| +comp +std (`io_comp_std_no_freq_slopes_seedlings`) | 1.42 [1.16, 1.76] | 2.72 | 0.94 [0.90, 0.97] | 7.78 | ~720 |

**Reading.**
- **+comp drops σ_α by 46%** (2.59 → 1.39): comp identifies
  `log α_i` independently of production noise within the same age
  window, removing a major chunk of parent-report measurement error.
  π_α moves from the 0.98 ceiling to 0.94 — still high, but no longer
  pinned.
- **+std drops σ_α by 63%** (2.59 → 0.96), but also collapses δ from
  8.05 to 3.77. The δ drop is suspicious: standardized scores at
  ~4;6 are far outside the CDI window (6–18 mo), so `log α_i` becomes
  pinned by a future-state measurement and the model has fewer
  degrees of freedom to fit early-window growth, soaking the gap into
  reduced δ. Best read as an identifiability artifact, not a
  substantive estimate. **Don't quote +std alone as the headline
  correction.**
- **+comp+std** matches +comp on σ_α and δ — comp dominates,
  std-channel contribution is absorbed. ESS for δ drops to ~720
  (from ~3 700), suggesting the joint fit has slight identifiability
  trouble but still yields reasonable posteriors.

**Headline correction.** Use **+comp** (not +std, not baseline) as
the SEEDLingS reading: σ_α = 1.39, π_α = 0.94. Cross-sample range
tightens from [0.84, 0.98] to [0.84, 0.94] — still efficiency-
dominated, no longer at the ceiling for any sample.

**Caveat for the abstract.** The "80–90% of between-child variation
unexplained by input quantity" phrasing remains accurate at the
sample level, with the SEEDLingS extreme rationalized post-correction
rather than dismissed. The bigger-picture reading is unchanged:
efficiency variance dominates input-rate variance robustly.

**LOO.** `log_irt_io.stan` does not currently emit `log_lik`, so the
four-way comparison is on scalar posteriors only, not LOO. Adding
log_lik would let us formally compare these four; for the parent-
report-noise question scalar-posterior comparison is sufficient.

**Artifacts.**
`fits/summaries/io_{,comp_,std_,comp_std_}no_freq_slopes_seedlings.{summary,draws}.rds`.

---

## 🟢 17. Peekbank LWL channel: M_best fit (no_freq variant)

**Setup.** Following the §14 finding that frequency contributes
nothing structural beyond per-word ψ_j, the production Peekbank fit
drops the frequency channel: `long_proc_no_freq_slopes` on the 62
Stanford-linked subjects with both LWL admins and item-level CDIs.
Same `log_irt_long_proc.stan` as §9 with `beta_c` pinned at 0.

**Headline scalars.**

| variable | median | 95% CrI |
|---|---:|---|
| σ_α | 1.64 | [1.34, 2.03] |
| σ_ζ | 6.30 | [5.24, 7.43] |
| π_α | 0.90 | [0.86, 0.94] |
| δ | 2.56 | [1.62, 3.44] |
| γ_rt | 0.082 | [0.046, 0.125] |
| μ_rtslope | −0.74 | [−1.09, −0.40] |
| σ_rtslope | 0.76 | [0.47, 1.07] |

**Reading.** π_α = 0.90 confirms the §9 `long_proc_slopes` finding
robustly (0.90 vs. 0.91 with frequency). σ_α drops slightly under the
no-freq simplification but otherwise the structural picture is the
same. The δ value (2.56) is much smaller than English `long_slopes`
(9.39), echoing the Norwegian pattern: when a second high-quality
readout on log α_i is available (LWL here, ζ_i in Norwegian), δ is no
longer the load-bearing parameter.

**γ_rt = 0.082** [0.046, 0.125] — bounded firmly above 0; LWL
processing-speed gain per unit of `log α_i` is real. Each unit of
`log α_i` corresponds to ~8% lower mean log-RT.

**Cross-readout correlation reminder (from earlier work).** ρ(ζ_i,
rtslope_i) = −0.024 [−0.33, 0.27] — null at the noise floor.
CDI-side and LWL-side growth rates do not share variance even at this
sample size. The two channels are independent maturation clocks.

**Artifacts.** `fits/summaries/long_proc_no_freq_slopes.{summary,draws}.rds`
(no LOO; `log_irt_long_proc.stan` doesn't currently emit log_lik).

---

## 🟢 18. σ_r sensitivity for M_best (analytical + one confirmatory refit)

**Why revisit.** §4b ran a 4-point σ_r sweep on the cross-sectional 2PL
model and found σ_ξ² ≈ 4.8 stable across σ_r priors, with π_α
following the analytical formula 1 − σ_r²/σ_ξ². With M_best now being
longitudinal slopes (different identification structure), and the
SEEDLingS comp-correction (§16) demonstrating that pinned model
assumptions can shift π_α, the question is worth re-asking. But the
§4b geometry plausibly extends to M_best, so we did the cheap
analytical extension first and queued one confirmatory refit at the
"challenging" σ_r = 0.8 value.

**Analytical extension to M_best.** For each draw of σ_ξ from the
fitted posterior, π_α(σ_r) = 1 − σ_r²/σ_ξ² is computed across a σ_r
grid; the §4b refit points at matched σ_r values are overlaid as
validation.

| fit | σ_ξ med | π_α(σ_r=0.30) | π_α(σ_r=0.534) | π_α(σ_r=0.80) | π_α(σ_r=1.20) |
|---|---:|---:|---:|---:|---:|
| English `long_no_freq_slopes` | 1.83 | 0.97 | 0.91 | 0.80 | 0.55 |
| English `long_slopes` | 1.83 | 0.97 | 0.91 | 0.81 | 0.57 |
| Norwegian `long_no_freq_slopes_norwegian` | 2.16 | 0.98 | 0.94 | 0.86 | 0.69 |
| Norwegian `long_slopes_norwegian` | 2.17 | 0.98 | 0.94 | 0.86 | 0.69 |
| §4b cross-sectional 2PL (refit) | 2.19 | 0.98 | 0.94 | 0.87 | 0.70 |

The §4b refit values land **on** the analytical curves (Norwegian σ_ξ
and §4b σ_ξ are nearly identical, so curves coincide). This validates
that the σ_ξ posterior is anchored by the data and the σ_r prior just
partitions σ_ξ² into σ_α² and σ_r².

**Reading.**
- At the externally pinned σ_r = 0.534, all four production fits give
  π_α ∈ [0.91, 0.94]. Robust.
- M_best (English) is **slightly more sensitive** than §4b at high σ_r
  (π_α drops to 0.55 at σ_r = 1.2 vs §4b's 0.70) because longitudinal
  σ_ξ ≈ 1.83 is smaller than cross-sectional σ_ξ ≈ 2.19. The
  longitudinal fit attributes more variation to σ_ζ (slopes).
- Norwegian σ_ξ ≈ 2.17 lands near the §4b cross-sectional value, so
  Norwegian's σ_r sensitivity matches §4b almost exactly.
- The qualitative structural claim — efficiency dominates input rate
  in any plausible σ_r regime — survives. At σ_r = 0.8 (the upper end
  of any defensible external estimate), π_α is still 0.80–0.86. Even
  at σ_r = 1.2, English M_best gives π_α ≈ 0.55, well above 0.5.

**Confirmatory refit (job 23881868, completed).** σ_r = 0.8 refit on
`long_no_freq_slopes` via the new `STAN_SIGMA_R_OVERRIDE` env-var
hook in `fit_longitudinal.R` (commit 88b2d4f). Output tag:
`long_no_freq_slopes_sigmaR_0p80`. **Refit and analytical prediction
agree to three decimal places on every relevant quantity:**

| | π_α | σ_α | σ_ξ |
|---|---:|---:|---:|
| Analytical prediction (from σ_r=0.534 fit) | 0.802 [0.762, 0.837] | 1.608 [1.430, 1.815] | 1.796 |
| Refit at σ_r = 0.8 | 0.801 [0.760, 0.839] | 1.607 [1.425, 1.824] | 1.795 [1.634, 1.991] |

σ_ξ is unchanged between the two priors (1.795 vs. 1.796). σ_ζ
similarly stable (3.51 vs. 3.48), δ moved trivially (9.65 vs. 9.39,
within CrI overlap). The σ_r prior just partitions a data-identified
σ_ξ² into σ_α² and σ_r², exactly as the §4b finding said for the
cross-sectional 2PL.

**Conclusion.** The analytical extension is now validated for M_best
as well. We can quote π_α(σ_r) curves from the analytical formula
without further refits; the rest of the §4b sweep is unnecessary on
M_best.

**Caveat.** Refit ESS_bulk for σ_α / σ_ξ / π_α was ~102 (Rhat = 1.03)
— borderline but adequate for the central-tendency validation. If
the abstract leans on tight CrIs at σ_r = 0.8 specifically (as
opposed to just the median), we should rerun with longer chains.
Currently the M_best CrI at σ_r = 0.534 is the headline, so this is
not blocking.

**Headline π_α robustness across σ_r ∈ [0.3, 1.2]:**

> Across plausible external estimates of σ_r (range [0.3, 1.2] from
> Sperry / Hart-Risley / Weisleder-Fernald sensitivity), π_α for
> M_best on English ranges from 0.97 (σ_r = 0.3) to 0.55 (σ_r = 1.2);
> at the externally pinned σ_r = 0.534, π_α = 0.91. For Norwegian
> the corresponding range is [0.98, 0.69] with π_α = 0.94 at the
> external pin. **In every plausible σ_r regime, efficiency variance
> dominates input variance** (π_α > 0.5), and at the external pin
> the dominance is strong (≥ 0.91 in both languages).

**Artifacts.**
- `model/scripts/sigma_r_analytical_sensitivity.R`
- `outputs/figs/longitudinal/sigma_r_analytical_sensitivity.png`
- `outputs/figs/longitudinal/sigma_r_analytical_sensitivity.csv`
- `fits/summaries/long_no_freq_slopes_sigmaR_0p80.{summary,draws}.rds`

---

## 🟢 19. Cross-sample π_α replication (post-correction)

| sample | n | σ_α | π_α | source fit |
|---|---:|---:|---:|---|
| English long_slopes | 200 | 1.83 [1.65, 2.04] | 0.91 [0.90, 0.93] | `long_slopes` |
| Peekbank long_proc | 62 | 1.64 [1.34, 2.03] | 0.90 [0.86, 0.94] | `long_proc_no_freq_slopes` |
| BabyView io | 20 | 1.13 [0.86, 1.61] | 0.84 [0.68, 0.93] | `io_no_freq_slopes` |
| **SEEDLingS io+comp** | 44 | 1.39 [1.13, 1.74] | **0.94 [0.89, 0.97]** | `io_comp_no_freq_slopes_seedlings` |
| SEEDLingS io baseline (uncorrected) | 44 | 2.59 | 0.98 | `io_no_freq_slopes_seedlings` |
| Norwegian long_slopes | 200 | 2.10 [1.90, 2.35] | 0.94 [0.93, 0.95] | `long_slopes_norwegian` |

**Headline.** π_α ∈ [0.84, 0.94] across **five samples** (English,
Peekbank-Stanford, BabyView, SEEDLingS-comp-corrected, Norwegian),
two languages, three input-observation channels, with parent-report
noise correction where available. Robust efficiency-dominated
decomposition; no sample sits at the ceiling once correction channels
are deployed.

---

## 🟢 20. s-identifiability check on M_best (free_s_no_freq_slopes)

**Motivation.** §2b and the explainer §"s near-non-identifiability"
documented an (s, δ) ridge: log_age and time-from-onset are
near-collinear over the 16–30 mo CDI:WS range, so the marginal
likelihood is nearly flat along a ridge in (s, δ). M_best
(`long_no_freq_slopes`) pins s ≈ 0 and reads κ_pop = 1 + δ off the
posterior. Two open questions: (a) does the same data prefer an
interior s when s is freed? (b) is the qualitative
"super-linear scaling" finding (κ_pop > 1) parameterization-invariant,
or does it depend on the s-pinning?

**Setup.** Variant `free_s_no_freq_slopes` (added to
[`model/R/helpers.R`](../model/R/helpers.R)): same as M_best
(`long_no_freq_slopes`) except `s_prior_mean = 4.5, s_prior_sd = 2`
on the truncated-normal prior s ∈ [0, 15]. Frequency channel pinned
at zero (`beta_c_prior_sd = 0.001`), per-child slopes free
(`sigma_zeta_prior_sd = 1`). 4 chains × 1500 iter on Sherlock
(job 23900133, completed).

**Headline result.**

| | M_best (s ≈ 0) | free_s_no_freq_slopes |
|---|---:|---:|
| s | 0 (pinned) | **14.999 [14.997, 15.000]** (boundary) |
| δ | 9.39 [8.79, 10.04] | **1.54 [1.49, 1.58]** |
| σ_α | 1.83 [1.65, 2.04] | 2.08 [1.88, 2.30] |
| σ_ζ | 3.51 [3.10, 3.95] | 1.22 [1.11, 1.35] |
| κ_pop = 1 + δ | 10.4 | 2.54 |
| ESS_bulk (s) | — | 3110 |
| Rhat (s, δ) | — | 1.001, 1.000 |

**Interpretation.**

1. **The posterior settles at the upper prior boundary.** s = 14.999
   with sd = 0.0008 — the free-s mode is *concentrated* at the
   boundary, not failing to mix (Rhat = 1.00, ESS = 3110). Chains agree
   on the s = 15 corner.
2. **Removing the frequency channel weakens the ridge constraint.**
   With frequency channel free (`long_free_s_slopes`, run earlier on
   the same English data), the posterior gives s ≈ 3.05 [2.14, 3.87]
   with δ ≈ 8.15 — a sensible interior solution. The frequency channel,
   even when LOO-null at the model-comparison level, was silently
   anchoring the ridge to the small-s end. With it pinned, the ridge
   tilts to the large-s, small-δ corner.
3. **Direction of acceleration is parameterization-invariant; magnitude
   is not.** Both modes give κ_pop > 1 (10.4 vs. 2.54), so the
   structural claim ("super-linear scaling") is robust. But the
   headline number κ_pop ≈ 10 is conditional on the s-pinning.
4. **Both modes fit the 16–30 mo data about equally well.** The swing
   in η across the empirical age range is similar; they differ
   primarily in (a) extrapolation to ages < 15 mo (the s = 15 mode is
   structurally unable to predict ages below s) and (b) the σ_ξ²
   partition between σ_α² and σ_ζ² (free-s pushes more variance into
   σ_α, less into σ_ζ).
5. **Free-s does not address the panel-4 floor effect.** Inspection
   of the predicted fan at ages 16 and 30 under s = 15 shows a *wider*
   fan at age 16 and a similar fan at age 30 — the wrong direction.
   The floor effect is not a missing-s problem; per-child onset
   variation (s_i random rather than global s) would be the natural
   next extension if we wanted to address it directly.

**Conclusions.**

- Stick with M_best (s ≈ 0) as the headline parameterization. κ_pop ≈
  10.4 is reported with the explicit caveat that the magnitude depends
  on s-pinning (now documented in the explainer §"s
  near-non-identifiability").
- The κ_pop > 1 finding is robust to s-pinning; the magnitude is not.
  This goes in the slide deck as a footnote/sensitivity remark, not a
  retraction.
- A per-child s_i extension is the natural next step to address the
  panel-4 floor effect — explicitly *not* what free-s does.

**Artifacts.**
- [`model/R/helpers.R`](../model/R/helpers.R) (`free_s_no_freq_slopes`
  variant added)
- `fits/summaries/free_s_no_freq_slopes.{summary,draws}.rds`
- Explainer §"s near-non-identifiability" updated with "The ridge
  depends on what else is free." paragraph and implications.

---

## 🟢 21. Per-child onset s_i: a fourth random effect

**Motivation.** §20 documented that global s is near-non-identifiable
and the s-free mode doesn't address the panel-4 floor effect (q10
not anchoring at 0 at age 16, q90 overshooting at age 30). The
natural extension proposed there: per-child onset variation
`s_i ~ Normal_+(0, sigma_s)`, with global s pinned tight near 0.

**Setup.** Two new variants on M_best
(`long_no_freq_slopes`, σ_α + σ_ζ free):

| variant | σ_α | σ_ζ | σ_s | (ξ,ζ,s) correlations |
|---|---|---|---|---|
| `long_no_freq_si_only` | free | pinned 0 | free | independent |
| `long_no_freq_slopes_si` | free | free | free | independent |
| `long_no_freq_slopes_si_corr` | free | free | free | trivariate LKJ(2) |

`si_corr` uses a new Stan file
[`log_irt_long_si_corr.stan`](../model/stan/log_irt_long_si_corr.stan)
with a 3×3 Cholesky on (ξ, ζ, s_lat) and Tobit clipping
`s_i = max(s_lat, 0)`. All variants: HN(0, 2) prior on σ_s.

**Sampling settings.** 4000 iter × 4 chains, 2500 warmup,
adapt_delta = 0.95-0.99. Centered parameterization (non-centered
was tried briefly and produced a bimodal posterior — see notes
in `model/stan/log_irt_long.stan` comments).

**Headline: adding s_i to M_best is a large LOO win.**

PSIS-LOO across the family on the English longitudinal bundle
(N = 145,350 observations; M_best was fit on N = 147,000 and is
shown per-obs-normalized for comparison):

| model | elpd_loo | Δ vs si_corr | SE_diff | p_loo |
|---|---:|---:|---:|---:|
| **slopes_si_corr** | **−36,982** | 0 | — | 636 |
| slopes_si (indep) | −37,044 | −62 | 8 | 616 |
| demo_kappa (ζ alone) | −37,274 | −292 | 23 | 527 |
| si_only (α + s_i) | −37,482 | −500 | 26 | 568 |
| **M_best** (per-obs proj.) | **~−37,820** | ~−535 | — | 577 |
| demo_pure / demo_alpha | ~−62,695 | −25,713 | 175 | ~340 |

**Two findings.**

1. **Adding s_i to M_best earns ~535 elpd by per-obs estimate** —
   far beyond the 90-effective-parameter overfit penalty. The s_i
   mechanism is doing substantive structural work for the model.

2. **LKJ correlations earn another 62 elpd over independent s_i**
   (SE_diff = 8, ~8σ separation). Statistically significant, but
   practically a smaller contribution.

**Posterior summaries (slopes_si independent, headline variant).**

| param | median | 95% CrI | Rhat | ESS_bulk |
|---|---:|---|---:|---:|
| σ_α | 1.86 | [1.64, 2.10] | 1.80 | 6 |
| σ_ζ | 3.10 | [2.85, 3.41] | 1.39 | 9 |
| σ_s | 4.14 | [3.28, 5.01] | 1.83 | 6 |
| δ | 8.06 | [7.73, 8.36] | 1.91 | 6 |
| π_α | 0.92 | [0.90, 0.94] | 1.80 | 6 |

Compared to M_best: δ drops from 9.4 → 8.1 (some of the
"population acceleration" gets redistributed into onset variation),
σ_α and σ_ζ are essentially unchanged. The added σ_s = 4.14 mo
means kid onsets distributed half-normally with scale ~4 — so 50%
of kids in [0, 2.8] mo, 95% in [0, 8.1] mo. Empirically plausible.

**Mixing diagnostic concerns.** Rhat = 1.5–1.9 on key params despite
4000 iter / 2500 warmup / adapt_delta 0.99. Per-chain medians agree
within ~5%, so the chains are exploring the same single mode but
slowly within it. Point estimates are stable; tight CIs are not yet
publication-trustworthy. Earlier non-centered parameterization had
true multimodality (chain 4 in a "sigma_s ≈ 36, sigma_alpha ≈ 0.9"
mode disconnected from the others); centered fixed that.

**Correlation findings (si_corr variant; sidebar).** Three pairwise
correlations on (ξ, ζ, s_lat):

| correlation | median | interpretation |
|---|---:|---|
| ρ(ξ, ζ) | −0.04 | drops to ~0 once s is in the model |
| **ρ(ξ, s)** | **+0.64** | high efficiency ↔ later onset ("catch-up learners") |
| **ρ(ζ, s)** | **−0.63** | high acceleration ↔ earlier onset ("early bloomers") |

The (ξ, ζ) correlation flips from +0.35 in M_best to ~0 here — its
old role is being absorbed by the two new s_i correlations.

These patterns connect to separately-documented Wordbank findings
on onset timing and later vocabulary composition (Wordbank book
[ch. 13.3.3](https://langcog.github.io/wordbank-book/grammar.html#results-3)
and [ch. 15.2](https://langcog.github.io/wordbank-book/style.html#variation-in-vocabulary-composition)).
But the LOO gain over independent-s_i is modest (62 elpd, SE 8),
and Pareto-k diagnostics are concerning (202 obs with k > 0.7).
**For headline reporting we lean on the independent-s_i variant.**
The correlation structure is real but a sidebar finding.

**Artifacts.**
- Stan models: [`log_irt_long.stan`](../model/stan/log_irt_long.stan) (independent s_i toggle), [`log_irt_long_si_corr.stan`](../model/stan/log_irt_long_si_corr.stan) (trivariate LKJ)
- Variants: `no_freq_si_only`, `no_freq_slopes_si`, `no_freq_slopes_si_corr` in [`helpers.R`](../model/R/helpers.R)
- LOO files: `fits/summaries/long_no_freq_*_si{_corr}.loo.rds` (gitignored, regenerable via `sherlock/extract_loo_thinned.R`)
- Figure: [`outputs/figs/longitudinal/quantile_demo.png`](figs/longitudinal/quantile_demo.png) (7-panel comparison)

---

## 🟢 22. Fixing s_i mixing: reparam + signed-normal

**Motivation.** §21's slopes_si family had Rhat 1.5–1.9 on key
parameters (σ_α, σ_ζ, σ_s, δ) despite single-mode posteriors. Two
overlapping ridges were diagnosed via per-chain median comparison:

1. **(σ_ζ, σ_s) variance-partition ridge**: the two between-kid
   variance channels traded off — when σ_s grew, σ_ζ shrank.
2. **(σ_s, δ) population-mean compensation ridge**: with half-normal
   s_i, E[s_i | σ_s] = σ_s × √(2/π) ≈ 0.8σ_s depends on σ_s. When σ_s
   grew, the average effective onset shifted later, requiring δ to
   compensate.

**Two fixes, applied sequentially.**

**Fix 1: (σ_total, p_zeta) reparameterization** [`log_irt_long_si_reparam.stan`]
- Sample axis-aligned coordinates: σ_total = sqrt(σ_ζ² + σ_s²) (well
  identified) and p_zeta = σ_ζ²/σ_total² ∈ (0, 1) (the loose direction).
- Back-transform: σ_ζ = σ_total × √p_zeta; σ_s = σ_total × √(1−p_zeta).
  Same interpretation as before; only sampling axes differ.
- Variant: `no_freq_slopes_si_reparam` (job 24959733; 9h52m runtime).

**Fix 2: signed-normal s_i with sum-to-zero centering**
[`log_irt_long_si_signed.stan`]
- Switch from `vector<lower=0>[I] s_i; s_i ~ normal(0, σ_s)` (half-normal)
  to `vector[I] s_i_raw; s_i = s_i_raw - mean(s_i_raw); s_i_raw ~ normal(0, σ_s)`
  (signed, sum-to-zero). Now E[s_i] = 0 by construction, regardless of σ_s.
- Interpretation shifts from "literal delay" to "developmental offset
  around population mean" (kids with s_i < 0 are "ahead" — their
  effective log-age at calendar age t is larger).
- Matches Mike's intuition that the half-normal felt artificial: empirical
  age-of-first-word varies symmetrically around a population mean.
- Variant: `no_freq_slopes_si_signed` (job 25093875; 2h53m runtime).

**Mixing fully solved.**

| param | direct (v2) | reparam only | signed + reparam |
|---|---|---|---|
| σ_α Rhat / ESS | 1.80 / 6 | 1.07 / 62 | **1.00 / 889** ✓ |
| σ_ζ Rhat / ESS | 1.39 / 9 | 1.04 / 59 | **1.00 / 3250** ✓ |
| σ_s Rhat / ESS | 1.83 / 6 | 1.21 / 16 | **1.02 / 295** ✓ |
| δ Rhat / ESS | 1.91 / 6 | 1.28 / 12 | **1.00 / 6590** ✓ |
| π_α Rhat / ESS | 1.80 / 6 | 1.07 / 62 | **1.00 / 889** ✓ |

Chain wall-clock times also collapsed from 9–17h spread (direct) → ~10h
synchronous (reparam) → **2:03–2:30 synchronous (signed + reparam)**.

**Posterior shifts compared to half-normal.**

| param | direct (v2) | reparam only | signed + reparam |
|---|---|---|---|
| σ_α | 1.86 | 1.98 | **1.56 [1.33, 1.79]** |
| σ_ζ | 3.10 | 3.30 | **3.51 [3.14, 3.92]** |
| σ_s | 4.14 | 4.62 | **1.40 [0.88, 1.88]** |
| δ | 8.06 | 7.93 | **9.62 [9.47, 9.79]** |
| π_α | 0.92 | 0.93 | **0.894 [0.862, 0.919]** |
| ρ(ξ,ζ) | +0.35 | +0.29 | **+0.43 [+0.25, +0.59]** |
| κ_pop = 1+δ | 9.06 | 8.93 | **10.62** |

The half-normal model was inflating σ_s because all kids were forced
into a single one-sided distribution; once kids can be signed,
true per-kid onset spread is only ~1.4 mo SD. δ snaps back to near
M_best (κ_pop ≈ 10.6 vs M_best's 10.4). σ_α and σ_ζ near M_best values.

**LOO comparison.**

| model | elpd_loo | bad pareto-k | p_loo | Rhat |
|---|---:|---:|---:|---:|
| slopes_si_corr (LKJ+Tobit) | −36,982 | 202 ⚠ | 636 | 2.5 |
| slopes_si (half-normal) | −37,044 | 3 | 616 | 1.91 |
| **slopes_si_signed** | **−37,123** | **0** ✓ | 577 | **1.02** ✓ |
| demo_kappa (ζ alone) | −37,274 | 0 | 527 | 1.00 |
| si_only (α + s_i) | −37,482 | 2 | 568 | 1.01 |
| M_best (per-obs proj.) | ~−37,820 | — | 577 | 1.04 |

Signed-s_i ranks third on raw elpd but is the **only mixed-cleanly
high-elpd variant**. Half-normal beats it by 79 elpd but has Rhat 1.91 —
those point estimates were not trustworthy. The 697-elpd gain over
M_best per-obs projection is the real bottom line: **adding signed s_i
to M_best is a substantial and well-mixed improvement.**

**Substantive interpretation.**

The signed s_i posterior says per-kid developmental onsets vary by
~1.4 mo SD around population mean. That's modest — about ±2.8 mo for
2 SD. Combined with the (nonlinear) log transform: a kid 2 mo "late"
at calendar age 16 gets log((16-0-2)/19) = -0.169 (vs population
log(16/19) = -0.172), so the log-age effect is small at the
population edges but compounds in the tails.

Most importantly: the signed model **restores the canonical M_best
picture** (κ_pop ≈ 10.6, σ_α ≈ 1.6, σ_ζ ≈ 3.5, π_α ≈ 0.89) with
s_i as a small extra channel. The half-normal had pushed κ_pop down
to 9.1 — that was an artifact of bad parameterization, not a real
finding.

**Artifacts.**
- [`model/stan/log_irt_long_si_reparam.stan`](../model/stan/log_irt_long_si_reparam.stan): (σ_total, p_zeta) reparam only
- [`model/stan/log_irt_long_si_signed.stan`](../model/stan/log_irt_long_si_signed.stan): reparam + signed s_i (headline variant)
- Variants `no_freq_slopes_si_reparam`, `no_freq_slopes_si_signed` in [`helpers.R`](../model/R/helpers.R)
- Figure: [`outputs/figs/longitudinal/quantile_demo.png`](figs/longitudinal/quantile_demo.png) (6-panel; panel 6 is signed s_i)

> **Update (2026-05-23).** Sections 23–28 below back out the §22 conclusion.
> After scaling to I=500, signed s_i developed multi-mode posteriors that the
> I=50 pilot had concealed. We excised both `s` and `s_i` from the headline
> model. M_best is now `α + ζ + δ` only.

---

## 🟢 23. Excising `s` and `s_i` from the headline model (2026-05-22)

**TL;DR.** After §22, we tried to scale up the signed-s_i model to I=500. It
broke (multimodal posterior, three chains in different modes). We diagnosed
the failure as a non-smooth boundary in `fmax(t - s - s_i, 0.01)` plus an
unanchored `(s, δ)` ridge, neither of which the smaller pilots had exposed.
Backed out both `s` and `s_i` from the headline model: M_best is now
`α + ζ + δ` only, with `s` pinned tight at 0 (`s_prior_mean=0,
s_prior_sd=0.001` in `DEFAULT_PRIORS`) and `sigma_s_prior_sd=0.001` keeping
`s_i ≈ 0`.

**Why this took so long to surface.** Three separate things lined up at I=500
that didn't at I=50:

1. **(s, δ) ridge re-emerged.** During the s_i campaign we had pivoted
   `s_prior_mean` from 0.5 (historical) to 6 for signed-s_i interpretability
   ("kids with s_i < 0 are ahead of population, kids with s_i > 0 are
   behind"). With `s = 6`, `fmax(t − s − s_i, 0.01)` clamps the youngest
   ~25% of admins (those with age < 6 + s_i). Clamped admins no longer
   inform `s` from the likelihood side, weakening the data's pull. The
   (s, δ) ridge — flat by construction in the model — then has nothing
   anchoring it.

2. **Multimodality only appears at I=500.** Pilot fits at I=50 landed all
   chains in one mode by luck (Bayesian shrinkage was strong enough to
   smooth out the geometry). At I=500 we found **three modes**: σ_s ≈ 55
   (chain 1, all variance in s_i), σ_s ≈ 5 (chains 2+4, the pilot mode),
   σ_s ≈ 2.8 with ρ_xi_zeta = 0.999 (chain 3, LKJ boundary).

3. **The `fmax` floor is the structural problem.** Even with `s = 0.5`
   (historical), adding signed-`s_i ≥ 0` with σ_s ≈ 5 made kids with large
   positive s_i ride the floor individually. The May-16 si_reparam fit at
   `s = 0.5` already showed Rhat 1.28 on δ, ess 12 — the symptom was real
   even before the `s = 6` mistake.

**The fix.** Pin `s` tight at 0 globally, pin `s_i` off (`sigma_s_prior_sd =
0.001`). The (s, δ) ridge disappears by construction; the floor never
activates because `t − 0` is always positive for admins ≥ 12 mo.

**Affected files.**

- [`model/R/config.R`](../model/R/config.R) — `s_prior_mean = 0,
  s_prior_sd = 0.001` in `DEFAULT_PRIORS`.
- [`model/R/helpers.R`](../model/R/helpers.R) — `no_freq_slopes` variant
  comment updated to reflect the new architecture.

**Open thread.** The si_signed family fits in §21–22 are now archival.
Their elpd_loo numbers don't carry over to the new regime: the (s, δ)
freedom they had isn't in M_best anymore, so a direct comparison would
mix model classes.

---

## 🟢 24. EN M_best refit at I=500 (the multimodality discovery)

**Data.** EN longitudinal Wordbank WS, stratified to I=500 kids, J=671
items, A=1829 admins, N=1.1M obs.

**First attempt (broken).** `long_no_freq_slopes_si_signed` at I=500 with
the inherited `s_prior_mean=6, s_prior_sd=0.05` priors and σ_s free. All
four chains finished, but:

- 100% of transitions hit `max_treedepth=10`.
- σ_α Rhat = 2.94, ess = 5.
- σ_s posterior was multimodal: mean 16.9, 95% CI [2.70, 58.2]; per-chain
  medians 54.7 / 5.0 / 2.8 / 5.0.

This was the failure that motivated §23.

**Second attempt (M_best, headline).** `long_no_freq_slopes` with the
s-excised regime, cold start, 4 chains × iter=2000 warmup=1000 ×
threads=16, adapt_delta=0.95, max_treedepth=10.

| param | mean | 95% CI | ess_bulk | Rhat |
|---|---:|---:|---:|---:|
| σ_α | 1.81 | [1.69, 1.94] | 109 | 1.02 |
| σ_xi | 1.89 | [1.78, 2.01] | 109 | 1.02 |
| σ_ζ | 3.81 | [3.58, 4.07] | 218 | 1.01 |
| ρ(ξ,ζ) | **−0.09** | [−0.18, −0.01] | 154 | 1.03 |
| π_α | 0.920 | [0.910, 0.930] | 109 | 1.02 |
| δ | **10.31** | [10.25, 10.37] | 3948 | 1.000 |
| s | 0.0008 | [0.00002, 0.0023] | 3458 | 1.000 |

All chains finished in lockstep (~106 min each, 3% spread) — strong
single-mode evidence. Compare to the prior (broken) si_signed run where
chain 1 finished in 98 min and chains 2-4 took 9+ hr (the "racer vs
laggard" pattern that indicated multimodality).

**δ shifted up vs the s=0.5 historical era** (was 9.36 → now 10.31), as
expected: with `s` pinned at 0 instead of 0.5, the youngest admins'
log-age term is slightly larger in magnitude, and δ has to push harder
to fit the same data. The 10% bump in δ matches the analytical prediction.

**ρ(ξ, ζ) flipped sign** (was +0.35 historical → now −0.09). The bivariate
prior structure changed slightly when `s_i` was excised and `s` was
re-anchored; the new fit reports kid intercept and slope as essentially
uncorrelated.

**Disk-full crash during save_object.** First time the fit landed,
cmdstanr's `read_cmdstan_csv` (called internally by `save_object`) tried
to materialize the 35 GB `log_lik` array, fread spilled to /tmp which
filled, and the R session died before the .rds could be written. **No
output saved.** Lessons (added to skill files):

- cmdstanr's default per-session tempdir gets auto-cleaned on R exit, so
  if `save_object` crashes, the raw draws are also gone.
- Patched [`fit_variant_cmdstanr`](../model/R/helpers.R) to set
  `output_dir = fits/csvs_<tag>/` (persistent) and wrap `save_object` in
  `try()`. Now if it crashes, raw CSVs survive and
  [`sherlock/recover_from_csvs.R`](../sherlock/recover_from_csvs.R) can
  extract scalars + delta_j without the full materialize.
- Freed ~80 GB of disk on the VM (EN CSVs + bad-prior fits) before
  relaunch.

**Tag:** `long_no_freq_slopes` (this overwrites the May-3-era same-named
fit; the I=200 architecture-demo refit goes to
`long_no_freq_slopes_english_I200`).

---

## 🟢 25. NO M_best at I=500 (cold + warm-start)

**Bundle upgrade.** Earlier NO bundle was I=200 J=197 N=310K, sized for
pilots. Wordbank Norwegian has 1676 kids with ≥2 admins available
(median 6 admins/kid — more than EN's 3.7). Regenerated bundle at
I=500 J=674: **A=3222, N=2.06M** — ~1.85× more observations than EN at
the same I. See
[`fits/long_subset_data_nor.rds`](../fits/long_subset_data_nor.rds).

**Cold-start fit.** Same config as EN (iter=2000, warmup=1000, threads=16,
adapt_delta=0.95). All chains finished in lockstep (8305-8754 sec, 5%
spread). But σ_α / σ_xi / π_α showed Rhat = 1.23 and ess = 12 — much
worse than EN's Rhat 1.02 / ess 109. δ and σ_ζ were fine (Rhat 1.005).

**Warm-start refit.** Used STAN_INIT_FROM with posterior medians from the
cold-start as initial values (jittered per chain). Same iter/warmup;
no change to treedepth or adapt_delta. Result:

| param | mean | 95% CI | ess | Rhat |
|---|---:|---:|---:|---:|
| σ_α | 2.05 | [1.93, 2.18] | 29 | **1.11** |
| σ_ζ | 4.79 | [4.49, 5.13] | 102 | 1.03 |
| ρ(ξ,ζ) | −0.13 | [−0.21, −0.04] | 89 | 1.03 |
| π_α | **0.936** | [0.929, 0.944] | 29 | 1.11 |
| δ | **11.47** | [11.40, 11.53] | 1380 | 1.005 |

σ_α / π_α cluster has Rhat 1.11 (improved from 1.23 cold-start) but still
above the 1.05 threshold. **Headline numbers stable** (π_α = 0.94 is the
key replication finding); only the posterior width on σ_α / π_α is
imprecise. Open question: more iterations might resolve it, or it's a
slow-mixing ridge particular to NO's denser longitudinal data.

**Same disk-full incident pattern.** First save_object attempt crashed.
Recovered scalars and delta_j from persistent CSVs via
`recover_from_csvs.R`. Second attempt with disk freed (deleted bad-prior
.rds + old si_signed remnants) succeeded.

**Tag:** `long_no_freq_slopes_norwegian`.

---

## 🟢 26. IO + Peekbank refits with the cleaned priors

**Motivation.** The earlier IO/proc fits (May 3–4 era) had been done with
`s_prior_mean=0.5, s_prior_sd=0.05` (the historical pin near zero —
*not* corrupted by the s=6 disaster). They were substantively correct
but used the pre-cleanup Stan code (β_c, λ_j still in the linear
predictor, even if pinned off). For consistency with the new headline
regime, refit all six.

**Variant list.**

1. `io_no_freq_slopes` — BabyView (I=20 head-mounted video data)
2. `io_no_freq_slopes_seedlings` — SEEDLingS LENA (I=44)
3. `io_comp_no_freq_slopes_seedlings` — + comprehension channel
4. `io_std_no_freq_slopes_seedlings` — + standardized test channel
5. `io_comp_std_no_freq_slopes_seedlings` — + both
6. `long_proc_no_freq_slopes` — Stanford-linked Peekbank LWL (I=62)

**Queued via [`gcp/queue_io_proc.sh`](../gcp/queue_io_proc.sh)** after the
EN fit landed. Each fit ~30-60 min (small bundles).

**Results — all clean:**

| variant | δ | σ_α | π_α | max Rhat | min ess |
|---|---:|---:|---:|---:|---:|
| io_no_freq_slopes (BabyView) | 6.42 [5.98, 6.87] | 1.13 [0.82, 1.56] | 0.83 [0.66, 0.93] | 1.007 | 844 |
| io_no_freq_slopes_seedlings | 8.19 [7.82, 8.57] | 1.39 [1.12, 1.74] | 0.94 [0.89, 0.97] | 1.003 | 710 |
| io_comp_no_freq_slopes_seedlings | 8.01 | 1.39 | 0.938 | 1.012 | 397 |
| io_std_no_freq_slopes_seedlings | 3.97 | 1.39 | 0.885 | 1.008 | 661 |
| io_comp_std_no_freq_slopes_seedlings | 8.01 | 1.39 | 0.939 | 1.005 | 510 |
| long_proc_no_freq_slopes (Peekbank) | 2.51 | 1.67 | 0.904 | 1.002 | 1506 |

**π_α landing 0.83–0.94 across all six datasets** (different N, different
input channels, different ages). The "efficiency dominates kid-level
intercept variance" finding is robust.

**Peekbank-specific extras.** Extracted γ_rt, μ_rt, μ_rtslope, σ_rtslope,
σ_lwl, and the three pairwise correlations (ρ_alpha_zeta,
ρ_alpha_rtslope, ρ_zeta_rtslope) into
`long_proc_no_freq_slopes.draws_full.rds`. Headline:
- γ_rt = **0.084** [0.044, 0.125] — strictly positive, confirming the
  level-coupling between log_α and log_rt-intercept.
- ρ(ζ, rtslope) = **+0.05** [−0.26, +0.34] — vocab acceleration and RT
  maturation are uncoupled.
- **Reframing**: "one clock for the level, separate clocks for the
  maturation rates." Replaces the earlier "two clocks, not one" framing
  which conflated levels and slopes.

---

## 🟢 27. Cross-sectional empirical via `bundle$df` + GAMLSS BEINF

**Two methodology bugs found and fixed** when building the new quantile
plots.

**Bug A: long_items.rds has no admin_id.** Multiple admins for the same
kid at the same age appeared as repeated rows (no way to tell them apart).
`group_by(child_id, age) %>% summarise(vocab = sum(produces))` then
*summed produces across admins*. Some kids had `n_items` = 2007 = 3 × 669
items — three separate admins collapsed. Vocab counts inflated beyond J.

**Fix.** Use the bundle's `$df` slot instead. It carries `admin_key`
which uniquely identifies an admin. Helper
[`build_xsec_empirical`](../model/R/empirical_xsec_helper.R) groups by
`(child_id, age, admin_key)` and then randomly samples one admin per
child for the cross-sectional reduction. Bonus discovery: ~3% of bundle
df rows are exact duplicates (same admin_key, same item, same produces);
`distinct(admin_key, item, .keep_all=TRUE)` handles that.

**Bug B: `quantregGrowth::gcrq()` silently breaks in a wrapper.** The
package's `ps(age, lambda = X)` term reaches back into the calling
environment to resolve `lambda`. When called from inside a function, this
produced "attempt to set an attribute on NULL" deep in gcrq internals
that `tryCatch` swallowed (returning NULL silently). Spent ~30 min
chasing this. Switched the empirical quantile smoother to **GAMLSS**
with `BEINF` family (beta inflated at 0 and 1), `pbm(age, lambda=10000)`
for the mean and `pb(age)` for the variance, per the MB-CDI manual
demo recipe (Mike's `wordbank/demo-vocab/gamlss_demo.R`).

`pbm`/`pb` are well-behaved in wrapper functions; `centiles.pred()`
needs the data exposed via a global symbol (workaround in the helper).
With this in place, the empirical quantile fans on the side-by-side
plots are smooth and stable.

---

## 🟢 28. Plot suite for the slide deck (2026-05-23)

The slide deck now uses the following figures and tables. Each is
generated by a single script reading from `fits/summaries/` and (where
needed) the bundle. Provenance map in [`outputs/PROVENANCE.md`](PROVENANCE.md).

**Vocab-space quantile fans (model vs empirical).**
- [`m_best_quantile_I500.png`](figs/longitudinal/m_best_quantile_I500.png)
  — EN single panel (`quantile_demo_mbest_I500.R`)
- [`m_best_quantile_NO.png`](figs/longitudinal/m_best_quantile_NO.png) —
  NO single panel (`quantile_demo_mbest_NO.R`)
- [`m_best_quantile_EN_NO.png`](figs/longitudinal/m_best_quantile_EN_NO.png)
  — EN + NO side-by-side, bundle-internal x-sec
  (`quantile_demo_mbest_EN_NO.R`)
- [`m_best_quantile_EN_NO_wordbank.png`](figs/longitudinal/m_best_quantile_EN_NO_wordbank.png)
  — EN + NO side-by-side, wordbank-wide x-sec (slide 20 in the deck;
  `quantile_demo_mbest_EN_NO_wordbank.R`)
- [`m_best_quantile_io_proc.png`](figs/longitudinal/m_best_quantile_io_proc.png)
  — BabyView | SEEDLingS | Peekbank 3-panel
  (`quantile_demo_io_proc.R`). LOESS empirical median in lieu of full
  fan because N is too small (20/44/62) for stable per-quantile estimates.

**Architecture build-up (4-panel comparisons).**
- [`quantile_demo_4panel_I200.png`](figs/longitudinal/quantile_demo_4panel_I200.png)
  — `pure → +α → +κ → M_best`, vocab-space, English I=200
  (`quantile_demo_4panel_I200.R`). Slide 18.
- [`theta_spaghetti_4panel_I200.png`](figs/longitudinal/theta_spaghetti_4panel_I200.png)
  — same build but in latent-θ space
  (`theta_spaghetti_4panel_I200.R`). Slide 19.

**Input/processing channel plots (IO + Peekbank).**
- [`m_best_input_quartile_io.png`](figs/longitudinal/m_best_input_quartile_io.png)
  — BabyView + SEEDLingS; kids binned by their measured log_r quartile;
  model-predicted trajectory per quartile (`quantile_demo_io_input.R`).
  Slide 29.
- [`m_best_rt_quartile_proc.png`](figs/longitudinal/m_best_rt_quartile_proc.png)
  — Stanford Peekbank; kids binned by their age-adjusted RT intercept
  quartile (`quantile_demo_proc_rt.R`). Slide 31.

**Exposure-to-learn (per-word view).**
- [`exposure_to_learn_EN.png`](figs/longitudinal/exposure_to_learn_EN.png)
  — for each word, predicted age of 50% production and cumulative
  exposures-of-that-word at that age, coloured by lexical class with
  per-class lm fits (`exposure_to_learn.R`). Slide 21.

**Parameter table.**
- [`outputs/param_table.csv`](param_table.csv) + `.md`, plus the .xlsx
  Mike maintains for the slide (`param_table.R`). Slides 22, 32.

**Supporting infrastructure.**
- [`model/R/empirical_xsec_helper.R`](../model/R/empirical_xsec_helper.R)
  — shared `build_xsec_empirical` + `fit_xsec_quantile_fan` used by all
  the model-vs-data plots.
- [`sherlock/recover_from_csvs.R`](../sherlock/recover_from_csvs.R) —
  scalar + delta_j recovery from persistent cmdstanr CSV output when
  `save_object` crashes.

---

## 🟢 29. glmer model-ladder across longitudinal CDI samples (2026-05-27)

Self-contained frequentist companion to the Bayesian M_best work, in
its own project-root directory [`glmer_ladder/`](../glmer_ladder/) (see
its [README](../glmer_ladder/README.md)). The point: make the §1
("what is the best model for vocabulary growth?") claims with plain
AIC/BIC on `glmer`, on as much longitudinal data as Wordbank has, with
*no* Bayesian machinery and *no* σ_r decomposition (those belong to the
input-share section). This was motivated by §0 of the Kachergis
re-read: M_best is exactly
`glmer(produces ~ 1 + log_age + (1 + log_age | child) + (1 | word),
family = binomial)` + the external σ_r pin (verified earlier: every
glmer point estimate inside the Stan 95% CrI; implied π_α = 0.922 vs
Stan 0.920).

**The ladder (7 nested models per language).** `log_age = log(t/a0)`,
`age_c = t − a0`, `a0` = median admin age. IRT reading:
`logit P = θ_i(t) − δ_j`, δ_j = −(word RE).

| Model | glmer formula | θ_i(t) |
|---|---|---|
| A     | `~ offset(log_age) + (1\|item)` | β₀ + log(t/a₀) (κ≡1) |
| B_log | `~ 1 + log_age + (1\|item)` | β₀ + κ·log(t/a₀) |
| B_lin | `~ 1 + age_c + (1\|item)` | β₀ + β₁·(t−a₀) |
| C_log | `+ (1\|child)` | + ξ_i (random intercept) |
| C_lin | `+ (1\|child)` | + ξ_i |
| D_log | `(1 + log_age\|child)` | β₀ + ξ_i + (κ+ζ_i)·log(t/a₀) = M_best |
| D_lin | `(1 + age_c\|child)` | β₀ + ξ_i + (β₁+ζ_i)·(t−a₀) |

Four design choices, one clean ΔAIC each: A→B (free κ vs unit
accumulator), B_lin↔B_log (exponential vs power-law growth), B→C
(per-kid intercept), C→D (per-kid slope = "is σ_ζ > 0 needed?").

**Data.** Wordbank survey ([`00_survey_languages.R`](../glmer_ladder/00_survey_languages.R)
→ [`outputs/glmer_ladder/00_language_survey.csv`](glmer_ladder/00_language_survey.csv))
kept languages with ≥100 kids with ≥2 admins (any form). WG + WS
combined at the item level, production only. 7 languages qualified:
English (American) 1840 kids, Norwegian 1676, Finnish 236, French
(Quebecois) 179, Japanese 178, Spanish (Mexican) 119, French (French)
111. Spanish-MX later dropped from figures (WS-only, narrow 17–30 mo
window → degenerate D fits, σ_slope blew to 43).

**Compute.** 49 cells (7 langs × 7 models) as a Sherlock SLURM array
([`sherlock/glmer_ladder.slurm`](../sherlock/glmer_ladder.slurm), 8
cpus × 64 GB × 12 hr). All `nAGQ=0`. D_log on the big languages is the
bottleneck: EN D_log 10.2 hr (1840 kids, 3.3M obs), NO D_log 6.8 hr
(1676 kids, 4.4M obs); A/B/C cells finish in seconds-to-minutes. glmer
is single-threaded for the IRLS but RcppEigen/BLAS auto-threads the
Cholesky updates (~4.5 cores even at `--cpus-per-task=1`); explicit
8-core pinning did not buy a real speedup, and re-fits are byte-identical
(deterministic). Per-fit RDS + per-kid BLUP CSV (`ranef_*.csv`) saved
for downstream demographic analysis.

**Findings (consistent across all comparable languages).**
- A → B_log: free κ beats the unit accumulator by a huge AIC margin
  (NO: ~940k). κ̂ lands 8.6–12.7 across languages (all ≫ 1).
- B_lin vs B_log: **log wins everywhere** (EN ΔAIC = 585 at matched df,
  N=1.1M). But the two are visually near-identical over the CDI age
  window — log(t) and t are near-affine across a ~2:1 age ratio, so the
  difference is fine-tail-structure, not gross shape. Honest framing:
  you need ~10⁶ observations to distinguish exponential from power-law
  growth over 12–30 months. Prior work using either parameterization
  wasn't making a detectable fit error, just an interpretive one
  (constant-rate accumulation vs efficiency gain).
- B → C: per-kid intercept is the single biggest rung (NO: ~2.2M AIC).
- C → D: per-kid slope adds another large chunk (NO: ~115k). σ_ζ > 0
  is real and needed in every language.

**Pipeline.** `00_survey` → `01_extract_{one,all}` →
`02_fit_one` (Sherlock array) → `03_aggregate` (ΔAIC table + figure) →
`04a_simulate` (BLUP-bootstrap predictions, slow, → `sim_cache.rds`) →
`04b_plot` (cache → figures, ~7 s; iterate here). The 04 split keeps
the 500-kid × 42-fit bootstrap out of the plot-iteration loop.

**Prediction figures** ([`outputs/figs/glmer_ladder/`](figs/glmer_ladder/)).
Per (lang, model): bootstrap 500 kids from the fit's *BLUP* distribution
(not MVN(0,Σ̂) — the unshrunken parametric draws produce extreme
intercept/slope combos absent from the data and blow up the upper
quantiles at thin-data ages), compute each kid's `Σ_j inv_logit(η_ij)`
over the main form's items, take 10/25/50/75/90 quantiles. Empirical
shown as per-kid spaghetti, both restricted to the largest form per
language (so the vocab ceiling matches; e.g. Finnish WS=111 not the
201-item WG+WS union).
- `main.png` — 4 well-powered langs (EN, NO, FR-CA, JA) × the log-only
  conceptual ladder (A→B→C→D). Main text.
- `supp_log.png` — 6 langs × log ladder.
- `mega.png` — full 6 langs × 7 models (lin + log).
- `deltaAIC.png` — bar/lollipop ΔAIC summary (`03_aggregate.R`).

Visual signatures: A under-fits; B's median hugs the population mean
(item REs + κ do most of the shape work); C gives parallel quantile
fans (baseline variance); D's fans **open with age** — the σ_ζ
signature.

---

## 🟢 30. Pooled IO κ-deflation diagnostic — the delta prior was the culprit (2026-06-02)

The pooled hierarchical IO model
([`log_irt_io_pooled.stan`](../model/stan/log_irt_io_pooled.stan), fit
via [`fit_io_pooled.R`](../model/scripts/fit_io_pooled.R)) had been
reporting `kappa_pop = 6.42` and wildly heterogeneous `kappa_study`:
**3.4** (BabyView), 6.6 (SEEDLingS), 8.2 (AM2018), 7.5 (FMW2013). That
implied either a real cross-study population difference (especially
BabyView) or a model artifact. Per-study spaghetti looked healthy with
BabyView reaching ~600 items / proportion 0.75+ by 25–30 mo, suggesting
the data wasn't the problem. A systematic diagnostic chain identified
the cause as a **misspecified default prior on `delta`**.

### The chain

1. **glmer C_log per IO dataset** ([`glmer_io_datasets.R`](../model/scripts/glmer_io_datasets.R)):
   κ ∈ [9.6, 11.7] across all four studies — homogeneous, matching
   longitudinal Wordbank's κ ≈ 11.3.
2. **glmer D_log per IO dataset** ([`glmer_io_datasets_Dlog.R`](../model/scripts/glmer_io_datasets_Dlog.R)) —
   the proper apples-to-apples comparison with random slopes — κ ∈
   [10.1, 10.85], same answer.
3. **Hypotheses ruled out**:
   - λ_j is pinned at 1 in no_freq_slopes (σ_λ = 0.001), so the
     multiplicative λ-slope decomposition isn't absorbing slope.
   - δ_j anchor relaxation
     ([`fit_io_pooled_unanchored.R`](../model/scripts/fit_io_pooled_unanchored.R)):
     refit with anchor SD = 5 for all items gave identical kappa_pop
     and kappa_study. Anchor wasn't deflating slope.
   - Input-channel age confound
     ([`plot_io_input_by_age.R`](../model/scripts/plot_io_input_by_age.R)):
     within-study age × log_r_obs correlations were all near zero
     (BabyView +0.08, SEEDLingS −0.08, AM2018 +0.06, FMW2013 NA).
     Input channel isn't pulling age-related variance into ξ.
4. **The smoking gun**: extracting per-kid Stan slope posterior
   medians (slope_i = 1 + δ + study_δ + ζ_i) and averaging across
   kids gave **9.3, 10.0, 10.6, 10.6** per study — matching glmer.
   The model was *correctly fitting* per-kid trajectories; only the
   `kappa_study` summary was misleading.
5. **The mechanism**: `zeta_i` posteriors were averaging **+2.5 to
   +5.9 per study** (systematically positive), violating their prior
   `N(0, σ_ζ)`. The model was compensating for a shrunken `delta`.
6. **The cause**: `DEFAULT_PRIORS$delta_prior = N(0, 0.5)`. This
   prior was set for the cross-sectional reference fit. At full EN
   Wordbank longitudinal (I=500, N=1.1M obs), the data overwhelms
   it (`delta = 10.3`). At pooled IO (I=183, N=404k), the data is
   strong but not strong enough — the prior wins on delta, and zeta
   absorbs the shift.

### The fix

[`fit_io_pooled_widedelta.R`](../model/scripts/fit_io_pooled_widedelta.R)
overrides to `delta_prior = N(0, 10)` (essentially uninformative over
the plausible κ range). Saved as `fits/io_pooled_widedelta.rds`. The
result:

| param           | anchored (old) | wide-delta (corrected) |
|-----------------|----------------|------------------------|
| delta           | 5.40           | **9.63**               |
| kappa_pop       | 6.42           | **10.62**              |
| kappa_study[BV] | 3.40           | **10.46**              |
| kappa_study[SS] | 6.59           | **10.50**              |
| kappa_study[AM] | 8.15           | **10.79**              |
| kappa_study[FMW]| 7.48           | **10.74**              |
| sigma_zeta      | 5.22           | **4.00**               |
| sigma_alpha     | 2.16           | 2.14                   |
| sigma_r         | 0.36           | 0.36                   |
| pi_alpha        | 0.972          | 0.972                  |
| Rhat            | 1.01–1.07      | 1.00–1.01              |

All four IO studies are now homogeneous at κ ≈ 10.5–10.8 (matching
longitudinal). σ_zeta shrunk because zeta no longer needs to carry
the +3 to +6 mean shift. The **intercept partition (`π_α`, `σ_r`,
`σ_α`) is unchanged** — that was correctly identified all along, so
the input-uptake variance partition story stands as previously
characterized.

Gamma variants were similarly refit with widened delta as
[`fit_io_pooled_gamma_widedelta.R`](../model/scripts/fit_io_pooled_gamma_widedelta.R)
→ `fits/io_pooled_gamma_widedelta_{add,mult}.rds`.

### A secondary finding: the multiplicative parameterization is fragile

The wide-delta refit produced a sharp, scientifically interesting
asymmetry between the two γ parameterizations:

|                          | γ-add (wide-delta) | γ-mult (wide-delta) |
|--------------------------|---------------------|----------------------|
| `gamma` Rhat / ess_bulk  | 1.01 / 424          | **1.05 / 68.6**      |
| `delta` Rhat / ess_bulk  | 1.01 / 609          | **1.12 / 515**       |
| `pi_alpha` Rhat          | 1.00                | **1.13**             |
| `sigma_r` ess_bulk       | 1734                | **50**               |
| `sigma_alpha` ess_bulk   | 910                 | **30**               |
| `sigma_zeta` ess_bulk    | 1405                | **31**               |
| `study_input_mean[2,4]`  | clean               | **ess 21–22, Rhat 1.12** |

The multiplicative form is `slope_i = (1 + δ + study_δ + ζ_i)·(1 + γ·log_r_dev_i)`,
so the joint effect of input on slope flows through `γ·A` (where
`A ≈ 1 + δ ≈ κ_pop`). γ and A are jointly identified ONLY via their
PRODUCT — there's a ridge in the posterior along `γ·A = constant`. The
**tight default delta prior was implicitly pinning A** (around 6.4),
which collapsed the ridge to a point and made γ_mult cleanly
identifiable. With wide delta letting A roam 9–11, the ridge stretches
and HMC can't traverse it cleanly.

The additive form `slope_i = (1 + δ + study_δ) + γ·log_r_dev_i + ζ_i`
has no such ridge — γ_add is identified directly as the slope-shift-
per-input. Its posterior is healthy under either prior choice.

**Methodological takeaway:** the additive parameterization is the
robust canonical γ form. The multiplicative form is interesting
conceptually (efficiency × input synergy / "Matthew on top of Matthew")
but is not statistically tractable without an informative prior on the
scale parameter. For analyses that depend on γ, use additive; treat
multiplicative as a sensitivity check only when paired with a
regularized prior on `delta`.

### Reproducibility note

The original `io_pooled.rds` / `io_pooled_gamma_{add,mult}.rds` are
preserved as the as-built reference. Going forward, the
`*_widedelta.rds` files are the production baseline for any analysis
or figure that depends on `delta` / `kappa_pop` / `kappa_study`. The
DEFAULT_PRIORS comment in `model/R/config.R` flags the gotcha for
future users.

---

## 🟡 31. Cross-sectional demographic decomposition: sex & maternal ed across many languages (2026-06-09, in progress)

**Motivation.** The longitudinal demographic analysis (paper §"Efficiency
and acceleration relate to separate demographic predictors") only has
maternal ed for ~3 languages (after the by-study split: Marchman + Norwegian,
~2) and sex for ~4. Wordbank has *no* further longitudinal CDI data with
maternal ed (surveyed: Danish has 6,112 kids but **0** with ≥2 admins; same
for German, Italian, etc. — the big norming samples are cross-sectional).
Strategy B: recover the efficiency-vs-acceleration split from **cross-sectional**
data via the population predictor×age interaction, expanding to ~13 languages
(maternal ed) and ~30 (sex), giving a parallel two-axis design (paper finding:
**sex → efficiency**, **maternal ed → acceleration**).

**Model.** Per language, one admin per child, item-level Rasch GLMM:
`produces ~ predictor * log(age/a0) + (1|item) + (1|child)`, `glmer`, `nAGQ=0`,
`a0 = median age`. `predictor`-main = effect on **efficiency** (intercept);
`predictor:log_age` = effect on **acceleration** (slope). Monolingual-TD =
exclude `dataset_origin_name ~ "Bilingual"` + clinical (`Edgin`, `Byers`);
**not** `is_norming` (which is study-membership, stricter than mono-TD — it
drops ~half of monolingual Norwegian). Québec French excluded (bilingual);
English (British) excluded (TEDS twins, short form).

**⚠️ Bug others should know about.** Linking item data to the sampled admin
**by `child_id` silently conflates a longitudinal child's multiple admins**
(item rows summed/duplicated across admins, all mislabeled with one age →
sumscores up to 8× form size; models include the extra rows). **Fix: link by
`data_id` (the admin).** Symptom: Norwegian scatter proportions capped at
~0.25; `max items/admin` ≫ form size. Cross-sectional-only languages (Danish)
are unaffected; heavily-longitudinal ones (Norwegian) were materially biased.

**Validation (Marchman + Norwegian, where we have both methods),** raw
log(age/a0) so coefficients match the longitudinal raw-BLUP regressions
(regress raw ξ_i, ζ_i on scale(matEd)):

| | x-sec eff | long eff | x-sec acc | long acc |
|---|--:|--:|--:|--:|
| Marchman | +0.29 | +0.26 | +0.66 | +0.41 |
| Norwegian | +0.12 | +0.08 | +0.52 | +0.22 |

**→ Efficiency channel validates (x-sec ≈ long); acceleration is
systematically inflated ~1.5–2.3× cross-sectionally** (sign/order right, CIs
overlap because longitudinal acceleration is noisy). Methodological reading:
the predictor×age interaction over-attributes to slope (cf. SI Fig 1 — per-child
slope variance isn't cleanly identified cross-sectionally). **Implication:**
the **sex → efficiency** cross-linguistic analysis is on solid ground;
maternal-ed *acceleration* loadings are directional/upper-bound.
NB: z-scoring ξ and ζ *separately* (sd ζ ≫ sd ξ) makes the longitudinal effect
look efficiency-dominant — an artifact; use raw BLUP units to compare.

**Maternal-ed results (13 languages, corrected, std log-age, per SD matEd).**
Real cross-national structure in *how* SES loads:
efficiency-only (German +0.30, Latvian +0.23, Portuguese +0.11);
acceleration-only (Norwegian +0.21, Spanish-Eur +0.14);
both (English +0.18/+0.36, Mandarin-Tw, Spanish-Arg); null (Danish, Czech,
Spanish-Mex). Anomalies to scrutinize: French-French acc = **−0.28**,
Spanish-Eur/Mex eff < 0.

**Sex result (the validated, headline channel).** Female **efficiency**
advantage in 26/31 languages; meta −0.52 [−0.60,−0.44], ~0 on acceleration.
Matches longitudinal closely: meta −0.52 (x-sec) vs −0.67 (long, k=4);
per-language Norwegian −0.86 (x-sec) vs −0.89 (long). MatEd meta: efficiency
+0.10 vs +0.14 (long); acceleration +0.28 both, but long k=2 (weak check) and
x-sec acceleration is the inflated channel. Spanish/French matEd anomalies =
skewed/ill-mapped education distributions (e.g. Spanish-Eur ~78% high-ed →
unstable efficiency), not data errors; sex (balanced) unaffected.

**Status / reproducibility.** Promoted from `/tmp` to a committed, reproducible
pipeline: [`cross_sectional_demographics/`](../cross_sectional_demographics/)
(`00_build.R` + `cross-sectional_demographics.qmd` + committed `cache/fits.rds`,
`cache/scatter.rds`; per-language frames/fits gitignored + regenerable). The
notebook produces the scatter data-checks, a cross-sectional forest+meta figure
(parallels `fig-demographics`), a combined cross-sectional+longitudinal figure
(paper candidate), and anomaly diagnostics. Wordbank pulls (not glmer) are the
bottleneck; Sherlock not used (compute nodes lack internet for `wordbankr`).

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
  [`data/peekbank/cdi_short_code_map_ws.csv`](../data/peekbank/cdi_short_code_map_ws.csv)
  and `cdi_short_code_map_wg.csv` to confirm the manual_disambig rows
  (especially for body-parts compounds, helping verbs with slashed
  forms, and place-names with `*` annotations) match what the form
  actually printed. Replace any wrong mapping in
  `manual_overrides` and rerun.


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
