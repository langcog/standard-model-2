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
