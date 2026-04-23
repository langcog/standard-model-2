# Experiments log

A running record of fits, findings, and backlog for the standard-model
project. For the durable model specification, see
[`model_explainer.pdf`](model_explainer.pdf).

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

## 🟢 4. σ_r sensitivity sweep (2PL, σ_r ∈ {0.3, 0.53, 0.8, 1.2})

**Setup:** Refit the 2PL variant on the same 500 × 200 subsample under
four values of the externally-pinned σ_r. 2 chains × 1500 iter each.

| σ_r | σ_α | π_α | σ_xi² |
|---:|---:|---:|---:|
| 0.30 | 2.18 | **0.981** | 4.84 |
| 0.53 | 2.13 | **0.941** | 4.82 |
| 0.80 | 2.05 | **0.868** | 4.84 |
| 1.20 | 1.83 | **0.699** | 4.79 |

**Key finding:** total child-level variance σ_xi² = 4.8 is tightly
identified by the data. Its **decomposition into input vs. efficiency
is entirely determined by σ_r**. π_α ranges from 0.70 to 0.98 across
plausible σ_r values. Even at the widest σ_r = 1.2, input explains at
most 30 % of child-level variance — a robust lower bound on the
efficiency share.

**Paper-ready phrasing:**

> Depending on how much input rates vary across the population the
> CDI sample represents (σ_r ∈ [0.3, 1.2]), individual differences in
> learning efficiency account for 70–98 % of between-child variance in
> vocabulary size at 16–30 months. The total child-level variance is
> tightly identified; its decomposition is prior-bound.

**Artifacts:** `model/fits/sensitivity_sigma_r_2pl.rds`,
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

## Backlog (⚪)

### Data / robustness
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
- **Class-specific frequency slopes β_c.** Extension to RQ1 — tests
  whether frequency matters less for function words.

### Instrumentation
- **Observable.js app** for interactive exploration of the fitted
  model once the posteriors stabilize. Sliders for σ_r, σ_α, σ_ζ, s,
  δ → live growth curves and distributions.
