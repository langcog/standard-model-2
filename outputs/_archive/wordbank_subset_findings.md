# First Wordbank fit — findings and diagnostic plan

## Fit details

- **Data**: 500 children × 200 items, stratified across age bins and lexical classes → 100,000 obs.
- **Fit**: `log_irt_v2.stan`, 4 chains × 1500 iter (750 warmup). 37 min wall clock. 0 divergences, 0 treedepth saturations.
- **External prior calibration** (from `standard_model/scripts/data/hourly_tokens_Sperry_HartRisley.csv`): `adult_child_tokens_hr` n=42 → μ_r = 7.34, σ_r = 0.53.

## Good news

### RQ1 is clean and interpretable

Posterior correlation between per-word ψ_j and log p_j: **r = 0.675, R² = 0.456**. 95% CrI on correlation: [0.667, 0.681].

![psi vs log_p](../model/wordbank_psi_vs_logp.png)

The positive correlation is informative. Since η = log p_j + ... − ψ_j, a positive r(ψ, log p) means the model is *compensating* for high-frequency words by giving them high thresholds. Translation: the pure-accumulator assumption (one token = one unit of learning) systematically *over-credits* frequent words. The model corrects for this by raising their ψ_j.

This matches Braginsky (2019) directly: frequency is only modestly related to acquisition difficulty.

By class:
| Class | r(ψ, log p) |
|---|---:|
| adjectives | 0.59 |
| function_words | 0.51 |
| **nouns** | **−0.25** |
| other | 0.50 |
| verbs | 0.51 |

Nouns are the exception — for nouns, higher frequency actually does correspond to lower ψ (easier). That's the "pure accumulator" regime. Other classes show the reverse.

### Class means are ordered correctly

Posterior medians of μ_c (log-thresholds by class):

| Class | μ_median | 95% CrI |
|---|---:|---|
| nouns | 1.58 | [1.08, 2.11] |
| other | 2.05 | [1.26, 2.92] |
| verbs | 3.49 | [2.96, 4.10] |
| adjectives | 3.70 | [3.16, 4.26] |
| function_words | 7.09 | [6.30, 7.81] |

Noun bias falls out cleanly: nouns have the lowest thresholds. Function words are 5.5 log-units harder (≈250× more tokens needed).

## Problems

### 1. Rhat warnings (chains not mixed)

Largest Rhat = 1.10 on μ_c parameters. n_eff ~70–200 across scalar parameters. **The posterior is multi-modal or very correlated.** Current fit is not trustworthy for publication.

### 2. Scalar parameters are extreme

| Parameter | Posterior median | Expected | Interpretation |
|---|---:|---|---|
| σ_α | 1.98 | ~0.5–1 | Ratio of e^(2·1.98) ≈ 54× between +1 and −1 σ kids. Too extreme. |
| π_α | **0.93** | 0.3–0.7 | Model claims 93% of child-level variance is efficiency. |
| s | 12.2 mo | 2–6 mo | Hits near upper prior bound (3.9 σ above prior mean). |
| δ | 2.31 | 0–0.5 | Age exponent = 3.3. A month at age 30 is 3.3× "more effective" than at 20. |

### 3. What's driving these values

Real Wordbank CDI growth is steep: typical words go from ~5% to ~95% over a 6-month age window (on the logit scale, a change of ~6 units). Pure linear accumulation with δ=0 can only generate log(a_max/a_min) ≈ 0.6 units over the same window. The model has three levers to close the gap:

- **Large δ**: amplify the age contribution.
- **Large σ_α**: let children transition at different ages (high-α kids cross threshold early, low-α kids late). Aggregate curve steepens.
- **Late s**: compress the log-age scale so the ratio a/(a_0) varies more across ages.

The posterior cranks all three. They trade off against each other — which is why Rhat is bad.

This is the **fundamental model-misspecification issue** Mitchell & McMurray flagged about pure accumulator models: without an accelerative mechanism (leveraged learning, changing rate, nonlinear difficulty, etc.), linear time just doesn't generate sharp enough curves.

## Diagnostic plan

I want to run three focused fits to pin down what's happening before scaling up. Each is a small perturbation of the current model:

1. **Fix δ = 0 and s = 2 mo** — force "pure accumulator + early start". Expect very poor fit, but shows how much variance the accelerative machinery is absorbing. Quick diagnostic.

2. **Fix δ = 0, s free** — does s alone explain the steepness by compressing the log-scale?

3. **Fix s = 2 mo, δ free** — does δ alone explain the steepness?

If any of these fits as well as the all-free model (via LOO-CV), we know that extreme parameter is not needed. If none do, the data really require the full flexibility, and we can move on.

Additionally:

- Rerun the current all-free model with **2× more iterations** and more aggressive **`adapt_delta = 0.98`** — should fix the Rhat if the posterior is just slow-mixing.
- Look at pairs plots between (σ_α, δ, s) to confirm the trade-off.

## Decision needed

Two questions for you:

**Q1.** Am I right to be alarmed by π_α = 0.93? Or is it plausible to you that most of the child-level variance in CDI vocabulary at this age is individual-difference, not input-driven? The extreme σ_α value is what makes me suspicious, but it's not impossible a priori.

**Q2.** Of the three diagnostic fits above, do you want me to run all three, or skip some? My default: run all three, plus a longer rerun of the baseline, in parallel. Maybe 1–2 hrs wall clock total.

Once we understand what's driving the extreme posterior, we'll know whether to scale to full Wordbank as-is or fix the model first.
