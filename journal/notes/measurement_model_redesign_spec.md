# Spec: input + processing as measurement models with literature priors

**Status:** spec only (2026-06-13). No fitting yet. For MCF review.

## The problem (verified in code)

The joint `proc_dp` model treats its two observed channels **asymmetrically**:

- **RT (processing) — done right.** `lwl_log_rt ~ N(tau_s + rt0_i + (psi_s+rt1_i)·log_age, sigma_lwl)`; `sigma_rt0` (between-child RT SD) is a **free parameter estimated from data**, noise-corrected by `sigma_lwl`.
- **Input — done wrong.** The bundle computes the observed noise-corrected between-child input SD (`sd_true` ≈ 0.45/study) and then **divides it out**: `z_lena = (lmean − mean)/sd_true`. The model forces `z_r ~ std_normal()` and sets the input contribution to `sigma_r * z_r` with **`sigma_r` pinned as data**. So `var_input_xi = sigma_r²` (pinned). We use the observed input only for *relative positions* (who is high/low), and **throw away the scale** — the actual between-child input variance we measured.

This asymmetry is an artifact: the σ_r pin was inherited from the **io-imputed** EN/NO model (no observed input → must pin). The joint model copied it even though it has ~190 kids with real input.

**Consequence:** the input share is a literature-pin (`σ_r²/σ_ξ²`), not a data estimate — and its "uncertainty" is the literature band [0.36, 0.58], which no amount of data fixes. This is the imprecision MCF flagged, and it mis-describes what we have.

## The fix: input as a measurement model, symmetric with RT

Make input mirror RT exactly. Latent per-child true log-input deviation `d_i`; observed per-recording log-input measures it with per-study noise:

```stan
// parameters
real<lower=0> sigma_r;                 // between-child input SD  (ESTIMATED)
vector[I]     d_raw;                    // per-child input deviation (raw)
vector[S]     mu_r_s;                   // per-study input level (like tau_s for RT)
vector<lower=0>[S] sigma_meas;          // per-study input measurement noise (headcam vs LENA)
// transformed
vector[I] d = sigma_r * (d_raw - mean(d_raw));     // centered input deviation
// model  (V_obs raw per-recording log-input observations)
d_raw     ~ std_normal();
sigma_r   ~ <input-rate literature prior>;         // anchor, not a hard pin
log_input_obs ~ normal(mu_r_s[study_of_rec] + d[rec_to_child], sigma_meas[study_of_rec]);
// xi uses d_i directly:
xi = mu_r + d + beta_xi*rt0 + log_alpha;           // var_input_xi = sigma_r^2 (now DATA-driven)
```

Structural parallel:

| | between-child scale | study level | measurement noise |
|---|---|---|---|
| **RT** | `sigma_rt0` (est.) | `tau_s` | `sigma_lwl` |
| **input (new)** | `sigma_r` (est.) | `mu_r_s` | `sigma_meas_s` |

**This also dissolves the whole "which σ_r?" debate** — instead of *picking* a sample-matched σ_r from the literature, we **estimate it from our samples' own observed input variance** (which IS the sample-matched quantity), with the literature as a prior. The apples-to-apples work (entry 36) becomes the *prior*, the data updates it.

## Priors

**Input `sigma_r`** — informative from the input-rate literature (entry 36, channel-matched US English): center ≈ **0.44**, e.g. `sigma_r ~ normal(0.44, 0.10)` truncated >0 (or lognormal matching the [0.36, 0.58] spread). Data can move it; prior keeps it identified given measurement noise.

**RT — frank_etal_2026 priors** (re-fit `log_rt ~ log_age + (1+log_age|child)` on the full Peekbank RT data, 327 kids / 281 longitudinal, age-centered 18 mo):

| param | frank2026 value | current joint D′3 | prior to set |
|---|---|---|---|
| pop. log-RT @18mo (`tau`) | 6.90 (≈991 ms) | — | `normal(6.9, 0.3)` |
| pop. slope `psi` (dlogRT/dlogAge) | **−0.53** | — | `normal(−0.53, 0.15)` |
| between-child level **`sigma_rt0`** | **0.174** | 0.142 (shrunk) | `normal(0.174, 0.05)` |
| between-child slope `sigma_rt1` | 0.261 | — | `normal(0.26, 0.08)` |
| residual `sigma_lwl` | 0.202 | 0.203 ✓ | `normal(0.20, 0.05)` |

(σ_lwl already matches → the measurement layer is consistent; the gain is constraining σ_rt0 up from the shrunk 0.142, which **raises and tightens the processing share**.)

## Bundle changes

`prepare_joint_io_proc_bundle.R`: stop standardizing. Pass **raw per-recording log-input** (`io_pooled`'s `log_r_obs`) + `rec_to_child` + `study_of_recording`, for all input-having kids (AM2018/FMW2013 16+18-mo LENA reads; BabyView ~261 head-cam clips/kid; SEEDLingS ~12 daylong/kid). The model estimates `sigma_meas_s` from the within-child replicate spread and `sigma_r` from the between-child spread — no pre-computed `sd_true`, no `z_lena`.

## Open design questions (for MCF)

1. **σ_r: one pooled value or per-study?** Pooled is the cleaner population quantity and matches the RT `sigma_rt0` treatment; per-study would let BabyView/SEEDLingS differ but fragments the estimate. **Lean: pooled**, with `mu_r_s` (level) per study absorbing study differences.
2. **σ_meas per study: estimate or pin?** Head-cam (≈260 clips/kid, low noise) vs LENA daylong (few/kid, higher) genuinely differ → **per-study, estimated** (the replicates support it). AM2018/FMW2013 have only the 2 (16/18-mo) reads — thin, may need a weak prior.
3. **Input as a trait or with an age slope?** Current = trait (`d_i`, no age term). RT has a developmental slope; input rate may also drift with age. **Lean: keep trait for now**, flag as an extension.
4. **Cross-language kids** (if Weisleder Spanish data arrives from the lit review): different `mu_r_s` + CDI, same `sigma_r`? Defer until data exists.

## Predicted effect

- **Input share → data-driven.** σ_r estimated from ~190 kids' observed between-child input variance (noise-corrected), prior-anchored. Point ≈ unchanged (we chose 0.44 well), but the **CI replaces the literature band with a data+prior posterior — likely tighter**, and the *interpretation* flips to "estimated from observed home recordings."
- **Processing share → up + tighter.** frank2026 σ_rt0 prior (0.174) corrects the shrunk 0.142 → processing share rises from ~3.1% toward ~4–5% and narrows.
- **Separation (β_xi, γ_in) → ~unchanged** — still bounded by the 97 both-channel kids (the lit review / TotLot kids are the only lever there).

## Implementation steps (when greenlit)

1. New model `log_irt_long_proc_dp_joint_mm.stan`: input block → measurement model (above); add the frank2026 RT priors as data-passed hyperparameters.
2. `prepare_joint_io_proc_bundle.R` (or a `_mm` variant): pass raw log-input + indices; add the RT + input prior hyperparameters.
3. Fit D′3 on Sherlock; compare σ_r posterior vs the 0.44 pin, σ_rt0 vs 0.142, and the share CIs vs the current fit.
4. If σ_r posterior ≈ 0.44 with a tighter CI → the figure's Panel A σ_r-band framing can become a single data-driven estimate.
