# Priors across the io-imputed and io-proc models (SI draft)

Draft of the priors table for the SI. Values are the *effective* priors used in
the reported fits (DEFAULT_PRIORS in `model/R/config.R` + variant overrides in
`variant_hyperpriors()` + bundle-pinned scalars + the D′-ladder rung toggles).
HalfN = half-normal (the parameter has `<lower=0>`).

| parameter | io-imputed (D, `long_no_freq_slopes`) | io-proc (joint, `proc_dp`) | role / source |
|---|---|---|---|
| **σ_r** (between-child input SD) | **pinned 0.53 → 0.44** | **pinned 0.44** | external input-rate literature (entry 36). *The measurement-model redesign frees this — estimate from observed input.* |
| **σ_alpha** (efficiency var.) | HalfN(0, 1) | HalfN(0, 1) | weakly informative; free |
| **σ_zeta** (acceleration var.) | HalfN(0, 1) | HalfN(0, 1) | free (longitudinal) |
| **δ** (accumulation exponent) | N(0, 0.5) | N(0, 0.5) | shrinks toward 0 unless data overwhelm (EN longitudinal does; small IO samples don't) |
| μ_c, σ_c (class mean ability) | N(8, 3) | N(8, 3) | per-lexical-class |
| **δ_j** (per-word difficulty) | hierarchical N(0, σ_δ) | same | word random effect |
| **μ_rt** (pop. log-RT @ a₀) | — | **N(6.5, 1)** ← vague | *→ frank2026: N(6.84, 0.2)* |
| **μ_rt-slope** (RT age trend) | — | **N(0, 1)** ← vague | *→ frank2026: N(−0.35, 0.1)* |
| **σ_rt0** (between-child RT level) | — | **HalfN(0, 1)** ← vague | *→ frank2026: N(0.143, 0.04)* |
| **σ_rt1** (between-child RT slope) | — | **HalfN(0, 1)** ← vague | *→ frank2026: N(0.26, 0.08)* |
| **σ_lwl** (within-child RT noise) | — | **HalfN(0, 1)** ← vague | *→ frank2026: N(0.24, 0.05)* |
| **σ_lena** (input meas. noise) | — | pinned per-study [.39, –, .39, –, .15, .31] | from within-child input replicates |
| **γ_in** (input → acceleration) | — | HalfN(0, 1) (on) | D′0+ |
| **β_xi** (rt0 → efficiency) | — | HalfN(0, 1) free / HalfN(0, .001) pinned | D′1 toggle |
| **β_k0, β_k1** (rt → acceleration) | — | HalfN(0, 1) free / HalfN(0, .001) pinned | D′2/D′3 toggles |
| s, σ_s (onset) ; σ_λ (2PL) ; β_c (freq) | pinned ~0 (.001) | pinned ~0 (.001) | off by default |

## Two observations this table makes plain

1. **The entire RT block is vague** (HalfN(0,1) / N(·,1)) — this is exactly what
   the **frank2026 priors** replace (right-hand column). The current proc/joint
   fits let the data carry the RT measurement model unaided; informative priors
   tighten it.
2. **σ_r is pinned** in both models — the **input-as-measurement-model redesign**
   (`measurement_model_redesign_spec.md`) frees it. So the SI priors table and the
   spec tell one story: the two least-justified choices are the σ_r pin and the
   vague RT priors, and both are being addressed.

(To finalize: confirm the δ_j hierarchical prior SD and the exact `variant_hyperpriors`
overrides per reported fit; cite the literature anchors for σ_r and the RT priors.)
