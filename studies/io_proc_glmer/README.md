# io-proc glmer ladder — the simple-measurement-model benchmark

A deliberately simple **glmer benchmark** for the input and processing channels,
built to sanity-check the Bayesian io-proc measurement model. The goal: the
Bayesian model should *reproduce/improve on* these estimates, not disagree with
them. (Origin: the Bayesian D′3 fit reported input explaining almost no
acceleration variance; MCF asked for a model-independent check.)

## Channels (one number per child, z-scored within study)
- **input_z** — mean observed log input rate per child (LENA / head-cam), from the
  bundle's `log_input_obs`.
- **proc_z** — from a simple RT measurement model `lmer(lwl_log_rt ~ log_age + (1|dataset))`;
  take each child's **mean residual**, flip sign, z-score → higher = **faster** processing.

cor(input_z, proc_z) = **0.24** (142 both-channel kids) → weak, so the dissociation
below isn't a collinearity artifact.

## Models (`produces ~ <chan>*log_age + (log_age|child) + (1|item) + (1|dataset)`, a0=21)
Fit on the **common** sample (142 kids with both channels) for the clean nested
comparison, plus **standalone** anchors on each channel's full sample.

| spec (sample) | input→level | input→accel | proc→level | proc→accel |
|---|--:|--:|--:|--:|
| input · full (193) | 0.351 | **0.903** | — | — |
| input · common (142) | 0.285 | **0.918** | — | — |
| proc · full (326) | — | — | **0.580** | 0.182 |
| proc · common (142) | — | — | **0.662** | 0.474 |
| **both · common (142)** | 0.13 (n.s.) | **0.845** | **0.625** | 0.25 (n.s.) |

### ⭐ Double dissociation (robust to joint fit + low collinearity)
- **Input → ACCELERATION** (slope 0.85–0.92; level ~0 once adjusted).
- **Processing → LEVEL** (intercept ~0.6; slope n.s.).

## Benchmark vs the Bayesian model
SM2 D′3 (approx θ≈logit scale, λ̄=1, input only): input→level ≈ 0.36 (forced by the
io identity `ξ = μ_r + 1·d_i`), input→accel ≈ 0.25 (γ_in·σ_r, wide CI). So the
Bayesian model **under-credits input→acceleration** (0.25 vs 0.85) and **over-credits
input→level** (0.36 vs the adjusted 0.13). See `figs/io_proc_glmer_coefs_vs_sm2.png`.
This is *not* mixing (r̂ 1.03) and *not* processing competition (γ_in flat across the
D′0–D′3 ladder) — it's the shared-latent / coeff-1-on-ξ structure. See
`../../journal/notes/glmer_ladder_benchmark.md` and `mm_run_state.md`.

**Modeling implication:** the data want input on the *slope* (a compounding /
Matthew-effect accumulator: input sets the growth exponent κ) and the *level* belongs
to processing — not the simple "input = constant rate multiplier ⇒ input in the level"
that coeff-1-on-ξ encodes.

## Reproduce
```r
Rscript studies/io_proc_glmer/fit_ladder.R    # fits 5 models -> cache/ + results/glmer_ladder_coefs.csv (~40 min)
Rscript studies/io_proc_glmer/plot_ladder.R   # -> figs/io_proc_glmer_coefs{,_vs_sm2}.png
```
`cache/` (fitted glmer objects) and `results/fit_log.txt` are gitignored.
