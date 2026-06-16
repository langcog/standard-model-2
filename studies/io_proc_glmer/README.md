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

## Benchmark vs the Bayesian model — RESOLVED (the cause is the γ_in prior)
The SM2 D′3 (canonical) looked like it **under-credited input→acceleration** (0.25 vs the
glmer's 0.85; share 0.8% vs ~4%). We traced the cause through three eliminations:
- **Not mixing** (r̂ 1.03). **Not processing competition** (γ_in flat ~0.7 across the D′0–D′3
  ladder). **Not the coeff-1-on-ξ efficiency identity** — the G0–G3 glmer morph showed pinning
  efficiency to 0.358 does *not* suppress acceleration (0.845→0.964), and latent-vs-observed is
  moot (cor 0.997). (So the io identity is **fine**; no re-spec needed.)
- **It's the N(0,1) prior on γ_in.** Widening to N(0,5) (Sherlock job 29719815): γ_in 0.70→1.69,
  input→accel 0.25→0.60, share 0.8%→3.0% — ~reconciling with the glmer (~4%). Also λ̄=1.001, so
  θ≈logit and the overlay scale is exact. The N(0,1) prior was unintentionally shrinking a
  **weakly-identified** parameter toward 0.

**Why weakly identified:** σ_ζ ≈ 4.3 dominates the slope variance, so *all three* slope-channel
coefs are barely pinned (γ_in [-0.49,1.89], β_k0 [-1.78,1.47], β_k1 [-1.98,1.30]) and prior-sensitive
— not just input. σ_ζ's own prior is innocent (half-normal(0,1), posterior 4.3, data-dominated).
Figure `figs/io_proc_glmer_coefs_vs_sm2.png` shows both priors (filled vs open triangles); notes in
`../../journal/notes/glmer_ladder_benchmark.md` and `mm_run_state.md`.

**Takeaway:** the double dissociation (input→acceleration, processing→level) **is recoverable inside
the existing io model** with an honest weakly-informative prior on the slope coefficients. The fix is
the prior, not the architecture. Pending: a systematic prior-sensitivity sweep on the slope-channel
coefs (+ σ_r's literature prior).

## Reproduce
```r
Rscript studies/io_proc_glmer/fit_ladder.R    # fits 5 models -> cache/ + results/glmer_ladder_coefs.csv (~40 min)
Rscript studies/io_proc_glmer/plot_ladder.R   # -> figs/io_proc_glmer_coefs{,_vs_sm2}.png
```
`cache/` (fitted glmer objects) and `results/fit_log.txt` are gitignored.

## Figures
- `figs/io_proc_glmer_coefs.png` — glmer specs, level | acceleration facets.
- `figs/io_proc_glmer_coefs_vs_sm2.png` — + SM2 D′3 overlay (N(0,1) filled, N(0,5) open triangles).
  **SM2 rows are input-only on purpose:** `proc_z` reliability ≈ **0.16** (per-child RT mean over
  ~8 trials is ~84% noise; child-RE RT model σ_child=0.069 vs σ_resid=0.385), so the SM2↔glmer
  processing scale-bridge is unstable (disattenuation implies implausible ~18/log-RT effects). The
  glmer processing rows stand on their own; we just don't overlay SM2 processing. (This low
  reliability also explains the wide processing→acceleration CIs.)
- `figs/io_proc_glmer_coefs_rescaled.png` — **ONE axis**: acceleration coefs rescaled to
  level-equivalent θ units (× log(t_ref/a0), t_ref=30, a0=21 ⇒ ×0.357), so input→acceleration and
  processing→level are directly comparable as "contribution to log-odds of production by 30 mo."
  (θ≈logit since λ̄=1.) The controlled rows (both-adjusted, SM2) are the clean comparisons; input
  rows' level effect is inflated by processing leakage when proc isn't in the model.
