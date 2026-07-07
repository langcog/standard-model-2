# bayes_long — Bayesian longitudinal accumulator ladder (by-dataset)

Publication fits for the acceleration story. Bayesian re-fit of the descriptive
accumulator ladder, **one fit per dataset** (matching the glmer by-study units:
thal / smith / marchman / norwegian / japanese), replacing the earlier pooled
English + Norwegian Stan fits.

## The ladder (four hard-structured models, one component per rung)

Each rung is its own minimal Stan model in [`stan/`](stan/) — no prior-SD toggles
(a near-zero variance component funnels HMC), so each fit is unambiguously the
model it claims to be.

| model | file | adds | between-child variation |
|---|---|---|---|
| **M0** | `m0_accumulator.stan` | κ = 1, item difficulties only | none (one global curve — the LLM analog / falsifiable null) |
| **M1** | `m1_acceleration.stan` | + free acceleration (κ = 1 + δ) | none |
| **M2** | `m2_efficiency.stan` | + per-child intercept ξ_i | efficiency only |
| **M3** | `m3_full.stan` | + per-child slope κ_i (correlated) | efficiency **and** acceleration (**headline**) |

Two comparisons carry the paper: **M0→M1** (κ > 1: children are not pure
accumulators — the LLM contrast) and **M2→M3** (acceleration *varies* between
children). M1 and M2 are included so each step moves one component.

## Model design (settled 2026-07)

- **Rasch core** `y ~ Bernoulli_logit(xi_i + log_H + kappa_i·log(age/a0) − delta_j)`,
  `reduce_sum` for within-chain threading (kept from the io-proc lean model).
- **Flat item-difficulty hierarchy** `delta_j ~ N(0, tau_delta)`, mean-centered
  (no lexical-class machinery — verified-equivalent for the headline params).
- **Free population intercept `mu_xi`** — no pinned `mu_r`, no `sigma_r`/input
  decomposition (this is a purely descriptive model).
- **Correlated (ξ, κ)** in M3 via a 2×2 Cholesky (glmer's `(log_age|child)`;
  earlier ρ ≈ −0.14 EN / −0.38 NO; ρ shapes the fan).
- `log_H` kept as an interpretive offset (accumulator units; used downstream).
- Stripped vs the lean base: LWL processing, input measurement, per-study index,
  2PL discrimination, frequency offset, developmental onset s/s_i.

## Data notes (from the by-dataset spaghetti check)

- **Longitudinal density is low** (median 2–3 admins/kid). M3's σ_b (acceleration
  variance) is the rung most exposed to this — watch its posterior / rhat per
  dataset.
- **Norwegian = full Wordbank longitudinal** (1,676 kids, sparse), not the dense
  200-kid subset the earlier pooled fit used.
- **Thal WG/WS**: combines WG (396 items, ~12–16 mo) and WS (680 items, 16–29 mo);
  the proportion "dip" at the transition is a scale artifact (the model handles
  multi-form via per-item difficulties). Two floating items — **"in" and
  "inside"** (a WS form-version quirk, on 641 of the WS admins) — dropped to the
  canonical 680-item WS set in bundle prep.

## Pipeline

- `00_prepare_bundles.R` — by-dataset Stan bundles (TODO)
- `01_fit.R` — fit one (dataset, model) → summaries / draws / loo (TODO)
- `02_compare.R` — LOO ladder table (TODO)
- `03_fan.R` — by-dataset BLUP-bootstrap fan (Fig 1 bottom panel) (TODO)
