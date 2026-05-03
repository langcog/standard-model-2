# Parameter-recovery test: results

**Setup.** I=250 children × J=150 words × C=3 classes → N=37,500 synthetic observations, age range 12–30 mo. Fit with 4 chains × 1500 iter (750 warm-up) in rstan. Total sampling: ~3.7 min. 1 divergence out of 3,000 (0.03%), Rhat ≤ 1.006 on all scalar parameters.

## Two parameterizations, only one works

**v1 (the original).** Both `log_r[i]` and `log_alpha[i]` as individual latents, with priors tying their SDs to `sigma_r` (fixed) and `sigma_alpha` (free). Since the likelihood only sees their sum, this is over-parameterized (2*I parameters where I are identified).

Result on the small sim: `sigma_alpha` posterior 0.24 [0.01, 0.46] for true 0.50 — **biased low and barely misses CI**. On the big sim, one chain got stuck in a bad region and never terminated (28 min and counting). Killed.

**v2 (what we use from here on).** Collapse to one latent per child, $\xi_i = \log r_i + \log \alpha_i$, with $\xi_i \sim \mathcal{N}(\mu_r, \sigma_r^2 + \sigma_\alpha^2)$. Posterior-expected `log_alpha[i]` reconstructed in `generated quantities` via the closed form for $\mathbb{E}[\log \alpha \mid \xi]$.

Mathematically equivalent, but removes the sampling pathology.

## v2 results (headline)

| Parameter | Truth | Posterior median | 95% CrI | In CrI? |
|---|---:|---:|---|---|
| **pi_alpha** (RQ4) | **0.610** | **0.624** | [0.501, 0.718] | ✓ |
| sigma_alpha | 0.50 | 0.516 | [0.401, 0.638] | ✓ |
| sigma_xi (derived) | 0.640 | 0.653 | [0.566, 0.753] | ✓ |
| s (RQ2) | 4.50 | 4.15 | [0.66, 7.88] | ✓ |
| delta (RQ3) | 0.10 | −0.16 | [−0.49, 0.20] | ✓ |
| mu_c[1] | 6.50 | 6.59 | [6.33, 6.84] | ✓ |
| mu_c[2] | 8.00 | 8.10 | [7.79, 8.44] | ✓ |
| mu_c[3] | 9.50 | 9.75 | [9.34, 10.28] | ✓ |
| tau_c[1] | 0.50 | 0.48 | [0.37, 0.63] | ✓ |
| tau_c[2] | 0.70 | 0.65 | [0.45, 0.93] | ✓ |
| tau_c[3] | 0.70 | 0.70 | [0.41, 1.12] | ✓ |

- **Word-level `psi_j`**: correlation(truth, posterior median) = 0.955; CrI coverage = 0.95 (exactly nominal).
- **Child-level `xi_i` (sum)**: correlation = 0.80.
- **Child-level `log_alpha_i` (decomposed)**: correlation = 0.60 — reasonable for cross-sectional data with only 150 items per child.

See `model/recovery_psi.png` and `model/recovery_log_alpha.png`.

## What this buys us

1. **RQ4 is answerable.** The variance-decomposition identification argument from the explainer doc works in practice: pinning $\sigma_r$ from external data lets us recover $\sigma_\alpha$ and $\pi_\alpha$ within a proper credible interval.

2. **`psi_j` recovery is excellent.** 0.955 correlation with nominal coverage means the per-word thresholds we extract from Wordbank will be meaningful absolute quantities in log-token units.

3. **`s` and `delta` are recoverable but uncertain with this sample size.** The CIs are wide (covering zero for delta, covering the prior for s). Real Wordbank (36× more children, 4.5× more items) should tighten them considerably, making RQ2 and RQ3 sharp.

## What's different from the explainer doc

**Link function.** The doc writes a probit link ($y \sim \text{Bernoulli}(\Phi(\eta))$); the code uses the canonical Rasch/IRT **logit** link ($y \sim \text{Bernoulli}(\sigma(\eta))$). These differ by a scale factor of ~1.7 on the latent predictor; `psi_j` is interpretable as log-threshold on the logit scale. I'll reconcile the doc when updating it.

**Parameterization.** The doc describes the v1 factoring of child-level latents into $\log r_i$ and $\log \alpha_i$. The v2 code collapses these into $\xi_i$ (identical likelihood, sampler-friendly) and reconstructs expected `log_alpha` in generated quantities. I'll note this in the doc.

## Next steps

1. Update `model_explainer.pdf` to reflect the logit link and the v2 parameterization.
2. Preprocess Wordbank + CHILDES data into Stan-ready format.
3. Fit v2 to a subset of Wordbank (say 500 children × 200 items) for a real-data shakedown.
4. Scale to the full CDI:WS data.
5. Planned model comparisons (variants B, C for RQ1; sensitivity on $\sigma_r$ for RQ4).
