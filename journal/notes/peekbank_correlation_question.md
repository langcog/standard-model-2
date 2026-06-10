# Peekbank correlation: $\rho(\zeta_i, \mathrm{rtslope}_i)$

## What I need from the Peekbank-Stanford fit

The posterior correlation between per-child slope-deviation parameters across CDI and LWL readouts:

- $\zeta_i$: per-child deviation from population log-age slope on CDI ability ($\delta + \zeta_i$)
- $\mathrm{rtslope}_i$: per-child deviation from population log-age slope on LWL log-RT

Per the existing model spec, $(\log\alpha_i, \zeta_i, \mathrm{rtslope}_i) \sim \mathrm{MVN}(0, \Sigma)$ with LKJ on the correlation matrix, so this correlation is already in the posterior of the joint fit; it just hasn't been pulled out and reported as a headline.

**Specifically, please extract and report:**

1. Posterior median + 95% CrI for $\rho(\zeta_i, \mathrm{rtslope}_i)$
2. Posterior median + 95% CrI for $\rho(\log\alpha_i, \mathrm{rtslope}_i)$ (already partially captured by $\gamma_{rt}$, but the random-effect-level correlation is the cleaner number)
3. Posterior median + 95% CrI for $\rho(\log\alpha_i, \zeta_i)$ — this is the within-CDI quality-vs-growth coupling, and is informative even without LWL

If the existing fit doesn't have all three (e.g., if the MVN was reduced to a bivariate $(\xi, \zeta)$ before LWL was joined in), then please refit with the full trivariate prior so the correlations are jointly estimated.

## Why this matters

The brainstorm with the side-quest session landed on a structural distinction between two between-child axes that the SM2 parameterization can in principle separate:

- **Quality** — age-invariant property of input or learner — sits in $\sigma_\alpha$ (intercept).
- **Efficiency-growth-rate** — age-varying maturation property — sits in $\sigma_\zeta$ (slope).

Empirically $\sigma_\zeta > \sigma_\alpha$ in the longitudinal fits, suggesting growth-rate variation dominates baseline-quality variation across kids. The next step is to ask whether this growth-rate variation is *one underlying process* (efficiency maturation, expressible in any sufficiently sensitive readout) or *two independent maturation processes* (lexical-knowledge growth and processing-speed growth, on different clocks). Cross-readout correlation is the discriminating quantity.

## Reading the result

| $\rho(\zeta, \mathrm{rtslope})$ sign+magnitude | Interpretation |
|---|---|
| Strongly negative ($\rho \lesssim -0.4$) | Single underlying efficiency-growth process; faster CDI growth $\Leftrightarrow$ faster RT decline. Supports treating $\sigma_\zeta$ as a unitary "maturation" axis. |
| Near zero ($|\rho| \lesssim 0.2$) | CDI-side and LWL-side growth are independent maturation clocks. The CDI exponent $\delta_{\mathrm{CDI}}$ inflation relative to $\delta_{\mathrm{LWL}}$ is then *not* attributable to IRT thresholding alone. |
| Strongly positive | Surprising; would indicate growth-rate trade-off or shared confound. Worth investigating before publishing either way. |

The sign is what matters most — it determines whether we have one efficiency axis or two — and the magnitude bounds how much of $\sigma_\zeta$ is shared variance vs. readout-specific.

## Optional follow-ups (only if cheap)

- Same three correlations conditional on residualizing per-child mean age, since the longitudinal Peekbank window (13–20 mo) and the CDI window (22+25 mo) overlap only partially.
- Posterior on $\sigma_{\mathrm{rtslope}}$ standalone, plus the $\zeta_i$ marginal SD restricted to the 62 Peekbank-linked kids — the $\sigma_\zeta \approx 3.48$ headline is on the full N=200 sample and may not hold on this subset.

This is meant to be a small ask on top of an existing fit, not a new modeling project. If the joint-MVN parameterization is already in the Stan file, the extraction is one `summary()` call on the right posterior draws.
