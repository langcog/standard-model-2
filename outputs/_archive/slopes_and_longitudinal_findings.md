# Per-child slopes + longitudinal evidence — findings

Three analyses bearing on the same question: do children genuinely differ
in vocabulary *growth rate*, or only in *level*?

## 1. Cross-sectional 2pl_slopes fit on CDI:WS subset (500×200)

Scalars (compared to 2PL without slopes):

| param | Rasch | 2PL | 2PL + slopes |
|---|---:|---:|---:|
| σ_α | 1.98 | 2.12 | 2.12 |
| π_α | 0.93 | 0.94 | 0.94 |
| s | 12.2 | 12.9 | 12.8 |
| δ | 2.31 | 2.20 | 2.28 |
| σ_λ | — | 0.275 | 0.276 |
| **σ_ζ** | — | — | **0.16 [0.007, 0.57], Rhat 1.15, n_eff 23** |

**Verdict: σ_ζ is not identifiable from cross-sectional data.** The
posterior hugs zero but has a long right tail; Rhat 1.15 and n_eff 23
mean the chains haven't found a stable mode. The heteroskedasticity-
based identification argument fails in practice with only a 14-mo age
range and 500 kids.

Adding slopes did *not* reduce σ_α. That variance is robust to every
extension we've tried.

## 2. Within-age ability SD PPC

New PPC panel (Plot 8): observed SD of logit(vocab-proportion) within
1-month age bins, vs. the model-implied $\sqrt{\sigma_\xi^2 + \sigma_\zeta^2 L(a)^2}$.

- **Observed**: rises from ~1.0 at 16 mo to ~1.6 at 30 mo (consistent
  with growing heterogeneity).
- **Predicted by the 2PL+slopes model**: nearly flat at ~2.2, well above
  observed.

So the model over-predicts within-age variation by ~30–50%. This is the
same σ_α inflation visible other ways, but now localized to a
testable prediction: the current model has σ_ξ too large to reproduce
observed within-age dispersion. It fits the *marginal* vocabulary
distribution because the Bernoulli likelihood averages over all that
excess latent variance, but the within-age SD of logit-vocab is a
sharper probe.

## 3. Longitudinal LMM on admin-level production counts

Wordbank `admins.feather` has 3,786 English WS admins from 1,653
children measured ≥2 times (median span 7 mo). Fit on logit of
production proportion (i.e., aggregate ability per admin, no item-level
data needed):

```
logit_prop ~ log_age + (log_age | child_id)
```

with `log_age` centered at 24 mo. Results:

- Fixed: intercept −0.26 (prop ≈ 0.44 at 24 mo), slope 7.14 per log-age unit.
- **SD(random intercept) = 1.49** (in logit units — compare to model σ_xi ≈ 2.2)
- **SD(random slope) = 2.43** per unit log-age
- **Correlation (intercept, slope) = +0.72**
- Residual SD = 0.92
- LRT for adding random slope: p < 5e-61

**Verdict: children genuinely vary in vocabulary growth rate**. The
random-slope variance is highly significant and substantial (≈34% of
the fixed slope). Unlike the cross-sectional 2pl_slopes result, this
signal is clear.

**And corr(intercept, slope) = +0.72**: kids who are ahead at 24 mo are
*also* growing faster. This is opposite to our model's assumption of
independent (ξ, ζ). It's consistent with a shared sigmoidal curve where
higher-ability kids at 24 mo are still in the steep middle of the
curve, but it's a strong enough correlation to suggest real coupling.

## 4. How this reconciles with Peekbank (your in-press paper)

Peekbank found faster processors have bigger vocabularies AND faster
vocabulary growth (longitudinal SEM, β = −0.13 for t₀ RT on vocab
slope). But growth-*in-speed* was uncoupled from growth-*in-vocabulary*.

Our longitudinal LMM shows real between-child variance in vocabulary
*growth rate* (σ(slope) = 2.43). This is what Peekbank's t₀-RT-on-growth
finding implies.

The intercept-slope correlation of +0.72 is consistent with the
Peekbank picture: stable processing differences manifest as both higher
current vocab AND faster current growth (because the child is further
along a shared learning curve).

## Implications for the standard model

1. **σ_α ≈ 2 is real and substantial.** None of (2PL, slopes, both)
   has moved it. The longitudinal LMM also shows SD(random intercept) =
   1.49 in the same logit units — a bit smaller than our 2.2 but same
   order of magnitude. Large genuine individual differences in ability
   at any given age.

2. **Genuine growth-rate differences exist**, but the current Stan
   model can't identify them from cross-section. This is the main
   motivation for the upcoming longitudinal-data fit.

3. **ξ and ζ are strongly correlated in the real data (+0.72)**. Our
   current model has them independent, which is a misspecification we
   should fix when we re-run on longitudinal data.

4. **Current π_α = 0.93 should be interpreted carefully.** Most of the
   child-level variance currently attributed to learning efficiency
   likely is real individual-differences in "processing ability" a la
   Peekbank. But some chunk of it is probably "the growth-rate variance
   the cross-sectional model can't separate cleanly" — and that chunk
   is arguably its own construct.

## Next moves

- **σ_r sensitivity sweep (running now)** on the 2PL variant, to check
  how robust π_α is to the externally-pinned input-variance estimate.
  Results will land in `fits/sensitivity_sigma_r_2pl.rds`.
- **Longitudinal accumulator fit** — need item-level longitudinal data,
  which we don't have locally. Admin-level totals + binomial likelihood
  is a middle ground we could pursue if item-level stays out of reach.
