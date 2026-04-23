# σ_r sensitivity sweep — findings

How robust is the RQ4 answer (π_α) to the externally-pinned
input-rate variance?

## Setup

Re-fit the 2PL variant on the same 500 × 200 subsample under four
values of σ_r (the SD of log input rate, pinned from external data):

| σ_r | 90% range of input (tokens/hr) | interpretation |
|---:|---|---|
| 0.30 | 730–2,100 | narrow — within a homogeneous demographic |
| 0.53 | 500–2,850 | Sperry/HR/WF pooled (default) |
| 0.80 | 320–4,450 | roughly 100k–1.5M words/month user range |
| 1.20 | 170–8,550 | very wide — cross-SES and cross-culture |

Each fit used 2 chains × 1,500 iter (750 warm-up). Runtimes 49–66 min.

## Results

| σ_r | σ_α | σ_xi² = σ_r² + σ_α² | **π_α** | CrI |
|---:|---:|---:|---:|---|
| 0.30 | 2.18 | 4.84 | **0.981** | [0.978, 0.984] |
| 0.53 | 2.13 | 4.82 | **0.941** | [0.931, 0.950] |
| 0.80 | 2.05 | 4.84 | **0.868** | [0.849, 0.887] |
| 1.20 | 1.83 | 4.79 | **0.699** | [0.657, 0.740] |

## What this shows

1. **Total child-level variance σ_xi² is pinned at ~4.82 across every
   fit**, even though σ_r varies 4× and π_α swings by 28 percentage
   points. **The data identifies the total child variance precisely; it
   does not identify its decomposition into input vs. efficiency.** The
   decomposition comes entirely from the σ_r prior.

2. **π_α is highly sensitive to σ_r.** Under the "narrow input
   variation within a homogeneous university-research sample" prior, π_α
   is 0.98. Under the "input varies 50× across SES" prior, π_α is 0.70.

3. **Even at the widest-credible σ_r = 1.2, input explains at most
   30% of child-level variance.** The robust finding — bounded below by
   the worst case — is that **individual-difference efficiency still
   dominates input rate in explaining vocabulary variation at this age
   range**. It just might be 70% rather than 94%.

## Implication for the paper

The RQ4 headline should be stated as a robustness interval rather than
a point estimate:

> Depending on how much input rates vary across the population the
> CDI sample represents (σ_r ∈ [0.3, 1.2]), individual differences in
> learning efficiency account for **70–98% of between-child variance in
> vocabulary size** at 16–30 months, and input rate for the complement.
> The total child-level variance on the logit scale is tightly
> identified (σ_xi² ≈ 4.8); its decomposition is prior-bound.

This is stronger than a single number because it honestly reports what
the data alone can and cannot identify. The Sperry/HR/WF σ_r = 0.53
gives our central estimate of π_α = 0.94, but readers who believe the
true population input variance is wider can plug in their own σ_r and
read off the corresponding π_α.

## What we need to pin this down better

The within-child (BabyView/Seedlings) plus between-child (Sperry/HR/WF)
decomposition in the backlog would give us a defensible σ_r. If
within-child variance is large, the *effective* between-child variance
for CDI-age purposes may be smaller than the Sperry pooled number
suggests, pushing π_α toward the upper end of our interval. If
within-child variance is small, σ_r is well-approximated by the pooled
estimate.
