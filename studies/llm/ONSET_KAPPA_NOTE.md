# "Early vs late-learned words" (onset→κ) on the CHILDES LMs — and why it's confounded

**Question.** Per word, does its *onset* (when it is acquired) predict its *κ* (how fast)? A
positive onset→κ, controlling frequency, = "later-learned words are learned faster" = the child
acceleration signature. Run on the cached CHILDES per-word sigmoid fits (no retraining):
`studies/llm/onset_kappa_childes.R` (checkpoint axis = `fits/llm/sigmoids/*.txt`; development axis =
per-word sigmoids fit to `fits/llm/ladder_bestval_finer.csv`); frequency from
`figs/longitudinal/psi_freq_regression_per_word.csv`.

**Raw result — looks like acceleration (both axes, significant):**
| axis | filter | n | partial cor(κ, onset \| freq) | onset β (p) |
|---|---|---|---|---|
| checkpoint (step) | all | 486 | +0.22 | +0.24 (1e-6) |
| checkpoint (step) | asymptote-reached | 408 | **+0.38** | +0.47 (2e-15) |
| development (budget) | all | 1485 | +0.15 | +0.49 (6e-9) |
| development (budget) | asymptote-reached | 1291 | **+0.32** | +1.46 (7e-32) |

This *contradicts* both the paper's thesis (LMs show no acceleration) and the aggregate
diminishing-returns finding (κ falls 1.25→0.47). Red flag → stress-tested it.

**The confound (proven by simulation).** Fitting a sigmoid on a *bounded* experience window
mechanically links onset and scale: a word whose transition starts late has less room before the
window ends, so its fitted transition is compressed → smaller scale → **inflated κ**, with no real
acceleration. Simulating words with **constant true κ (=1.24, no onset→κ by construction)** on the
real 18-budget grid, fit with the same nls, recovers:

- cor(κ, onset): **+0.07 (all)**, **+0.42 (asymptote-reached)**  ← reproduces the CHILDES +0.32

So the observed positive is essentially entirely the bounded-window fitting artifact — and the
"asymptote-reached" truncation filter *amplifies* it (0.07→0.42) rather than removing it. The real
onset→κ does not exceed the constant-κ null → **no evidence of genuine word-level acceleration.**
(This is exactly the truncation caveat flagged in the original `onset_kappa.R`.)

**Conclusion.** The per-word onset→κ ("early vs late-learned words") cut is **not a valid measure of
learner acceleration for LMs** — the sigmoid-on-bounded-axis geometry manufactures a positive slope
from flat-κ data. **Do not put the naive version in the paper**; it would spuriously suggest the LMs
accelerate.

**Robust alternative for the same scientific point.** The aggregate returns / κ-decline analysis
(`ladder_returns.R`, `LADDER_RETURNS.md`) measures improvement-per-decade *directly* (not via
per-word midpoints), so it is immune to this artifact: on the converged 10-seed ladder, κ_eff falls
**1.25 → 0.47** across the budget range (quadratic curvature b₂=+0.77; 84% of words diminishing) —
i.e. the model **decelerates** as it accumulates experience, the mirror image of children. That is
the artifact-free way to make the "later ≠ faster" point, and a natural sharpening of the datasize κ
already in Fig 6 (currently summarized as a flat ~1.18, which hides the decline). Sim:
`/tmp/sim_onset.R`.
