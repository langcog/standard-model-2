# Model-over-development returns: the level-matched κ comparison

**Question.** Children's acceleration is κ_c at the *child* level — one learner's
competence over *its own development* (κ_c ≈ 10 = increasing returns to experience).
The paper compared that to LMs' *per-word* slope κ_w (word level, within one fixed
model), which is silent about whether the *learner* changes with experience. The
level-matched LM quantity is **the model's competence over DATA** — one learner over
its own accumulating experience. The CHILDES ladder (converged models at increasing
distinct-input budgets) measures exactly this.

**Result (reproduces from `fits/llm/ladder_bestval_finer.csv`, 10 seeds × 18 budgets):**
models show **diminishing returns to data — they decelerate as they develop.**

- Improvement per 10× more data shrinks monotonically: ~2.9 nats/decade early
  (0.5–1M) → ~1.1 nats/decade late (16–24M).
- Aggregate curvature: quadratic coef on log₁₀(budget) = **+0.77** (convex =
  decelerating). **84%** of individual words show the same diminishing sign
  (median per-word b₂ = 0.70).
- For these low-probability CDI words logit-competence ≈ −surprisal, so the slope on
  log-experience is an effective κ. It **falls from ≈1.25 (low data) to ≈0.47 (high
  data).** Children's κ_c ≈ **10, and flat.**

**Interpretation.** This is the clean, level-matched contrast: *children show
increasing returns over development (κ_c≈10); models show diminishing returns over
data (κ falls ≈1.3→0.5).* Same level, converged models (no truncation confound). It
also corrects our earlier per-word development figure of κ≈1.19 — that single sigmoid
was forced through the whole budget range and silently averaged a κ that actually
*declines*; the curvature was hidden by the constant-slope fit.

**Caveat.** Measured on the surprisal scale; surprisal has a floor, but the floor *is*
diminishing returns (nearing what the data can buy), not an artifact. Direction holds
on the logit scale (logit ≈ −surprisal here).

**Bears on the paper.** Suggests reframing the LM-vs-child figure as *model-over-data
(diminishing) vs child-over-development (increasing)* — a stronger, level-matched
claim than the per-word κ comparison.

## Reproduce
```sh
# from repo root; needs fits/llm/ladder_bestval_finer.csv (built by
# studies/llm/extract_ladder_bestval.py from the ccn2 CHILDES-ladder runs)
Rscript studies/llm/autoresearch_pilot/ladder_returns.R
# -> console stats + studies/llm/autoresearch_pilot/ladder_returns.png
```
