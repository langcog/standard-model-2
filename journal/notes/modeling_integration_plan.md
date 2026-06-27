# Modeling integration plan — Spanish + sumscores, on a simplified Stan model

Status: **converging, NOT final.** More data work pending (HABLA item-level digitization,
more BabyView kids) and an open anchoring decision. The fits below are interim; they de-risk
the architecture before the final data lands.

## Guiding principle
We're on the last few runs, so every modeling phase should also be a SUBTRACTION: strip
complexity that isn't doing work in the paper. Build the Spanish extension on a clean base.

---

## Phase 0 — Simplify the (English) Stan model FIRST
Each removal is verified by a before/after fit that must reproduce κ_pop / γ_in / β_xi /
the variance partition. Audit from `log_irt_long_proc_dp_joint_lean.stan` + the bundle:

- **Rip λ (discrimination).** `sigma_lambda_prior_sd = 0.001` ⇒ λ_j ≈ 1 already. Remove
  `log_lambda_raw`, `sigma_lambda`, `lambda`, and the multiply in `accumulate()` → pure Rasch
  (`eta = admin_base + item_offset`). Bonus: makes "sumscore is a sufficient statistic for θ"
  EXACT, which the count likelihood (Phase 2) relies on.
- **Hard-code β_c = 1** (pinned). Drop the `beta_c` parameter; decide whether the structural
  `log_p` (frequency) offset stays in the paper or goes too.
- **Test the lexical-class difficulty hierarchy.** Currently C=4 (`mu_c`/`tau_c` per class) —
  still ON, contrary to assumption. Refit with C=1 (single global difficulty distribution);
  if κ_pop/γ_in/β_xi don't move, collapse it.
- **Audit + likely remove:** the std (CELF/SEEDLingS) channel (`mu_std`/`sigma_std`), `beta_age`,
  and any vestigial "dp" machinery — keep only if they appear in the paper.
- **Keep the population δ prior WIDE** (the real past backfire: `N(0,0.5)` shrank acceleration,
  ζ compensated — experiments.md ~entry 33). Validate δ_j vs production rate, not frequency.

Output: a lean, paper-faithful Rasch-accumulator io-proc model.

## Phase 1 — Spanish + sumscore data prep
- **Parse SLENA → item-level** (Spanish WG+WS, `data/raw/WF2013/SLENA_*`) → Spanish item bank,
  the same long schema as the English cohorts. ParticipantId == peekbank `lab_subject_id`.
- **Count-likelihood inputs** for sumscore cohorts (ELENA WS `VOCAB`, HABLA WordsProd/CDIVoc):
  per admin a (child, age, form, k_produced, n_on_form) row + the form's item set (so the
  Poisson-binomial knows which difficulties to sum over).
- **Wordbank Spanish (Mexican)** cross-sectional pulled for the Phase-3 scale check (and the
  optional anchor).

## Phase 2 — Bilingual model + count likelihood
- **Non-overlapping item sets per language**, each with its own hierarchical difficulties
  (same machinery). **Shared population geometry**: μ_r, κ_pop, and the channel couplings
  γ_in (input→κ), β_xi (proc→ξ), the σ's — sharing these IS the cross-linguistic claim.
- **Add the count/sumscore likelihood branch** (Poisson-binomial over a form's item set) so
  ELENA-WS + HABLA enter via their totals. ELENA gains a 2nd timepoint → acceleration.
  HABLA enters ONLY this way (no item responses).
- Fit **anchor-free** first.

## Phase 3 — Spanish scale check → anchoring decision
- The real Spanish risk is **scale identification** (~29 SLENA kids weakly tie the Spanish
  difficulty/θ scale to English through the shared priors).
- Check Spanish ξ/κ land in a plausible range vs English. If the scale drifts, add Wordbank
  anchors **symmetrically** (EN→Wordbank-EN, ES→Wordbank-ES) so scales stay comparable —
  per-item anchors were tested benign historically (SD=5 relaxation → identical κ). Else stay
  anchor-free.

---

## Still open (why these are interim fits)
- **Data:** HABLA item-level (digitize?), more BabyView kids.
- **Anchoring:** decided empirically in Phase 3.
- So lock the *architecture* now; re-fit for the paper once data + anchoring settle.

---

## Phase 0 RESULT — verified + adopted (2026-06-26)
50x50 channel-balanced subsample; simplified model reproduces the prior `_mm` lean:
- **Rasch** (rip lambda, delete s, rip log_p): **21/21** key scalars match (max |z|=0.15).
- **C=4 -> C=1** collapse: **18/18** match; well-identified science (delta, gamma_in,
  eff_input_k, partition, sigmas) identical. beta_xi/beta_k0 swing but are prior-dominated
  on 50 kids (subsample noise; item-class hierarchy is structurally independent of the RT
  regressions). Verified-equivalent but **NOT adopted** — keep C from the bundle; revisit on
  full data.

**ADOPTED:** the simplified Rasch content IS now the lean model
(`model/stan/log_irt_long_proc_dp_joint_lean.stan`); temp `log_irt_long_proc_rasch.stan`
removed. delta prior widened to N(0,10) in the driver. Production full-data ladder D'0-D'2
launched (Sherlock array 31622975). Next: Phase 1 (Spanish + sumscore data prep).
