# Autoresearch acceleration arc — standing plan & status

## The question
Children **accelerate**: as they accumulate experience they get *more* out of each
new bit — **increasing returns** (κ_c ≈ 10). Standard LMs are scaling-law learners:
more data helps *less and less* — **diminishing returns**. This arc asks: **can any
simple training intervention make a small from-scratch LM show increasing returns
instead of diminishing — i.e. behave like a child, not a scaling law?**

## The metric
Per model we log per-CDI-word surprisal across training (the probe) → per-word
acquisition curves. From those, the **learner's returns curve**: how competence-gain
per unit experience changes as it trains.

> **J = (gain rate, late) − (gain rate, early).**  J<0 = diminishing (default);
> J≈0 or >0 = the kid-like increasing-returns direction we hunt for.

J (returns curvature) replaced the per-word slope κ_w, which measures "fast at
learning a word" (wrong target). J is the level-matched analogue of children's κ_c.
Scorer: `returns_score.R`.

## Two model setups (don't conflate)
1. **CHILDES ladder** — the paper's real data: GPT-2 to convergence on CHILDES at 18
   budgets. DONE → models show **diminishing returns** (effective κ falls ~1.3→0.5;
   kids ~10, flat). See `LADDER_RETURNS.md`. This is the level-matched "LMs don't
   accelerate" result and bears on the paper's framing.
2. **autoresearch / ClimbMix** — cheap exploratory sandbox (tiny GPT-2 on web text,
   fast single-GPU runs on the ccn2 A40 box). Where we hunt interventions. Pipeline:
   `adapt_a40.py` (FA3→SDPA), `extract_cdi_contexts_ar.py` (365 single-token CDI
   words), `cdi_probe.py` + `apply_probe_patch.py`, `single_pass_patch.py`.

## Experiment ledger
- **Pilot (DONE):** can we measure per-word κ from 5-min / single-pass ClimbMix runs?
  Yes — 96% of CDI words resolve transitions.
- **Bank 1+2+3 (RUNNING):** repetition-depth sweep. Hypothesis (grokking research,
  `research_menu.md`): single-pass-over-fresh keeps words in a smooth data-surplus
  regime → diminishing; **repeating a fixed corpus** (Critical-Data-Size / grokking
  regime) might flip J toward 0. Variants:
  - flat LR: 1,2,3,4,6,8 passes (`patch_rep_constlr.py`, CDI_PASSES env)
  - decay schedule: 1(baseline),4,8 passes (`patch_rep_decay.py`) — schedule control
  - Question: does J climb toward 0 as repetition deepens?
- **Bank 4 (queued):** complete the decay×repetition grid (decay 2,3,6) so flat vs
  decay is matched at every depth.

## Harness
`bank_dispatch.py <manifest>` fans variants across only-currently-free GPUs (one
run/GPU, polite to other users), per-variant patch + env. Score afterward by pulling
each `runs_bank/<label>/word_surprisal.csv` and running `returns_score.R`.

## Next steps (after the repetition grid, roughly ranked)
1. Score the full flat×decay repetition grid → the J-vs-repetition table.
2. If repetition moves J: push deeper / smaller distinct corpus (true Critical-Data-
   Size); if not: that's an informative null.
3. Seed replication (error bars on J vs depth) — needs a seed patch.
4. Flavor-(ii) **late-held-out frequency-matched lexicon** (the decisive learning-to-
   learn test) — needs dataloader work.
5. Grokfast (slow-gradient amplification, ~20 LOC), curriculum, burstiness/Zipfian —
   from `research_menu.md`.

## Why it matters either way
An intervention that flips J positive = a candidate *mechanism* for acceleration
(Discussion material). A thorough null = constructive evidence LMs structurally can't
— which strengthens the paper.

## Compute notes
ccn2 A40 box (`ssh ccn2-14`), conda env `/data2/mcfrank/ladder/condaenv` (torch
2.4.1 + tiktoken/pyarrow/rustbpe), repo at `/data2/mcfrank/autoresearch`. SSH
ControlMaster drops under load (long-lived local watchers are flaky) → prefer on-box
dispatchers + periodic check-ins. ~24 min per single pass on an A40.
