# LLM / GPU experiments log

> **COPY — shared history.** Entries below up to and including the split of the
> `standard_model_2` repository (experiments 1–41, L1–L7) are a frozen shared record,
> copied verbatim into both `acceleration` and `standard_model_2`. They are the history of
> one research program that produced two papers, so neither repo owns them alone.
>
> **New entries go in whichever repo the work happens in**, and are not back-copied. If you
> need the other repo's later entries, read them there.
>
> Orientation for a cold start: [`STANDARD_MODEL_CONTEXT.md`](STANDARD_MODEL_CONTEXT.md).

A running record of the **language-model word-acquisition** experiments —
the LM side of the children-vs-LMs comparison. Kept separate from
[`experiments.md`](experiments.md) (the Stan / CDI psychometric-model log)
because this arc runs on GPUs (Sherlock, Marlowe), uses a different pipeline
([`model/scripts/feng_eval/`](../model/scripts/feng_eval/)), and answers a
distinct question.

**The through-line.** Children's per-child sigmoid slope on log-experience is
κ ≈ 10 (English M_best); LMs sit near 1 on the same C&B statistic. We are
(1) checking whether that gap is an input-distribution artifact, (2) putting an
*apples-to-apples* between-instance variance number on the LM side (kids'
σ_κ ≈ 3.5), and (3) re-asking the acceleration question on a **development-matched
axis** (distinct input, not training steps).

Background docs: [`journal/notes/llm_variability_plan.md`](notes/llm_variability_plan.md),
[`journal/notes/marlowe_pilot_results.md`](notes/marlowe_pilot_results.md),
[`journal/notes/feng_evaluation_report.md`](notes/feng_evaluation_report.md). Cluster
how-to: the [`gpt2-childes-training`](../.claude/skills/gpt2-childes-training.md)
skill.

---

## Status key
- 🟢 completed
- 🟡 running / active
- ⚪ queued / backlog

---

## 🟢 L1. Sherlock — C&B per-word sigmoids on CHILDES-trained GPT-2, 3 seeds (2026-05)

**Question.** Is the C&B kid-vs-LM per-word-slope gap (≈10 vs ≈1) partly an
input-distribution artifact? C&B trained on BookCorpus/WikiText (adult written
text); children hear child-directed speech. Retrain GPT-2 on CHILDES and refit
the same per-word 4-PL sigmoid, so only the structural axis remains.

**Path.** Feng et al. 2026 (`styfeng/babyscale-LM`) released no checkpoints or
per-step logs, so we retrained from scratch (**Path B**) using the
`styfeng/TinyDialogues` pipeline (the code Feng used for the CHILDES condition),
logging per-CDI-word surprisal at log-spaced steps.

**Setup.** GPT-2 small (124M, 12L×768d×12h), `GPT2_CHILDES` BPE (52K vocab),
`CHILDES_train_ordered.txt` (~24.5M tokens), 20 epochs, AdamW LR 1e-4 linear /
no warmup, batch 8, seq 1024, no in-epoch shuffle (mirrors Feng §B). Seeds
{42, 0, 123}. Surprisal at 73 log-spaced steps, 50 occurrences/word, 128-token
left context. **114,520 steps** total (45,807 1024-token blocks × 20).

**Compute.** 1× **L40S** 48GB per seed, Sherlock `gpu` partition; **7h21m /
8h12m / 7h57m** wall. (Seed 0 first landed on a V100, ~3× slower, cancelled at
11% and resubmitted with `--constraint=GPU_GEN:AMP|LOV|HPR`.)

**Coverage.** 611/611 C&B CDI words are single tokens in GPT2_CHILDES; 609 have
≥1 val occurrence, 578 have ≥50.

**Result.** Per-word slope = 0.434/ParamScale:

| population | N | median | IQR |
|---|---|---|---|
| Children κ_i (English M_best) | 5000 | **10.3** | [8.0, 12.6] |
| GPT-2-CHILDES seed 42 | 609 | **0.74** | [0.43, 1.11] |
| GPT-2-CHILDES seed 123 | 609 | **0.74** | [0.45, 1.16] |
| GPT-2-CHILDES seed 0 | 609 | **0.72** | [0.45, 1.14] |
| GPT-2-BookCorpus (C&B) | 604 | 0.81 | [0.45, 1.54] |
| BERT / BiLSTM / LSTM (C&B) | ~600 | 0.76 / 0.87 / 0.96 | — |

**Finding.** CDS-matched training does **not** move the per-word slope (0.72–0.74,
inside the BookCorpus cluster 0.76–0.96); seed-to-seed SD ≈ 0.01. **Input
distribution accounts for ~0 of the 10× gap — it is structural.**

**Partial-fit bias (important methodological note).** Mid-training fits
*overestimate* ParamScale (underestimate slope); the median evolved
3.22 → 2.36 → 1.85 → 1.35 → 0.84 → 0.69 → **0.74** as training reached
convergence. **Only fully-trained fits are the meaningful comparand.**

**Artifacts.** [`journal/notes/feng_evaluation_report.md`](notes/feng_evaluation_report.md);
`data/feng_2026/gpt2_childes_seed{0,42,123}_sigmoids.txt`;
`figs/longitudinal/feng_chang_bergen_slope_comparison.png`; pipeline in
[`model/scripts/feng_eval/`](../model/scripts/feng_eval/); SLURM
`sherlock/feng_train_gpt2.slurm`.

---

## 🟢 L2. Marlowe — data-variance pilot: 2 disjoint 10M-word CHILDES chunks (2026-06-08)

**Question.** Does the *identity* of the CHILDES training data move the per-word
slope, holding architecture/seed/tokenizer/eval/epochs fixed? Gates the main
study's σ_data axis (and the CHILDES-vs-TinyDialogues choice).

**Setup.** 2× GPT-2 small on **two random-disjoint 10M-word** CHILDES samples
(19.00M BPE each, equal to 0.002%), seed 42 both, shared tokenizer + shared
held-out CDI eval set, 20 epochs. Split via
[`make_disjoint_chunks.py`](../model/scripts/feng_eval/make_disjoint_chunks.py).

**Compute.** ~**1h12m / run** on 1× **H100** (Marlowe `preempt`). Coverage 609/611.

**Result.** Per-word slopes across the two disjoint samples: **Pearson 0.76,
Spearman 0.88**. Paired Δ(A−B): median **+0.013**, mean +0.036; marginal medians
0.643 / 0.570 (gap 0.073). The per-LM summary slope shifts only **~0.01–0.07** vs
children's σ_κ ≈ 3.5 — **40–270× smaller**.

**Finding.** **Data identity is a negligible variance source.** CHILDES is
usable for the σ_data axis (overlap would barely bias it) → **no forced move to
TinyDialogues**. Training *amount* (10M ~0.6 vs 24.5M ~0.73) looks like a bigger
lever than data *identity* — which motivated the developmental ladder (L3).

**Caveats.** n=2 (a single A-vs-B difference, no CI); 1−r² ≈ 42% unshared
per-word variance conflates true data-effect with fit noise; 10M absolute slopes
are depressed by plateau bias, so only the within-scale A-vs-B comparison is valid.

**Artifacts.** `data/feng_2026/gpt2_childes_chunk{A,B}_seed42_sigmoids.txt`;
[`journal/notes/marlowe_pilot_results.md`](notes/marlowe_pilot_results.md);
`figs/longitudinal/marlowe_data_variance_pilot.png`;
`model/scripts/feng_eval/pilot_data_variance_plot.R`. **PR #20 (merged).**

---

## 🟢 L3. Marlowe — developmental ladder, sweep 1 (1 seed × 6 budgets) (2026-06-08/09)

**Motivation.** The C&B slope is over *training steps* — re-seen passes of a
fixed corpus (optimization convergence), **not development** (accumulating
*distinct* input). Train separate models to convergence at increasing
distinct-input budgets; competence vs log(budget) is the development-matched
acceleration. Each seed is an "individual" with its own nested input stream.

**Setup.** 1 seed (42) × **6 nested** CHILDES budgets [1, 2, 4, 8, 16, 24]M words
(cumulative per seed; 24M ≈ full corpus), fixed tokenizer + eval. Per-budget
epoch caps [20,20,20,20,15,10] to stay under preempt's 4h cap. Competence read at
the **best-val** point (min per-epoch eval_loss). Sampler:
[`make_ladder_samples.py`](../model/scripts/feng_eval/make_ladder_samples.py).

**Compute / wall-time** (1× H100, preempt): 1M **24m**, 2M 29m, 4M 40m, 8M 60m,
16M 79m, 24M 78m. A **~20-min fixed overhead floor** (72 surprisal evals +
per-epoch val + 461MB context load) dominates small budgets — not a budget-scaling
artifact.

**Convergence (best-val epoch).** Falls with budget: 1M ~20 (still creeping),
2M 12, 4M 10, 8M 8, 16M 7 — but **24M was NOT converged at 10 epochs** (val still
dropping). → fixed epochs are a guessing game; use early-stopping (see L4).

**Result — the developmental trajectory** (best-val mean held-out CDI surprisal):

| budget | 1M | 2M | 4M | 8M | 16M | 24M |
|---|---|---|---|---|---|---|
| surprisal (nats) | 6.89 | 6.03 | 5.43 | 4.94 | 4.60 | 4.34 |
| Δ per e-fold | — | −1.25 | −0.87 | −0.71 | −0.51 | −0.63 |

**Finding.** Competence improves with **diminishing returns — concave in
log-input** (the neural scaling-law shape, surprisal ≈ budget^−β). That is the
**opposite curvature from children**, who *accelerate* (super-linear vocab
growth, κ ≈ 10). On a properly development-matched axis, LMs decelerate.

Rough κ-mapping (0.434/ParamScale on surprisal vs log10-budget): aggregate
**≈ 0.63** — in line with the C&B *training*-axis 0.74 (so the shallow slope is
not an axis artifact), ~15× below kids. Per-word κ_w at 6 rungs is unstable
(median 1.23, IQR [0.75, 2.25], 548/609 fit) → **need more rungs**.

**Caveats.** Raw surprisal (nats) ≠ kids' logit-ability scale — the curvature-vs-kids
claim wants the unit-matched sigmoid readout, which needs ≥10 rungs. 6 points
can't fit a per-word 4-PL.

**Decisions → L4.** Finer ladder (~11 rungs), early-stopping (auto-converge each
rung, fixes 24M), several seeds for σ_κ. Within-training trajectory is
unnecessary for the ladder (only best-val competence matters) → trim eval to cut
the overhead floor.

**Artifacts.** `model/scripts/feng_eval/marlowe/train_ladder.slurm`; surprisal
CSVs in `runs_ladder/` on Marlowe scratch (`/scratch/m000102/mcfrank/llm_var_pilot`).
(Not yet PR'd — bundling with L4.)

---

## 🟢 L4. ccn2 A40s — developmental-ladder grid (5 seeds × 11 rungs) (2026-06-10)

**Design.** 5 seeds [42, 0, 123, 7, 99] × **11 nested** rungs
[0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24]M words = 55 runs. Each rung trained to
convergence via **early-stopping** (patience 4 on eval_loss, load best-val model),
ceiling 40 epochs. Seed = init + that individual's nested input stream (the total
between-individual factor matched to kids' unpartitioned σ_κ).

**Where it ran.** Moved off Marlowe (its `preempt` queue stalled under a hero run)
to the lab **ccn2 A40 box** (8× A40, no SLURM). Round-robin dispatcher over free
GPUs: [`ccn2/run_ladder_grid_ccn2.sh`](../model/scripts/feng_eval/ccn2/run_ladder_grid_ccn2.sh).
~38 A40-GPU-hours overnight; **55/55 clean** (exit 0).

**Pipeline changes** to [`train_gpt2_childes.py`](../model/scripts/feng_eval/train_gpt2_childes.py)
(both additive, default-off): `--early_stopping_patience` (per-epoch checkpoint +
load best-val, so the on_train_end surprisal eval = competence at convergence —
verified `load_best_model_at_end` loads the best checkpoint) and `--max_eval_blocks`
(subsample the per-epoch val eval used for early-stopping; 99s→6s/epoch, 16×).

**Result** (reproduce: [`ladder_analysis_final.R`](../model/scripts/feng_eval/ladder_analysis_final.R);
data [`fits/feng_eval/ladder_bestval.csv`](feng_eval/ladder_bestval.csv)):

- **Trajectory decelerates** across all 5 seeds (concave in log-input; mean
  held-out CDI surprisal 7.8 → 4.6 nats over 0.5M → 24M words) — opposite to
  children's acceleration, now replicated across individuals.
- **Per-word dev-axis κ ≈ 1.17** (per-seed medians 1.05–1.24; 0.434/ParamScale,
  unit-matched to the training-axis 0.74 and kids' ~10). The shallow LM slope is
  **not** an axis artifact — it survives the move from training-step to
  distinct-input scaling.
- **Between-seed σ_κ ≈ 0.08 (CV 7%)** vs children's σ_κ ≈ 3.5 (CV ~35%). LM
  "individuals" develop near-identically — between-instance variability is
  categorically absent even with seed + data varied jointly.
- **Individuals converge with input**: between-seed competence SD shrinks from
  ~0.12 nats (0.5M) to ~0.02–0.03 (12–16M) — the *opposite* of children's
  age-divergence. (Overlap-confound check → L5.)

**Headline.** The LM-vs-kid disanalogy now stands on three matched-axis legs:
**rate** (κ ~1.2 vs ~10), **variability** (CV 7% vs 35%), **shape** (decelerating
vs accelerating).

**Artifacts.** `ladder_analysis_final.R`, `fits/feng_eval/ladder_bestval.csv`,
`ladder_kappa_summary.csv`, `figs/longitudinal/ladder_development_final.png`
(figure regenerates from the script).

---

## 🟢 L5. ccn2 A40s — disjoint-CHILDES-halves ladder control (2026-06-10)

**Question.** Does L4's "individuals converge with input" survive *genuinely
disjoint* data? In L4 two seeds shared ≈ B/24.5M of their nested data (~25% at 6M,
~98% at 24M), so top-rung convergence was partly mechanical. This removes the
confound.

**Design.** CHILDES split into **2 disjoint random halves** (poolA/poolB,
~12.2M words each; A∩B = ∅ at every budget) via `make_disjoint_chunks.py
--target_tokens 0`. 2 pools × 3 seeds × 8 rungs (0.5M–12M) = **48 runs**, same
pipeline. Pool = fixed disjoint split; seed = within-pool init+shuffle. Dispatchers:
[`ccn2/run_disjoint_ladder_ccn2.sh`](../model/scripts/feng_eval/ccn2/run_disjoint_ladder_ccn2.sh);
analysis [`disjoint_analysis.R`](../model/scripts/feng_eval/disjoint_analysis.R)
(data [`fits/feng_eval/disjoint_bestval.csv`](../fits/feng_eval/disjoint_bestval.csv)).

**Result (48/48 clean). The overlap confound is ruled out — convergence is real.**
Per-individual developmental slope: **between-pool gap = 0.021** (Pool A −0.962 vs
Pool B −0.983, *zero shared data*) ≈ **within-pool between-seed SD = 0.026** ≈ the
**L4 nested-grid between-seed SD = 0.021**. So genuinely disjoint training data
produces *no more* between-individual variability than same/overlapping data. The
between-pool competence gap shrinks with input just like the within-pool seed
spread (0.08→0.001 nats over 0.5M→12M), tracking each other throughout. Both
pools' six trajectories intermingle (`figs/longitudinal/disjoint_control.png`).

**Takeaway.** L4's "individuals converge with input" is not a B/24.5M overlap
artifact — models trained on disjoint halves of CHILDES develop at the same rate
and converge to the same competence. Caveat retained: only 2 pools, and CHILDES
halves are distributionally similar (low power to *find* divergence) → the
affirmative leg is the BabyLM register-mix control (backlog).

---

## 🟢 L6. ccn2 A40s — finer ladder + more seeds (10 seeds × 18 budgets) (2026-06-11)

**Goal.** Firm up two CIs: per-word κ_w (between-WORD — L4's 11-point sigmoids were
fit-noisy) and σ_κ (between-SEED — 5 seeds is a thin SD). Cheap given the proven pipeline.

**Design.** Extended L4's grid to **18 distinct-input budgets** (0.5M…24M; 7 rungs
inserted — 0.75/1.25/1.75/2.5/3.5/5/7M) and **10 seeds** (added 1–5 alongside the
original 0/7/42/99/123). Per seed: nested best-val competence at each budget, same
early-stopping pipeline. New work beyond L4's 55 runs: 5 new seeds × 18 + 5 original
× 7 inserts = **125 runs**; 180 (seed,budget) cells total. Dispatcher
[`cluster/ccn2/run_finer_ladder_ccn2.sh`](../cluster/ccn2/run_finer_ladder_ccn2.sh);
best-val extraction — max-step per word = the converged / `load_best_model_at_end`
checkpoint, **verified to reproduce L4's `ladder_bestval.csv` byte-for-byte** —
[`studies/llm/extract_ladder_bestval.py`](../studies/llm/extract_ladder_bestval.py).

**Result — L4 confirmed, CIs tightened.**
- **Per-word dev-axis κ**: median **1.19** (10 seeds, n≈5.8k word-fits; L4 5-seed: 1.17).
  This is the development-axis density (green) in the paper's `fig-llm-acceleration`.
- **Between-seed**: κ_pop **1.19**, **σ_κ 0.10**, **CV 8.4%** (L4 5-seed: 1.17 / 0.082 / 7%).
  σ_κ nudged *up* with more seeds — 5 underestimated it — but still ~4× tighter than
  children (σ_κ ~3.5, CV ~35%). Per-seed κ_med range 1.05–1.38, all near 1.
- No qualitative change: rate ≈1.2, variability CV ~8%, decelerating/converging — the
  three-leg disanalogy holds, now with firmer numbers.

**Paper port.** Canonical data behind `fig-llm-acceleration` (PR #47):
`paper/build_cache.R` §6 refits per-word κ from `fits/llm/ladder_bestval_finer.csv`
into `fig6_llm_slopes.rds`. Re-render after any future fits = overwrite that CSV +
re-run `build_cache.R`.

**Artifacts.** `fits/llm/ladder_bestval_finer.csv` (10×18, 109.6k rows),
`studies/llm/ladder_analysis_final.R` (now reads the finer file),
`fits/llm/ladder_kappa_summary.csv`, `figs/llm/ladder_development_final.png`
(regenerates from the script).

---

## 🟡 L6-archive. ccn2 — re-train + share the L6 ladder for reuse (2026-07-17, running)

**Why.** L4/L6's 180 checkpoints were discarded after surprisal extraction
(`--save_total_limit 1` then `rm -rf checkpoint-*`; `/data2` was 100% full). They're
needed for multiple downstream projects + HF sharing → re-train the **exact** L6 ladder
(byte-identical `chunks/`, same config) and stream each model to HF.

**What.** 10 seeds × 18 rungs = 180 GPT-2-small, re-trained from the intact chunks with the
L6 config (40 epochs, ES patience 4, batch 8, LR 1e-4). Only change vs L6: instead of
deleting the best checkpoint, export **bf16 safetensors → HF** (`mcxfrank/childes-gpt2-ladder`,
private) + save surprisal CSVs. Atomic-claim per-GPU dispatcher (`worker.sh`), 8-wide.

**Fidelity** (smoke, seed42/0.5M): re-trained per-word surprisal vs published
`ladder_bestval_finer.csv` — **r = 0.988**, aggregate mean within **0.04 nats**. Faithful
re-training (GPU nondeterminism → not bit-identical; published CSVs remain ground truth).

**Status.** ~140/180 uploaded, 0 failures. ~207 GPU-hr / ~1 day 8-wide.
**Artifacts.** `mcxfrank/childes-gpt2-ladder` (HF, private); `retrain_dev/surprisal_*.csv`;
`worker.sh` + `export_and_upload.py` on ccn2.

---

## 🟡 L7. ccn2 — composition control: BabyLM vs ClimbMix (2026-07-17, queued)

**Question — the affirmative leg for L5.** L5 split CHILDES halves (distributionally near-
identical, low power to *find* divergence). Do LM "individuals" still converge in developmental
rate on *maximally different* corpora? Between-corpus (child-oriented vs web) is the strong lever;
disjoint subsets within each corpus give the sample-variance floor.

**Design.** 6 datasets = 3 disjoint **BabyLM-2024** subsets (CHILDES source **excluded** to avoid
overlap with our ladders — ~23.4M words each from the 70.2M non-CHILDES: gutenberg/open_subtitles/
simple_wiki/bnc/switchboard) + 3 disjoint **ClimbMix** subsets (24M each,
`karpathy/climbmix-400b-shuffle`). **8 seeds × 12 rungs [0.5…24M] = 576 runs.** Identical
measurement to L6 (`GPT2_CHILDES` tokenizer, CHILDES-val CDI probe, ES pipeline). Variance model:
`κ ~ 1 + corpus + (subset|corpus) + (seed|subset)` → σ_corpus / σ_subset / σ_seed vs children σ_κ≈3.5.

**Zero-GPU go/no-go — CDI-word frequency across corpora** (surprisal ∝ log-freq, so freq gates
learnability *and* previews the contrast):

| corpus | words | CDI present | ≥1/M | ≥10/M | median/M |
|---|---|---|---|---|---|
| CHILDES (ref) | 24.5M | 611/611 | 609 | 588 | 104 |
| BabyLM non-CHILDES | 70.1M | 611/611 | 602 | 493 | 62 |
| ClimbMix (web) | 124M | 611/611 | 599 | 489 | 49 |

**GREEN**: all 611 CDI words present in every corpus (probe not floor-bound on web). The
composition contrast is **modest** — CDI log-freq profiles correlate r = 0.76 (CHILDES↔ClimbMix),
0.84 (↔BabyLM), 0.89 (BabyLM↔ClimbMix). Since LM surprisal ≈ log-freq, this **predicts convergence**
even child-vs-web, and is a mechanistic account of the L4/L5 null (cross-corpus convergence =
corpus-invariant frequency structure). MCF committed to genuine interpretation either way.

**Compute/infra.** ~5 days 8-wide. Auto-queued after L6-archive by a self-healing `supervisor.sh`
(smoke-gate → 8 workers → dead-worker restart + stale-claim reap; resumable). Models →
`mcxfrank/gpt2-composition-control` (HF, private, bf16). Builder `make_register_data.py`; spec
[`notes/babylm_register_control_spec.md`](notes/babylm_register_control_spec.md); CDI-freq table
`cdi_freq_corpora.csv`.

---

## Infrastructure & environment

- **Data provenance.** All self-contained in the public `styfeng/TinyDialogues`
  repo: CHILDES corpus (`data/CHILDES_data.zip`, Git-LFS), `GPT2_CHILDES`
  tokenizer, GPT-2-small config. No Sherlock rsync / private code needed.
- **Marlowe.** `preempt` partition (4h cap, evictable), `#SBATCH -G 1`,
  `-A marlowe-m000102`, group scratch `/scratch/m000102`, `conda` module +
  torch cu124, **no git-lfs** (fetch corpus via `media.githubusercontent.com`),
  SSH ControlMaster for the password+Duo login. H100s ~2× L40S/step.
- **Sherlock.** L40S/Ampere via `--constraint=GPU_GEN:...`, venv-on-python-module,
  Kerberos/GSSAPI. ~7–8h per full 24.5M-token run.
- **ccn2 (lab A40 box).** 8× A40 48GB, **no SLURM** — pick GPUs via
  `CUDA_VISIBLE_DEVICES`, round-robin dispatcher keeps 1 run/GPU (GPT-2-small ~31GB,
  so 1/GPU). Shared multi-tenant: reclaim parked GPUs with `sudo kill` (lab owns it).
  Self-contained env on `/data2` (miniconda via conda-forge channel + `ensurepip`;
  `/` is near-full). SSH ControlMaster reseeded by an interactive login; the socket
  drops under heavy node load, so monitoring tolerates hiccups while jobs run
  detached. ~38 GPU-hours for the full 55-run ladder, overnight, free.
- Full cluster how-to and gotchas: the
  [`gpt2-childes-training`](../.claude/skills/gpt2-childes-training.md) skill.

## Open questions / backlog (⚪)

- ~~**Finer ladder.**~~ ✅ Done — **L6** (18 budgets × 10 seeds): κ_w median 1.19,
  σ_κ 0.10 (CV 8.4%). Confirmed L4, tightened CIs.
- **BabyLM-100M disjoint control (Option 3, the affirmative leg for L5).** 4 disjoint
  ~24M subsets of the BabyLM mix (heterogeneous register: web + spoken), 2 seeds each.
  Genuinely *compositionally* different disjoint data — if individuals still converge,
  the L4/L5 convergence claim is strong. Hold the `GPT2_CHILDES` tokenizer + CHILDES-val
  CDI probe fixed; split at the document level (unstratified) so subsets differ in mix.
- **Architecture × budget.** Cross the ladder with a within-GPT-2 arch sweep
  (LR / width / depth) — a later addition, not needed for seed-variance-in-development.
- **Compute-controlled C&B control (from L1).** Does the slope shift with many
  more passes of the small corpus (matched-compute)?
- **Larger CHILDES models (from L1).** Does slope shift with model size at fixed
  input distribution?
- **Evanson, Lakretz & King (2023)** 48-seed GPT-2 (WikiText103): potential
  large-n seed-variance anchor *if* per-step checkpoints can be obtained — would
  let us run our per-word pipeline at n=48. Email TBD (see session notes).
