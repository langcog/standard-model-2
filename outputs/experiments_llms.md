# LLM / GPU experiments log

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

Background docs: [`outputs/llm_variability_plan.md`](llm_variability_plan.md),
[`outputs/marlowe_pilot_results.md`](marlowe_pilot_results.md),
[`outputs/feng_evaluation_report.md`](feng_evaluation_report.md). Cluster
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

**Artifacts.** [`outputs/feng_evaluation_report.md`](feng_evaluation_report.md);
`data/feng_2026/gpt2_childes_seed{0,42,123}_sigmoids.txt`;
`outputs/figs/longitudinal/feng_chang_bergen_slope_comparison.png`; pipeline in
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
[`outputs/marlowe_pilot_results.md`](marlowe_pilot_results.md);
`outputs/figs/longitudinal/marlowe_data_variance_pilot.png`;
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

## 🟡 L4. Marlowe — developmental-ladder grid (5 seeds × 11 rungs) (2026-06-09, running)

**Design.** 5 seeds [42, 0, 123, 7, 99] × **11 nested** rungs
[0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24]M words = 55 runs. Each rung trained to
convergence via **early-stopping** (patience 4 on eval_loss, load best-val model),
ceiling 40 epochs. Seed = init + that individual's nested input stream (the total
between-individual factor matched to kids' unpartitioned σ_κ).

**Goal.** (1) Per-word developmental sigmoid κ_w — fittable at 11 rungs, unit-mapped
vs kids; (2) **between-seed σ_κ** of the developmental trajectory — the LM analog
of σ_κ ≈ 3.5. The acceleration *and* its individual variability, from one grid.

**Pipeline change.** Added `--early_stopping_patience` to
[`train_gpt2_childes.py`](../model/scripts/feng_eval/train_gpt2_childes.py)
(additive, default-off; forces per-epoch checkpoint + load-best so the final
surprisal eval = competence at convergence). Branch `claude/ladder-experiment`.
Grid SLURM: `model/scripts/feng_eval/marlowe/train_ladder_grid.slurm`.

**Status (2026-06-09).** Smoke = job **341454** (array idx 0 = seed 42 / 0.5M,
real settings) **queued on preempt**, validating the early-stop edit before the
full `--array=1-54` launch. Samples for all 5 seeds × 11 rungs generated.

**Expected (priors from L1–L3).** Between-seed σ_κ likely small (seed SD ≈ 0.01
at full data in L1; data identity negligible in L2), so 5 seeds should bound it;
extend if unstable.

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
- Full cluster how-to and gotchas: the
  [`gpt2-childes-training`](../.claude/skills/gpt2-childes-training.md) skill.

## Open questions / backlog (⚪)

- **Unit-matched curvature.** Fit per-word developmental sigmoids on the L4 11-rung
  grid; map to κ-equivalent and overlay vs kids — is the LM trajectory
  decelerating where kids accelerate, in matched units?
- **Architecture × budget.** Cross the ladder with a within-GPT-2 arch sweep
  (LR / width / depth) — a later addition, not needed for seed-variance-in-development.
- **Compute-controlled C&B control (from L1).** Does the slope shift with many
  more passes of the small corpus (matched-compute)?
- **Larger CHILDES models (from L1).** Does slope shift with model size at fixed
  input distribution?
- **Evanson, Lakretz & King (2023)** 48-seed GPT-2 (WikiText103): potential
  large-n seed-variance anchor *if* per-step checkpoints can be obtained — would
  let us run our per-word pipeline at n=48. Email TBD (see session notes).
