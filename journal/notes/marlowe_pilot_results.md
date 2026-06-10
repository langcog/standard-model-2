# Data-variance pilot — results

**Run 2026-06-08 on Marlowe (Stanford DGX H100).** Companion to
[`journal/notes/marlowe_pilot_runbook.md`](marlowe_pilot_runbook.md) and
[`journal/notes/llm_variability_plan.md`](llm_variability_plan.md).

## Question

Does the *identity* of the CHILDES training data move the per-word
Chang & Bergen sigmoid slope (κ-equivalent), holding architecture, seed,
tokenizer, eval set, and epoch count fixed? This gates the main study's
σ_data axis: if data identity injects little variance, CHILDES subsamples
are usable; if it injects a lot, we'd need disjoint samples CHILDES can't
supply (→ TinyDialogues).

## Design

Two GPT-2 small (124M) models trained on **two random-disjoint 10M-word
CHILDES chunks** (19.00M BPE each, equal to 0.002%), seed 42 for both,
shared `GPT2_CHILDES` tokenizer, shared CDI-context eval set drawn from the
held-out val split (so no leakage asymmetry). 20 epochs, ~46.4k steps,
~1h12m each on one H100. Per-word 4-PL sigmoids fit to the surprisal
trajectories; slope = `0.434 / ParamScale`. 609/608 words fit cleanly
(0 skipped) under the C&B filters (`0.01 < ParamScale < 10`,
`ParamUpper − ParamLower > 1`).

## Result

![data-variance pilot](figs/longitudinal/marlowe_data_variance_pilot.png)

*(Figure regenerates via `Rscript model/scripts/feng_eval/pilot_data_variance_plot.R`;
PNGs are not version-controlled per repo convention.)*

| Metric | Value | Reading |
|---|---|---|
| Pearson r(slope_A, slope_B) | **0.76** | per-word slope is mostly intrinsic to the word |
| Spearman ρ | **0.88** | which words are learned fast is highly stable across data draws |
| Paired Δ (A−B) median | **+0.013** | the typical word's slope barely shifts |
| Paired Δ mean | +0.036 (≈2 SE) | small marginal offset |
| Paired Δ SD | 0.452 | per-word wobble (data effect + fit noise, inseparable at n=2) |
| Marginal median A / B | 0.643 / 0.570 (gap 0.073) | inflated by the right-skewed tail |
| — for scale — | | |
| Seed-only floor (24.5M, 3 seeds) | ~0.01 | |
| Children σ_κ | ~3.5 (median ~10) | 40–270× the data-identity shift |

(Reproduce: `Rscript model/scripts/feng_eval/pilot_data_variance_plot.R`;
numbers in `fits/feng_eval/pilot_slope_summary.csv`.)

## Interpretation

**Data identity is a negligible source of between-LM slope variance.**
Per-word slopes reproduce strongly across two fully-disjoint data samples
(ρ = 0.88), and the per-LM summary slope moves by only ~0.01–0.07 — a few
percent — versus children's σ_κ ≈ 3.5. This is the *"LM between-instance
variance stays small even with data variation activated"* outcome: it
**strengthens the categorical LM-vs-child difference** rather than
undermining it.

**Design implication — CHILDES is fine; the disjointness worry was moot.**
Because data identity injects almost no variance, overlapping CHILDES
subsamples would barely bias a σ_data estimate. There is **no forced move to
TinyDialogues**.

**A larger lever surfaced: training *amount* > data *identity*.** Both 10M
chunks sit at median ~0.57–0.64, below the 24.5M three-seed runs' ~0.73 — a
~0.1–0.15 shift from *quantity*, larger than the ~0.01–0.07 shift from
*identity*. (Partly the known partial-training bias: fewer steps → lower
plateau → underestimated slope.) The most informative main-study axis may be
**token budget**, not data subsample — and it maps onto Coffey & Snedeker's
"input quantity" share on the human side.

## Caveats (n = 2)

- One pair → the paired Δ is a single realization, not an SD; no CI on σ_data
  from two runs. But ρ = 0.88 is strong evidence the eval is stable and data
  identity is minor.
- 1 − r² ≈ 42% of per-word slope variance is unshared, conflating a true
  data effect with fit noise — inseparable without more draws.
- Absolute slopes are depressed by the 10M scale (plateau bias); only the
  within-scale A-vs-B comparison is valid. Do not compare the 10M chunk
  medians directly to the 24.5M seed medians as if same-scale.

## Artifacts

- `data/feng_2026/gpt2_childes_chunk{A,B}_seed42_sigmoids.txt` — per-word 4-PL fits
- `fits/feng_eval/surprisal_gpt2_childes_chunk{A,B}_seed42.csv` — raw surprisal trajectories (72 log-spaced steps × 609 words)
- `fits/feng_eval/pilot_split_manifest.json` — the random-disjoint split (line indices, token counts)
- `fits/feng_eval/pilot_cdi_coverage.csv` — CDI single-token / occurrence coverage
- `fits/feng_eval/pilot_slope_summary.csv` — the statistics above
- `figs/longitudinal/marlowe_data_variance_pilot.png` — the figure
- `model/scripts/feng_eval/pilot_data_variance_plot.R` — reproducible stats + figure

## Next (deferred)

Main-study axis choice parked pending discussion. Two candidates noted by this
pilot: (a) a **token-budget** sweep (5M/10M/24.5M) to trace the
quantity→slope relationship, which looks like the bigger lever; (b) **more
disjoint draws at smaller scale** (e.g. 3M × 6–8) to put an actual SD on the
~0 σ_data rather than a single A-vs-B difference. The cluster is committed to
the hero run through July 1; Cycle 6 opens then.
