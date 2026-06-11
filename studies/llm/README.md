# llm — large language model comparison  ·  Fig 5

Comparison of GPT-2 word learning to children's, on a development-matched axis.

**Start here:** [`llm_experiments.qmd`](llm_experiments.qmd) — the analysis walkthrough
(every claim, which use the Chang & Bergen sigmoid, mapping to paper claims).
Running log: [`/journal/experiments_llm.md`](../../journal/experiments_llm.md).

**Layout**
- **Pipeline (training + scoring):** `train_gpt2_childes.py`, `surprisal_callback.py`,
  `extract_cdi_contexts.py`, `make_ladder_samples.py`, `make_disjoint_chunks.py`,
  `fit_per_word_sigmoid.py`, `cdi_words.txt`.
- **Analysis:** `ladder_analysis_final.R` (developmental ladder, L4),
  `disjoint_analysis.R` (overlap control, L5), `pilot_data_variance_plot.R` (L2),
  `chang_bergen_comparison.R` / `feng_chang_bergen_comparison.R` (training-axis slopes, L1).
- **Cluster scripts:** `cluster/marlowe/`, `cluster/ccn2/`, `cluster/sherlock/feng_*`.
- **Fits / data:** `fits/llm/` — `ladder_bestval.csv`, `disjoint_bestval.csv`,
  per-word `sigmoids/`, raw `surprisal/` trajectories, run `provenance/`.

**Headline:** on a matched axis, LMs differ from children in *rate* (κ ≈ 1.2 vs ~10),
*variability* (between-instance σ_κ ≈ 0.08 vs ~3.5), and *shape* (decelerating vs
accelerating) — and the convergence is not a data-overlap artifact (L5).
