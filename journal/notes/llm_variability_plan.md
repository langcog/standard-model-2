# LLM variability experiment: handoff document

**Authored 2026-05-19** during the SM2 GCP-fitting session. Self-contained
handoff for a new parallel session that will scope and run the experiment.

## What this is

A planned follow-up to the SM2 paper's LM/kid comparison. The current
paper compares the **acceleration** parameter (per-instance scaling
slope `κ`) between children and LMs and finds children ~10× steeper.
This experiment extends the comparison to the **variability**
parameter (between-instance `σ_κ`) — which the current paper claims is
~0 for LMs vs ~3.5 for kids, but on a structurally asymmetric basis
(within-LM-between-word vs between-child).

The proposal: train a population of GPT-2-class LMs that varies along
multiple dimensions, fit Chang & Bergen per-word sigmoids to each,
measure between-LM `σ_κ`. Compare to kids' `σ_κ ≈ 3.5`.

## Why it's worth doing

The current claim has a real structural weakness. For kids:
`σ_κ_kids = SD across children of per-child scaling exponent`. For LMs,
the current comparable quantity is `within-LM, between-word slope spread`
— a different construct. Fixing this makes the apples-to-apples LM-vs-kid
disanalogy claim cleaner.

Two possible outcomes, both publication-worthy:

1. **LM σ_κ stays small even with all sources of variation activated.**
   Strong categorical-difference claim with a clean source-decomposition
   story (data variation contributes X, seed Y, architecture Z, none of
   them enough to reach kids' 3.5).
2. **LM σ_κ reaches kids' ~3.5 with some axis of variation.** The
   categorical claim is wrong as stated, but you get a richer story:
   "LMs match human variability IF you vary the right axes." That's
   arguably more interesting.

Either result strengthens the SM2 follow-up.

## Design

Three variation axes, each with n ≈ 30 LM training runs:

1. **Seed-only** (n=30): identical data, arch, hyperparams; n random
   seeds. Establishes the "irreducible randomness floor" of LM
   variability.
2. **Data subsample** (n=30): identical arch and hyperparams; n
   different ~5–10M-token subsamples from CHILDES (or a synthetic
   alternative; see below). Analog of human "input quantity variation"
   meta-analyzed by Coffey & Snedeker.
3. **Architecture / hyperparam variation** (n ≈ 36 from a small grid):
   vary learning rate (3 values: 1e-4 / 3e-4 / 1e-3), model size
   (small / medium / large GPT-2), and seed within each cell. Analog
   of "constitutional / maturational variation" in kids.

Then `σ_κ_LM_total ≈ √(σ_κ_seed² + σ_κ_data² + σ_κ_arch²)` as a
generous upper bound. Compare to `σ_κ_kids ≈ 3.5`.

The decomposition itself is the interesting answer — it tells you
which axes LMs vary on, and you can map those to human-side variance
sources (Coffey & Snedeker's input share, processing share, etc.).

## The token-amount calibration matters

Empirical CDS rates from `input_estimation/validation_set.csv`:

| Rate (Sperry within-pool) | tokens/hr | tokens/mo (H=365) |
|---|---|---|
| –1 SD (low) | 880 | 321K |
| Mean | 1,500 | 547K |
| +1 SD (high) | 2,620 | 956K |

Months of child input equivalent to a given LM training-corpus size:

| Corpus size | Months at Sperry mean | Months at –1 SD | Months at +1 SD |
|---|---|---|---|
| 5M tokens | 9 mo | 16 mo | 5 mo |
| 10M tokens | 18 mo | 31 mo | 10 mo |
| 29M tokens (CHILDES) | 53 mo | 90 mo | 30 mo |

**Takeaways:**
- 5M ≈ early-CDI window
- 10M ≈ middle/late CDI window (median 24 mo CDI subject has heard ~13M tokens)
- 29M ≈ post-CDI for most kids

**Recommended primary training scale:** start with **10M** as the
default — it's the empirically-grounded "CDI-completion equivalent"
and is a meaningful developmental anchor. Drop to 5M if the pilot
shows it's needed for subsample independence (see next section).

## The subsample-independence problem

CHILDES has only ~29M tokens. n=30 disjoint 10M subsamples is
impossible (only ~3 fit). Overlapping subsamples bias `σ_κ_data`
downward.

**Two paths:**

- **CHILDES with overlapping subsamples**: simpler, but requires
  verifying that data-identity effects are small relative to seed
  effects at 10M scale (pilot step A below).
- **TinyDialogues (Feng et al. 2024)**: synthetic GPT-4-generated
  CDS-style data, much larger volume → easy n=30 disjoint samples.
  Cost: not real CDS, so the "input variation" you sample isn't
  quite the variation real kids experience. Weaker scientific framing
  but cleaner statistics.

Pilot results determine the choice. **Default preference is CHILDES**
for cleaner scientific framing; fall back to TD if pilot shows it's needed.

## Pilot (DO THIS FIRST)

8 LM training runs total. ~80 H100-hours. Cheap relative to the main experiment.

### Pilot A: data-identity vs seed sensitivity (4 runs)

Train at 10M CHILDES:
- 2 LMs on **disjoint** 10M CHILDES samples (split 29M → 3 chunks, use 2)
- 2 LMs on the **same** 10M CHILDES sample, different random seeds

Fit Chang & Bergen sigmoids on the ~611 CDI words for each. Compare:
- `σ_κ` across the 2 disjoint-sample LMs
- `σ_κ` across the 2 same-sample-different-seed LMs

Decision tree:
- If `σ_κ_disjoint >> σ_κ_seed`: data identity matters a lot at this
  scale → move data axis to TinyDialogues (or drop to 5M training).
- If `σ_κ_disjoint ≈ σ_κ_seed`: data identity barely matters →
  proceed with CHILDES overlapping subsamples.

### Pilot B: sigmoid stability at small scale (4 runs)

Same as A but at 5M training scale. Check how many CDI words have
well-fit C&B sigmoids (r² > 0.9 against the 4-param logistic):

- If <500/611 fit cleanly at 10M: scale up or accept lower coverage
- If <300/611 at 5M: 5M is too small; stay at 10M
- If 5M gives 400+ usable fits: 5M is a viable backup if subsample
  independence forces it

### Pilot deliverable

Short report (1–2 pages) with:
- Per-LM C&B sigmoid fit statistics (coverage, mean κ, slope SE)
- Cross-LM `σ_κ` estimates for the 4 conditions
- Recommendation: training scale (5M vs 10M) and data source
  (CHILDES vs TD)

## Main experiment (after pilot)

n=90 LM training runs at the scale and data source recommended by the
pilot. Three axes (seed, data, arch) × 30 runs each. Then:

1. Per-LM × per-word sigmoid fits via Chang & Bergen method
2. Per-LM `κ_LM` = mean of per-word slopes (or median, depending on
   noise structure)
3. Cross-LM `σ_κ_axis` for each axis
4. Combined upper bound `σ_κ_LM_total = √(Σ σ_κ_axis²)`
5. Compare to kids' `σ_κ ≈ 3.5`

**Stretch goal** (paper would be stronger with this): fit the SM2
hierarchical IRT-accumulator model directly to LM data
(LM-instance × word × checkpoint) instead of taking sample-variance
of point-estimates. Apples-to-apples with `log_irt_long` family.
Would require minor model adaptation.

## Compute and resources

**Estimate** (rough order of magnitude):
- GPT-2 small (~124M params) on 10M tokens: ~5–10 H100-hours per run
- Main experiment: ~90 runs × 10 = ~900 H100-hours
- Pilot: ~80 H100-hours
- Total: ~1000 H100-hours

**Marlowe** (Stanford's GPU cluster) is the natural target. Application docs:
- https://datascience.stanford.edu/marlowe/project-application-guide
- https://datascience.stanford.edu/marlowe/recharge

At Stanford-internal H100 rates (~$1–2/GPU-hour), ~$1500–3000 for the
full campaign. Apply for project allocation rather than ad-hoc usage.

**Storage**: ~200 training checkpoints per run × 500 MB × 90 runs ≈ 9 TB.
Marlowe has multi-TB scratch allocations; verify size limits in
application.

**Analysis side**: ~55K sigmoid fits at 4 params × ~200 datapoints
each = tens of CPU-hours total. Trivial.

## Existing references in this repo

For the new session to ground itself:

- **`reports/slides/standard_model.qmd`**: the talk deck. The
  "Connecting to LLMs" section (currently slides ~32–36) shows the
  current LM comparison via the C&B method.
- **`reports/model_explainer.tex`**: full mathematical specification
  of the IRT-accumulator. The `κ_i` parameter is the focus.
- **`journal/experiments.md`**: log of all SM2 experiments. The §
  with current LM comparison numbers is the baseline.
- **`input_estimation/`**: empirical input-rate data and validation
  set. `validation_set.csv` is the canonical reference for the
  child-input rate distribution.
- **`papers/feng_et_al_2026.pdf`** (if present): the
  CHILDES-trained GPT-2 work that this extends. Same architectures
  to use.
- **`papers/chang_bergen_2022.pdf`**: original per-word sigmoid method.
- **`papers/feng_tinydialogues_2024.pdf`** (if present): the TD
  synthetic data alternative.
- **`papers/input_estimation/`**: source literature for the
  empirical input rates above.

## Files this experiment would create

- `llm_variability/` — top-level directory for the experiment
- `llm_variability/training/` — LM training code (probably nanoGPT-style)
- `llm_variability/checkpoints/` — `<axis>/<condition>/<seed>/`
- `llm_variability/analysis/` — C&B sigmoid fits, per-LM κ tables
- `llm_variability/results/` — final cross-LM σ_κ estimates by axis

## Open questions to resolve in the new session

1. **Training infrastructure**: nanoGPT, custom Hugging Face, or fork
   Feng et al.'s repo? Probably easiest to fork Feng et al. if their
   code is available; matches the established baseline.
2. **Checkpoint frequency**: Chang & Bergen used ~200 checkpoints over
   training. Match this. Verify Feng et al.'s checkpoint frequency.
3. **Random seed control**: confirm GPT-2 training is fully
   reproducible at fixed seed (data ordering, dropout, init). Some
   training setups have non-determinism even at fixed seed (CUDA
   nondeterministic kernels). Want to either fix this or measure it
   as part of σ_κ_seed.
4. **Coffey & Snedeker mapping**: the human-side input share is 4–7%
   of variance in *raw vocab outcome* (R²). To compare with our
   LM σ_κ_data, we'd need to translate. Probably easier to compare to
   the SM2 model's σ_α² / σ_ξ² (efficiency-vs-input share) instead.
5. **Cross-architecture sanity check**: should we include any
   non-transformer LM in the population (e.g., a small LSTM) to make
   the architecture-variation axis broader? Or stick with
   GPT-2-variants only?
6. **Order of operations**: this experiment is post-SM2-paper. But
   the pilot can run in parallel with the SM2 paper writeup. Decide
   whether to delay pilot until SM2 lands.

## Relationship to in-flight SM2 GCP fits

The SM2 GCP family fits are running independently on `sm2-fit-01`
and `sm2-fit-02`. They use the s=6 pinning and won't be affected by
this experiment. Expected to land in ~3 days from 2026-05-19. After
they land, the LM σ_κ comparison can use refreshed kid-side numbers.

The σ_κ_kids ≈ 3.5 figure used here is from s=0.5 fits and may shift
slightly with s=6 refits (probably σ_κ stays similar in magnitude;
κ_pop drops from ~10 to ~5–7). Check `fits/summaries/` for the
`long_no_freq_slopes_si_signed.summary.rds` once refits land, and
update the comparison target accordingly.

## Skills to load

A fresh session on this task should look at:

- `~/.claude/skills/sherlock-stan-fitting.md` (Sherlock has GPUs but
  shares Mike's fairshare with the SM2 campaign — Marlowe is preferred)
- Marlowe-specific skill doesn't exist yet; this experiment could
  generate one as a byproduct.

The standard_model_2-project skills don't apply directly here (they're
about Stan, not LM training), but the general `remote-vm-gotchas`
patterns do (detached launches, disk-fill traps).

## Suggested first session steps

1. Read this doc + Feng et al. 2026 paper
2. Verify CHILDES corpus access and exact token count (~29M is
   approximate; could be ±5M depending on filtering)
3. Verify TinyDialogues access and size
4. Pull Feng et al.'s training repo if available; otherwise pick
   training framework
5. Apply for Marlowe project allocation with 80 H100-hour pilot
   request (this can be approved fast, often within a week)
6. Run pilot A + B (8 runs)
7. Write pilot report
8. Discuss results with Mike before committing to main experiment
