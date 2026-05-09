# Handoff: extract CDI word-acquisition sigmoid slopes from Feng et al. (2026) CHILDES-trained GPT models

## Context (read this first)

This is a side quest off the **Standard Model 2** project (`github.com/langcog/standard-model-2`). The main project fits a Bayesian psychometric model of early word learning (Wordbank CDI data) and characterizes two signatures: **acceleration** (per-child scaling exponent $\kappa_i$ on cumulative input) and **variability** (between-child spread).

The current external-facing deck (`outputs/slides/standard_model_external.pdf`) ends with a **comparison to LLM word-acquisition slopes** drawn from Chang & Bergen (2022). On a common axis (sigmoid slope of $P(\text{word acquired})$ on $\log$-experience, in nat-log units), kids cluster at $\kappa_i \sim 10$ while LMs cluster near $1$ — about a 10× gap.

**Issue:** Chang & Bergen trained their LMs on BookCorpus + WikiText (adult written text). That conflates two factors:

1. The structural difference in scaling between humans and LMs.
2. The input-distribution mismatch (adult text vs. child-directed speech).

Feng et al. (2026) train GPT-style models on CHILDES (paper PDF in `papers/feng_etal_2026.pdf`). These models are input-distribution-matched to what kids actually hear. **Replicating Chang & Bergen's per-word sigmoid analysis on Feng et al.'s CHILDES-trained models is the next empirical step.**

## Your task

Produce a new version of `outputs/figs/longitudinal/chang_bergen_slope_comparison.png` that adds Feng et al.'s CHILDES-trained models alongside (or in place of) Chang & Bergen's BookCorpus-trained models. Same axis, same sigmoid-slope statistic, same CDI vocabulary set.

**Hypothesis to evaluate:** the structural-difference hypothesis predicts the kid-vs-LM gap will largely *remain* even when LM training is CDS-matched. If it shrinks substantially, the input-distribution explanation deserves more weight.

## What you need to read first

1. **`papers/feng_etal_2026.pdf`** — the source paper. Figure out:
   - What models do they release? (Architecture, parameter counts, training tokens.)
   - Do they save training-step-by-step checkpoints / surprisal trajectories? Or only end-of-training models?
   - If they don't release per-step surprisal, we'll need to retrain following their setup.
2. **`papers/chang_bergen_2022.pdf`** — the original LM-word-acquisition pipeline we're extending.
3. **`outputs/chang_bergen_derivation.tex`** — derivation of why per-word $1/\text{ParamScale}$ on $\log$-steps is the right statistic. Also documents the unit conversion (`× ln(10)`) we use.

## What's already in the repo

- **`data/chang_bergen_2022/{bert,bilstm,gpt2,lstm}_sigmoids.txt`** — Chang & Bergen's per-(LM, word) sigmoid fit parameters (`ParamUpper, ParamLower, ParamXmid, ParamScale`), ~600 CDI words × 4 architectures. Same format you'll want to produce for Feng et al.'s models.
- **`model/scripts/chang_bergen_comparison.R`** — current pipeline that reads those files, computes $0.434/\text{ParamScale}$ (logit per nat-log step), and plots against our $\kappa_i$ posterior. Read this to understand the target output format.
- **`outputs/figs/longitudinal/chang_bergen_slope_comparison.png`** — the existing figure your work would update.
- **`fits/summaries/long_no_freq_slopes.draws.rds`** — our English M_best posterior. The kid-side $\kappa_i$ distribution is sampled from this in the existing comparison script (lines 88–110 of `chang_bergen_comparison.R`).

## What you need to do

### Phase 1 — figure out what Feng et al. provides

Open `papers/feng_etal_2026.pdf`. Extract:
- Model architectures + sizes
- Training data scale (CHILDES tokens used)
- Whether per-step checkpoints are released (HuggingFace? GitHub? other?)
- Whether they release per-word surprisal trajectories during training (unlikely but check)

Decide between two paths based on what they release:

**Path A (cheap):** They publish per-step model checkpoints OR per-step surprisal logs.
- Pull the checkpoints / logs.
- For each (model, word) pair, fit Chang & Bergen's 4-parameter logistic to the surprisal trajectory.
- Produce `data/feng_2026/<model>_sigmoids.txt` in the same format as `data/chang_bergen_2022/`.

**Path B (expensive):** Only final-trained models are released.
- Replicate their training pipeline, logging per-step surprisal for the CDI vocabulary.
- Note: Chang & Bergen logged at ~200 training steps, sampling more heavily early.
- Total training cost depends on model scale.

If Path B is the only option, **stop and report back** before starting training — we should discuss whether to proceed.

### Phase 2 — extract per-word sigmoid slopes

For each Feng et al. model:

1. Pull the **same 611 CDI words** Chang & Bergen used (their list is implicit in `data/chang_bergen_2022/bert_sigmoids.txt`'s `Token` column).
2. For each word $w$ in each model: fit
   $$S_w(x) = \mathrm{ParamLower}_w + \frac{\mathrm{ParamUpper}_w - \mathrm{ParamLower}_w}{1 + \exp((x - \mathrm{ParamXmid}_w)/\mathrm{ParamScale}_w)}$$
   where $x = \log_{10}(\text{training steps})$ and $S_w$ is mean per-word surprisal at that step (averaged over $\sim$512 occurrences in held-out data, per Chang & Bergen).
3. Compute slope on natural-log experience: `slope_natural = 0.434 / ParamScale_w`.
4. Save outputs in `data/feng_2026/<model>_sigmoids.txt` with the same columns as Chang & Bergen.

### Phase 3 — produce the comparison figure

Adapt `model/scripts/chang_bergen_comparison.R` to add Feng et al.'s models. Two options:

- **Option A:** Add Feng models alongside Chang & Bergen's (so the figure shows kids + 4 BookCorpus LMs + N CHILDES LMs).
- **Option B:** Replace Chang & Bergen's models with Feng's, since CDS-matched is the cleaner comparison.

Mike will likely want both versions for the deck. Save as:
- `outputs/figs/longitudinal/feng_slope_comparison.png` (Feng-only)
- `outputs/figs/longitudinal/feng_chang_bergen_slope_comparison.png` (combined)

### Phase 4 — write a brief report

In `outputs/feng_evaluation_report.md`, write 1–2 pages covering:

1. What models you used / which path (A or B) you took.
2. Per-model median + IQR of `slope_natural`, in the same format as the existing Chang & Bergen table in `outputs/chang_bergen_derivation.tex`.
3. The headline: does the kid-vs-LM gap persist with CDS-matched training?
4. Caveats specific to Feng's setup (vocab size differences, tokenization, CHILDES preprocessing, sample sizes).

## Things to watch for

- **Tokenization mismatch.** Chang & Bergen restricted to the 611 of 651 CDI words their tokenizers covered as single tokens. Feng et al.'s tokenizer might differ; figure out what fraction of CDI items are single tokens.
- **CHILDES is much smaller than BookCorpus.** ~5–10M CHILDES tokens vs. 25M+ sentence pairs in Chang & Bergen. Their scaling won't reach the same training-step range; check whether sigmoid fits are well-conditioned at lower step counts.
- **What's the unit on the x-axis?** Chang & Bergen use $\log_{10}(\text{steps})$. Feng might use steps, tokens, or epochs. Convert consistently to natural-log $\text{steps}$ (or to natural-log $\text{tokens}$, but then convert ours too — currently kids' axis is in $\log$-time-as-experience-proxy, which Chang & Bergen also implicitly does via training steps).
- **Multiple checkpoints vs. one.** Per-word sigmoid fitting needs ~10+ training-step samples to converge. If Feng only released a handful of checkpoints, the sigmoid fits will be unstable.
- **Random seeds / replicates.** If Feng et al. release multiple seeds per model, fit each separately and treat each (model, seed) as a replicate. This is the closest thing on the LM side to between-instance variance, which we currently report as zero.

## Deliverables

- `data/feng_2026/<model>_sigmoids.txt` (per-word sigmoid params, one file per Feng model)
- `outputs/figs/longitudinal/feng_slope_comparison.png` (and/or combined version)
- `outputs/feng_evaluation_report.md` (1–2 page report)
- A small commit to the standard-model-2 repo with all of the above, plus any fitting scripts under `model/scripts/`

## What NOT to do

- Don't refit our kid-side model. The $\kappa_i$ posterior is already extracted and stable.
- Don't change the unit conversion. Stick with $0.434/\text{ParamScale}$ when ParamScale is in $\log_{10}$ units. If Feng's $x$ axis is in different units, document the conversion.
- Don't modify any scripts in `model/scripts/` other than to add new ones for Feng-specific fitting and plotting. The existing pipeline works.

## Open questions to resolve early

1. Does Feng et al. release training-step checkpoints, or only final models?
2. What's the largest model they trained? Smallest?
3. Did they ablate over CHILDES subsets (e.g., by age range of the children in CHILDES)?

If you can answer these from a quick read of the paper's methods section, do that first. If not, get back to Mike before committing to either Path A or Path B.

## Contact

When you have results, post them on the standard-model-2 repo with a PR or branch. Mike will fold them into the slide deck and the eventual paper.

---

*This handoff was generated by Claude (Standard Model 2 main session) on 2026-05-09.
The session that produced this lives at `~/.claude/projects/-Users-mcfrank-Projects-standard-model-2/`.*
