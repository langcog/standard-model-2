# Scoping: joint processing + input model (maximal-precision Panel E)

**Status:** scoping only (2026-06-12). No fitting yet. Decision pending MCF.

## Goal

Replace the two *separate, noisy* Panel-E direct-measurement points (io input
model on 4 datasets; proc RT model on 3 datasets) with **one joint model** that
fits input **and** processing on **all** candidate datasets simultaneously.

Two payoffs:

1. **Precision.** Pool every child that has *either* channel; partial-pool across
   studies → tighter input and processing variance shares.
2. **Mechanism (the real prize).** Separate io and proc models cannot
   disentangle input from processing because the two are correlated (Weisleder
   "rich get richer": input → faster processing → more learning). A joint model
   with *both* per-child predictors estimates input's effect on efficiency
   **controlling for** processing, and vice versa — i.e. whether input acts on
   vocabulary **directly** or **via** processing speed. The current panel can't
   answer this.

## Key finding: the Stan model already supports this

`model/stan/log_irt_long_proc_dp.stan` uses **ragged channel indexing**:

- LWL/RT side: `N_lwl`, `lwl_to_child[N_lwl]` — a child appears only if it has RT.
- Input side: `V`, `rec_to_child[V]`, `z_lena[V]` — a child appears only if it
  has observed LENA/headcam input.

So **input-only** children (BabyView, SEEDLingS) and **RT-only** children
(FM2012) are handled *natively* — they simply contribute to whatever channels
they have. **No masking change to the core model is needed.**

The single place that currently excludes input-only kids is the **bundle**, not
the model: `prepare_proc_dp_bundle.R:95`

```r
kids <- intersect(unique(cdi$lab_subject_id), unique(lwl$lab_subject_id))  # CDI ∩ RT
```

This must become `CDI ∩ (RT ∪ input)`.

## Dataset × channel inventory

| dataset | CDI | input | LWL RT | role |
|---|:--:|:--:|:--:|---|
| AM2018 (adams_marchman_2018) | ✓ | ✓ LENA | ✓ | **both** (66 kids) — identifies separation |
| FMW2013 (fmw_2013) | ✓ | ✓ LENA | ✓ | **both** (31 kids) — identifies separation |
| FM2012 (fernald_marchman_2012) | ✓ | imputed | ✓ | RT-only (119 kids) — RT precision |
| BabyView (Long 2024) | ✓ | ✓ headcam | — | input-only (~44 kids) — input precision |
| SEEDLingS (Egan/Bergelson) | ✓ | ✓ LENA | — | input-only (~51 kids) — input precision |

Identification of the input|processing split rests on the **97 both-channel
children** (AM2018 + FMW2013). FM2012 adds RT precision; BabyView + SEEDLingS
roughly **double the input channel** (97 → ~190 input-having kids). Total ≈ 320
children, S = 5 studies.

## The work = bundle assembly (the model is ~ready)

`prepare_proc_dp_bundle.R` already assembles CDI + LWL-RT + observed-LENA for the
3 LWL datasets. Extending it:

1. **Relax the child set** to `CDI ∩ (RT ∪ input)` (line 95).
2. **Add BabyView + SEEDLingS** to the CDI side and the input side. Their CDI +
   per-child input already live in the `io_pooled` bundle
   (`fits/io_pooled_subset_data.rds`) and the per-dataset io bundles
   (`io_{babyview,seedlings,...}_subset_data.rds`). Need to:
   - reconcile child IDs across the two sources (io uses an internal `ckey`; the
     per-dataset io bundles carry `subject_id` → map to `lab_subject_id`),
   - reconcile the **two CDI sources** for the overlapping AM2018/FMW2013
     (proc_dp's `stanford_cdi_items_long.csv` vs io_pooled's CDI) — pick one
     canonical CDI per child; they should agree post-QC but must be checked.
3. **Pool input across channels.** Input is LENA for AM2018/FMW2013/SEEDLingS and
   **headcam** for BabyView. The bundle already **standardizes input within
   study** (z_lena), so the *units* are comparable. The open issue is the
   **measurement-noise pin** `sigma_lena` (currently one scalar from the
   AM2018/FMW2013 16/18-mo replicates). Headcam noise ≠ LENA noise.
4. **Consistent child index `ii`** across CDI / RT / input arrays.

## Likely small Stan change

The model pins a **single** `sigma_lena` (input measurement noise, sd units).
With 5 studies and two instruments (LENA vs headcam), this should become
**per-study**: `vector[S] sigma_lena` indexed by `study_of_recording`. That's a
~3-line change (scalar → vector + index). Everything else (the ragged channels,
the D'0–D'3 ladder toggles, the variance partition) is unchanged.

`sigma_lena` per study needs an empirical value each: AM2018/FMW2013 from the
existing replicate analysis; BabyView/SEEDLingS from their own test-retest if
available (else a conservative pin + sensitivity).

## Open modeling decisions (for MCF)

- **σ_r pin still applies.** ξ_i = μ_r + σ_r·z_r + β_ξ·rt0 + log_α: σ_r scales the
  per-child input deviation, observed (via z_lena + sigma_lena) for most kids,
  imputed (z_r ~ N(0,1)) for FM2012. So the joint model **still depends on the
  σ_r pin** — feed it the apples-to-apples σ_r ≈ 0.44 we just derived (the
  validation pins confirm the curve). The input channel is mostly *observed*
  here, so it's far less σ_r-sensitive than the pure imputed EN/NO analysis.
- **FM2012 input:** keep imputed, or drop FM2012 (it adds only RT and may pull
  the separation if its imputed input is mis-scaled). Recommend keep + check.
- **Item subsample** must be reconciled to one shared J across all 5 datasets
  (proc_dp uses a class×difficulty stratified subsample; io has its own).
- **Rung to report:** the D'1 architecture (rt0→ξ selected in entry 33), now
  re-estimated on the larger sample; optionally re-run the D'0–D'3 ladder.

## Compute

Small: ≈320 children, vocab N on the order of 1e5 (item subsample), N_lwl ≈ 1k,
V ≈ 190. This is **much** smaller than the Wordbank longitudinal fits — expect
well under an hour per rung on one node. No disk/LOO concerns at this scale
(LOO still skippable; we only need the variance-partition scalars).

## Risks

- **Child-ID reconciliation** across the io and proc_dp sources is the main
  failure point (silent mis-joins → duplicated/dropped kids). Validate counts
  per dataset against the table above before fitting.
- **Two CDI sources** for AM2018/FMW2013 must not double-count or disagree.
- **Headcam vs LENA noise** — per-study `sigma_lena` mitigates, but BabyView's
  headcam input is a sparser proxy; watch its leverage on the input share.
- The separation (direct vs processing-mediated) is identified by only **97**
  both-channel kids — adding input-only data tightens the *marginal* input share
  but **not** the separation. State this honestly.

## Proposed sequence (when greenlit)

1. Extend `prepare_proc_dp_bundle.R` → `prepare_joint_io_proc_bundle.R`
   (union child set; add BabyView/SEEDLingS CDI+input; per-study sigma_lena).
2. Validate bundle counts vs the inventory table.
3. `sigma_lena` scalar → `vector[S]` in `log_irt_long_proc_dp.stan` (or a new
   `_joint` variant to avoid disturbing the committed proc_dp).
4. Fit D'1 (+ optional ladder) on GCP with σ_r = 0.44 pinned; extract the
   input/processing × efficiency/acceleration partition.
5. Replace Panel E's two separate points with the joint partition.
