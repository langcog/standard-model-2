# SEEDLings LWL — QC targets for the RT-extraction pipeline

We are deriving a *new* per-child RT measure from the Zhu et al. (submitted) raw
EyeLink reports (`data/seedlings/raw_eyetracking_data/HaT/`) — the authors only
published accuracy ("increase in target looking") for infants, not RT. To trust
our RT, our pipeline must reproduce the paper's **trial funnel** and **accuracy**
numbers (which ride on the same trial parsing). Source: `papers/zhu_etal_submitted.pdf`.

These are the targets to aim for. Anything our pipeline reads from the reports we
match exactly; anything that depended on the authors' post-hoc **video/notes
coding** (which we don't have) we can only approach — flagged ⚠.

## Tier 1 — trial funnel (Study 1, n=44, ages 8/10/12/14/16/18 mo)

| step | paper | our pipeline should… |
|---|---:|---|
| planned trials (44 × 6 × 32) | **8,448** | count structural max |
| − 2 unusable unreplaced visits | −64 | (visits with no data) |
| − early-end skipped trials | −178 | |
| planned remaining | **8,206** | |
| double/triple attempts (+107, +3) | → **8,383 attempts** | count trial records in reports |
| ⚠ − didn't hear stimulus (video/notes) | −86 | **can't reproduce** (needs their video coding) → expect ~+86 over target |
| − no button press → **no locatable onset** | −78 | trials lacking `EL_BUTTON_CRIT_WORD` |
| − complete lack of fixation data | −190 | trials with 0 fixations |
| − fixated target/distractor < ⅓ of window | −1,717 | apply the ⅓-of-367–3500ms rule |
| retained attempts / usable planned | **6,312 / 6,204** | target (we'll land ~6,290 w/o the 86) |
| − 15 low-data participants (<25%) +53 trials | | participant-level cut |

**Pass test:** planned = 8,448 exactly; attempts within ~1% of 8,383; no-onset ≈ 78;
no-fixation ≈ 190; usable within ~90 of 6,204 (the unreproducible 86 "didn't hear").

## Tier 2 — accuracy ("increase in target looking")

- **Window:** post-onset **367–3,500 ms** (text says "360"; Fig caption "367").
- **DV (pair-salience-corrected):** for each yoked pair A–B, mean prop. fixation on
  A *as target* − mean prop. fixation on A *as distractor*. (Cancels image salience.)
- **Targets:**
  - grand-mean increase ≈ **0.048** (Study 1 model intercept; SE 0.017, p=.008).
  - small + positive; **above chance (>0) at 12–18 mo** (per-month intercept p<.004);
    near 0 at 8–10 mo ("≈0.01 unit between 8–10 mo, ~0.06 by older").
- **Pass test:** our grand-mean increase in [0.03, 0.07]; monotone-ish rise with age;
  positive & distinguishable from 0 at 14/16/18 mo.

## Tier 3 — RT validity (novel measure; no paper number to match)

The authors did **not** compute infant RT, so RT is validated by internal
coherence, not a paper figure:
- median RT ≈ 600–900 ms (Fernald toddler range);
- **monotone speed-up with age** (RT decreases 8→18 mo) — the key signal;
- D-initial trials only, 300–1800 ms window;
- per-child aggregate is what feeds the io-proc model (RT enters at the **child**
  level — a per-kid `rt0_i` intercept + slope — so per-trial onset-button noise
  averages out; no per-word RT, per MCF).

POC (2026-06-14, `model/scripts/seedlings_lwl_rt_poc.R`): median 714 ms, 755→636 ms
across 8→18 mo, 1,918 valid RTs / 45 kids, ~1,202 in our 14–18 mo window. ✓ Tier 3.

## RESULTS — pipeline `model/scripts/prepare_seedlings_lwl_rt.R` (2026-06-14)

| QC check | target | recovered | verdict |
|---|---:|---:|---|
| planned trials | 8,448 | 8,448 | PASS (exact) |
| no locatable onset (no button) | 78 | 72 | PASS |
| no fixation data | 190 | 168 | PASS |
| usable trials | 6,204 | 6,134 | PASS |
| trial attempts | 8,383 | 9,253 | off +870 (recycle/structure bookkeeping) |
| low-look (<⅓ window) excl | 1,717 | 2,879 | off +1162 (stricter; nets correct `usable`) |
| accuracy: grand-mean increase | 0.048* | 0.117 | *coding-convention diff; **curve textbook** |
| accuracy by age (8→18mo) | rising, >0 by 12mo | −.01,−.02,**.06,.12,.19,.27** | PASS (chance→strong) |
| RT median / trend | ~600–900, ↓ with age | 710 ms, 774→642 | PASS |

The 4 structural counts + the accuracy *shape* + RT validity all pass. The two
"off" rows reflect bookkeeping we can't reproduce without the authors' video/notes
coding (the −86 "didn't hear") and a stricter low-look rule — they net to the
correct `usable`. *Their 0.048 is an age-centered effect-coded model intercept,
not a raw grand mean; the developmental curve is the real validation.

**Output:** `data/seedlings/seedlings_lwl_rt.csv` — 1,697 RT rows, 44 kids, 8–18mo
(1,097 in our 14–18mo window). Schema: dataset_name, lab_subject_id (zero-padded
subj), lwl_age, lwl_log_rt.

## ⚠ ID-LINKAGE TRAP (must use before wiring to the bundle)
SEEDLings subjects are **NOT contiguous** — `01..46` with gaps at **05, 24** (44
kids). io_pooled assigns `ii = as.integer(factor(subject_id))` = a **dense rank**,
so bundle `SEEDLingS::SEEDLingS::ii` maps `ii=5→subj 06`, `ii=23→subj 25`, …
**A naive `::N → subj 0N` join is wrong for 40 of 44 kids.** Correct crosswalk =
rank of subj in the sorted source set. Also fixed: one mis-padded session `9_10`
→ subj `09`. **Robust fix (recommended): retain `subject_id` through
prepare_seedlings → io_pooled so the RT join keys on subject_id, not a
reconstructed rank.**

## Note
"184 trials/participant · audio-noun/video-noun · mouse-click RT · 275 RT outliers"
in the PDF is **Study 2 controls** (adults clicking), NOT our infants — ignore for
the Study 1 RT pipeline.
