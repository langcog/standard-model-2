# Marchman data issue — glmer cross-check confirms it's an extraction artifact (2026-07-17)

**TL;DR (MCF):** The Marchman dataset looks anomalous in the raw glmer ladder
because of a **data-extraction defect, not a model-structure problem**. The
three targeted cleaning steps already implemented in
`studies/bayes_long/00_prepare_bundles.R` recover the data; a standard-glmer
refit on the cleaned frame confirms the anomaly is a data artifact. **The
bayes_long bundles already apply this cleaning**, so the Bayesian Marchman fit
should already be on clean data — see "Implication for the Bayesian session."

This note is a handoff for the Bayesian-fits session.

## The defect (raw Wordbank → glmer extraction)

`studies/glmer_ladder/01_extract_one.R` keys children on Wordbank `child_id` and
keeps kids with ≥2 administrations. For Marchman this breaks two ways:

1. **`child_id` fails to link a child's WG and WS administrations.** A two-wave
   child (WG at ~8–16 mo, WS at ~16–30 mo) is split into two *unlinked*
   single-form records. Each then has only 1 admin, so the ≥2-admin longitudinal
   filter **deletes them wholesale** — and the entire WG (younger-age) arm is
   lost. Symptom: raw Marchman spans only **age 16–30 mo** (the WS window) with
   **314 kids**.
2. **Mis-keyed WG comprehension craters.** Some Marchman WG records have
   comprehension mis-scored as production: vocabulary spikes then craters to ~0,
   which a "child cannot un-produce words" QC rule should drop.

## The fix (already in `bayes_long/00_prepare_bundles.R`)

Three targeted steps, all data-cleaning, no model change:

1. **Child key = `dataset::study_internal_id`** (not Wordbank `child_id`) — relinks
   WG↔WS admins into one child. Requires
   `get_administration_data(..., include_study_internal_id = TRUE)`.
2. **Item harmonization by `uni_lemma`** — a WG item and WS item sharing an
   unambiguous `uni_lemma` are the same latent item (handles in/inside → inside/in).
3. **QC crater drop** — drop a child whose final vocabulary collapses far below its
   running peak (`end < peak·(1−0.25)`, with peak≥0.10 and drop≥0.05 floors).

## Evidence: standard-glmer refit on cleaned data

I ported those exact three fixes into a clean glmer extraction (at MIN_ADMINS=2,
matching the raw pipeline's ≥2-admin filter, so raw-vs-clean isolates *only* the
cleaning) and refit the standard `glmer_ladder` model unchanged.

**Extraction — the WG arm comes back:**

| | raw | clean |
|---|---:|---:|
| kids | 314 | **2060** |
| age range | 16–30 mo | **8–30 mo** |
| obs | 485,520 | 2,397,435 |
| QC craters dropped | — | 134 |

Child-key integrity checked: no NA-keyed mega-child, max 4 admins/kid (median 2),
real subject codes; **934 kids genuinely span WG→WS** (≤14 to ≥20 mo). The 314→2060
jump is legitimate — the raw defect was deleting two-wave kids entirely, not just
trimming the WG arm.

**B_log (fixed κ, `produces ~ 1 + log_age + (1|item)`) — slope normalizes:**

| dataset | a0 | log_age slope (≈κ) |
|---|---:|---:|
| Marchman **raw** | 24 | **8.08** (steepest of English) |
| Marchman **clean** | 20 | **6.74** (now in-band) |
| Thal | 16 | 8.05 |
| Smith | 22 | 7.21 |

The raw steepness was an artifact of fitting only the fast-accumulating 16–30 mo
window plus craters; cleaned, Marchman sits squarely in the English band.

**Raw D_log intercept anomaly:** raw Marchman `D_log` fixed intercept ≈ **−0.02 at
a0=24 mo** (≈50% production at 24 mo) vs Smith −1.65 at a0=22 — implausibly high,
consistent with the craters inflating production. (a0 differs, so compare at a
common age; the clean B_log intercept −1.51 at a0=20 is far more sensible.)

## Open item — clean D_log (M_best) not fit locally

The M_best model `produces ~ 1 + log_age + (1 + log_age | child) + (1 | item)`
did **not** complete locally: on the clean frame (2.4M rows, 2060 kids, random
intercept+slope) `glmer`/bobyqa ran **5+ hours at 100% CPU without converging**
and was killed. The raw D_log fit took minutes — the clean data is 5× the rows
and 6.5× the kids, and random-slope models scale badly in the kid dimension. If a
glmer M_best cross-check is wanted, run it on the **Sherlock `glmer_ladder` slurm
pipeline** (see `studies/glmer_ladder/README.md`), not locally. It is not needed
to establish the point above.

## Implication for the Bayesian session

`studies/bayes_long/00_prepare_bundles.R` **already applies all three fixes**, and
`fits/bayes_long/bundle_marchman_a3.rds` was built from it (with MIN_ADMINS=3). So:

- The Bayesian Marchman fit is **already on cleaned data** — it does *not* share the
  raw-glmer extraction defect. If Marchman looked odd in a Bayesian fit, check that
  the bundle in use came from `00_prepare_bundles.R` (study_internal_id key + QC
  drop), not from the glmer `data_marchman.rds`.
- The takeaway from the glmer cross-check: **cleaning, not model structure, is the
  lever** for Marchman. The accelerating-accumulator structure is fine; the raw
  glmer anomaly was extraction.
- Note the threshold difference: the glmer cross-check used MIN_ADMINS=2 (to match
  the glmer pipeline); the a3 Bayesian bundle uses MIN_ADMINS=3.

## Artifacts / repro

- `studies/glmer_ladder/01b_extract_marchman_clean.R` — clean extraction (the three
  fixes at MIN_ADMINS=2), writes `fits/glmer_ladder/data_marchman_clean.rds`.
- `fits/glmer_ladder/summary_marchman_clean_B_log.csv` — clean B_log fit (done).
- Raw comparators: `fits/glmer_ladder/summary_marchman_{B_log,D_log}.csv`,
  `summary_{thal,smith}_{B_log,D_log}.csv`.
- Repro:
  ```
  Rscript studies/glmer_ladder/01b_extract_marchman_clean.R
  Rscript studies/glmer_ladder/02_fit_one.R marchman_clean B_log   # fast
  # D_log: run on Sherlock, not locally
  ```
