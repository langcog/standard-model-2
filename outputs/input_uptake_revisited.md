# Input-uptake revisited: data + model audit before the γ analysis

**Status:** planning. Nothing here is run yet. The goal is to clean up
the input-observed ("IO") datasets and models *before* testing whether
input rate predicts the per-child growth slope (the "input-on-slope"
coefficient γ).

## Why we're here

The accumulator puts input rate `log r_i` in the **intercept** only —
more input is a level shift, not a change in the growth exponent κ. The
alternative (bootstrapping / Matthew effect) says more input also makes
you *accelerate faster*: `log r_i` should predict the per-child slope
deviation ζ_i, i.e. γ > 0. This is only identifiable where input is
*measured* (the IO datasets), not in CDI-only data where `r_i` is
latent.

**Reduced-form pre-check (already run, `model/scripts/input_slope_check.R`):**
regressing posterior ζ_i on posterior `log_r_true_i` across kids gave

| Dataset | N | γ (ζ-SD per +1 SD log r) | cor(log r, ζ) | P(γ>0) |
|---|---|---|---|---|
| BabyView | 20 | +1.19 [+0.72, +1.66] | +0.43 | 1.00 |
| SEEDLingS | 44 | +0.37 [−0.09, +0.83] | +0.09 | 0.94 |

Suggestive positive signal in BabyView, not replicated in SEEDLingS.
Tiny N, the two disagree, and SEEDLingS' input signal is weaker
(`cor(log r, ξ)` = 0.18 vs BabyView 0.53). **Not trustworthy until the
data + models below are cleaned up.**

---

## Part 1 — Data audit & fixes (must happen first)

### 1.1 The datasets (CONFIRMED crosswalk, 2026-05-28)

Five datasets. Going forward, refer to the three Fernald/Marchman-lab
cohorts by their reduced paper titles: **AM2018, FMW2013, FM2012**.
Retire "Stanford" / "Peekbank" as dataset names (Peekbank is the *LWL
data source*, not a study).

| Name | Internal | Lab IDs | CDI | Input (LENA) | LWL (Peekbank) | Roles |
|---|---|---|---|---|---|---|
| **BabyView** | — | — | WG+WS, 9–32 mo, 20 kids | head-cam AWC | — | IO |
| **SEEDLingS** | — | — | WG, 6–18 mo, 44 kids | LENA daylong | — | IO |
| **AM2018** | TL3 / totlot3 | 11xxx | WG (13–18) + WS (20–32) | LENA @ 16,18 mo (66) | `adams_marchman_2018` (711) | **IO + processing** |
| **FMW2013** | TLO | 20xxx | WG@18 + WS@24 + WS@30 | LENA @ 18 mo (51) | `fmw_2013` (178) | IO |
| **FM2012** | TL2 / totlot2 | — | WG+WS, 409 admins | **none** | `fernald_marchman_2012` (679) | **processing only** |

- **4 IO datasets** (CDI + measured input): BabyView, SEEDLingS,
  AM2018, FMW2013 → the γ analysis.
- **2 processing datasets** (CDI + LWL): AM2018, FM2012 → the separate
  processing study. FM2012 has no LENA, so it's processing-only.
- AM2018 is the only dataset with all three channels.

CDI source files (Marchman lab, `data/peekbank/`): `TL3_compiled_WS.csv`
+ `TL3_compiled_WG.xlsx` (AM2018); `TLO_18m_WG.xlsx` + `TLO_24_WS.xlsx` +
`TLO_30m_WS.xlsx` (FMW2013); `totlot2/TL2_W{G,S}_compiled.xlsx` (FM2012).

**Parsing notes for the new files:** the TLO files carry the CDI age in
the *filename* (no age column in the sheet) and include some
`study == "misc"` rows to filter out (keep `tlo`). TL3 WG has a real
`age_cdi` column (13–18). The existing `parse_stanford_cdi.R` →
`stanford_cdi_items_long.csv` pipeline must be extended to ingest these
(it currently only emits totlot2 + totlot3, WG only for totlot2).

### 1.2 The big fix: use WG + WS everywhere

`prepare_stanford_linked.R:43` restricts to WS only ("All TL3 admins
are WS; restrict for purity"), producing the narrow 20–27 mo window
that crushed κ and inflated σ_ζ (see experiments §17, and the
glmer_ladder Spanish/Finnish form-window problems). **This was a
mistake.** Every dataset's CDI side should combine WG + WS at the item
level (production only), exactly as the glmer_ladder extractor now does
(`glmer_ladder/01_extract_one.R`):
- AM2018 (TL3) genuinely has WG+WS across 16–30 mo — pulling both
  widens the window from 7 months to ~14 and should recover a sane κ.
- SEEDLingS is WG-only by age (6–18 mo) — fine, but confirm no kids
  aged into WS.
- BabyView already does WG+WS — no change.
- Watch the per-form item-count ceiling issue: a kid's vocab is bounded
  by the items on the form they took. For *fitting* this is fine
  (item REs handle it); for *plotting* predictions, sum over a single
  reference form's items (glmer_ladder lesson).

### 1.3 Deliverable for Part 1

Re-prepared bundles, one per dataset, all WG+WS, with three optional
channels marked present/absent: `cdi`, `input`, `lwl`. Naming:
`io_babyview`, `io_seedlings`, `io_am2018`, `io_fmw2013`. A short
`data/README` crosswalk from Peekbank dataset names → study.

---

## Part 2 — Model architecture decisions

### 2.1 The reactivity / observation multiplier (needs review)

The IO model has `beta_react`: observed input is modeled as
`log_r_obs ~ N(log_r_true + beta_react, sigma_within)`. Intended to
absorb the gap between *observed* tokens (during a recording) and
*habitual* daily input. Fitted values: BabyView β_react ≈ 0.28 (head-cam,
short reactive sessions), SEEDLingS ≈ 0 (LENA daylong, less reactive).

**Questions to resolve:**
- Is a single additive log-scale shift the right structure, or should it
  be per-dataset (different recording regimes) and/or interact with age?
- Does β_react trade off against the γ estimate? If reactivity were
  *rate-dependent* (talkative families both higher-input AND more
  reactive), a mis-specified β_react could leak into γ. Worth a
  sensitivity check: fit with β_react fixed at 0 vs free, see if γ moves.
- For head-cam (BabyView) the sampling is also *spatially* partial (only
  what the camera sees), a different bias than reactivity. Consider
  whether the multiplier is conflating two things.

### 2.2 Share word difficulties δ_j across datasets — **recommended**

These are the **same English CDI items** across BabyView, SEEDLingS,
AM2018, FMW2013 (and the 1840-kid Wordbank longitudinal). Right now each
IO fit re-estimates ~200 δ_j from 20–60 kids — data-starved, and it
wastes the small samples' information on nuisance parameters that the
big data already pins precisely.

**Two-stage approach (recommended first pass):**
1. Estimate δ_j precisely from the big EN longitudinal fit (1840 kids).
2. In each IO fit, fix δ_j (or use tight informative priors centered on
   stage-1 values). The 20–60 IO kids then spend all their data on the
   per-child input / efficiency / slope params — exactly what we need
   for a clean γ.

This is modular, debuggable, and captures most of the benefit. It
directly addresses why the reduced-form γ is so noisy.

### 2.3 Separate fits vs one big hierarchical model — **my advice**

Three options:

- **A. Separate fits (status quo).** Simple, self-contained, easy
  π_α-across-samples comparison. But shares no strength, can't formally
  test dataset differences, and re-estimates δ_j four times.
- **B. Full joint hierarchical model.** Datasets as a grouping factor;
  partial pooling on σ_α, κ, (and γ); shared δ_j. Principled, lets us
  ask "does γ differ across datasets?" directly. But the channels are
  heterogeneous (CDI-only vs +input vs +input+LWL), so one monolithic
  Stan file is a big, brittle build.
- **C. Middle path (recommended).** Keep *per-dataset* fits for the
  per-child structure, but **share δ_j** via the two-stage anchor
  (§2.2). Compare γ across datasets descriptively first (like the
  5-sample π_α table). Only graduate to the full joint model (B) if the
  descriptive γ's look like they differ and we want a formal test.

**Recommendation: C now, B later.** The δ_j sharing is the high-value,
low-risk move; the full hierarchy is the publication-grade follow-up.

### 2.4 Sequencing: IO first, processing second

- **IO analysis (all 4 datasets):** CDI + input channel only. Uniform
  model across BabyView / SEEDLingS / AM2018 / FMW2013. This is where γ
  gets estimated.
- **Processing analysis (AM2018 + FMW2013 only):** add the LWL/RT
  channel as a *separate* second study, treating AM2018 and FMW2013 as
  distinct datasets (not pooled with the no-LWL ones). This is the
  γ_rt coupling work from experiments §9, redone on cleaned, WG+WS,
  wide-age data.

Do not entangle the two: get the input→slope question settled on clean
IO fits before layering processing back in.

---

## Part 3 — The γ analyses (the payoff, after Parts 1–2)

### 3.1 Reduced-form robustness (quick)

Re-run `input_slope_check.R` with the **raw observed per-child token
rate** (mean `log_r_obs` from the bundle — pure data, not the latent
`log_r_true`) as the predictor of ζ_i. If BabyView's γ survives with the
data-side predictor, the signal isn't latent-coupling. Can be done on
BabyView/SEEDLingS as soon as their bundles are confirmed WG+WS clean;
full version after all four datasets are re-prepped.

### 3.2 Model-based γ (the real test)

Add an input-on-slope term to the IO Stan model:

  η_ij(t) = log r_i + log α_i + [κ_pop + γ·(log r_i − μ_r) + ζ_i]·log(t/a₀) − δ_j

- `log r_i`: measured (anchored by observed tokens)
- γ: estimated input-on-slope coefficient (the new parameter)
- ζ_i: residual per-child slope variation after accounting for input
- γ = 0 → input is purely a level effect (clean accumulator)
- γ > 0 → bootstrapping: more input → steeper power law

Identifiability: γ is identified in the IO datasets because `log r_i` is
a measured covariate (between-kid input variation × within-kid age
variation; not collinear). Power is the constraint (N = 20–60), which is
exactly why the shared-δ_j anchor (§2.2) matters. Fit per dataset first,
compare γ̂; consider the joint model (§2.3B) for a pooled estimate.

### 3.3 Reconcile with "minimal input effects"

The slide-deck input-quartile plot is a *level* story (quartile
trajectories nearly overlap). If γ > 0 they should *fan out with age*.
Re-make that plot from the cleaned fits and check whether the high-input
quartile pulls away at older ages — a direct visual of γ.

---

## Open questions / risks

- **Naming crosswalk**: which Peekbank LWL dataset (`fernald_marchman_2012`,
  `fernald_totlot`, `fmw_2013`) maps to which CDI/LENA study? `fmw_2013`
  is presumably FMW2013/TLO; `fernald_totlot` (TOTLOT) may be the same
  or a precursor. Resolve before splitting.
- **β_react × γ confound** (§2.1) — sensitivity check needed.
- **Head-cam spatial sampling** (BabyView) is a different bias from
  daylong-LENA reactivity; the single multiplier may conflate them.
- **Small N**: even cleaned, BabyView N=20 limits γ precision. The
  shared-δ_j anchor is the main lever; pooling (joint model) is the
  fallback.
- **Does the reduced-form ζ-vs-log_r correlation survive** using the
  raw observed rate? If not, the model-based γ is unlikely to be real.

## Data inventory (audited 2026-05-28)

What we actually have on disk, per dataset × channel:

| Dataset | CDI (forms, ages, source) | Input (LENA) | LWL (Peekbank) | IO-ready? |
|---|---|---|---|---|
| **BabyView** | WG+WS, 9–32 mo, 20 kids (`babyview_subset_data.rds`) | head-cam AWC ✓ | — | **yes** |
| **SEEDLingS** | WG only, 6–18 mo, 44 kids (`seedlings/cdi_items_long.csv`) | `lena_data.csv` 560 rec / 44 kids ✓ | — | **yes** (WG-only is age-appropriate) |
| **AM2018 (TL3 / totlot3)** | **WS only**, 14–32 mo, 65 kids (`stanford_cdi_items_long.csv`) — **WG NOT in our files** | `TL3TLO_LENA.csv` TL3: 66 kids, AWC @ 16+18 mo ✓ | `adams_marchman_2018` 711 | **partial** — missing WG |
| **FMW2013 (TLO)** | **MISSING entirely** (no TLO/FMW rows in `stanford_cdi_items_long.csv`) | `TL3TLO_LENA.csv` TLO: 51 kids, AWC @ **18 mo only** | `fmw_2013` 178 | **blocked** — no CDI |
| **TotLot2 (TL2 / totlot2)** | WG (15–19, 97 kids) + WS, 409 admins ✓ | **not in `TL3TLO_LENA.csv`** | `fernald_totlot`? `fernald_marchman_2012`? | input unclear |

The parsed `stanford_cdi_items_long.csv` contains only `study ∈ {totlot2,
totlot3}` — **414 WS admins (181 kids) + 97 WG admins (97 kids), all
totlot2/3**. There is no FMW2013/TLO CDI in it, and no totlot3 (AM2018)
WG.

### Gaps that block / limit the IO analysis

1. **AM2018 (TL3) WG CDI is missing.** Only `TL3_compiled_WS.csv` exists;
   the WG-equivalent (TL3 ~16–18 mo) was never compiled into the long
   file. Without it, AM2018 stays in the WS-only narrow window we're
   trying to escape. *Need: the TL3 WG compiled export from the
   Marchman lab.*
2. **FMW2013 (TLO) CDI is missing entirely.** We have its LENA input (51
   kids @ 18 mo) and its LWL (`fmw_2013`, 178), but no item-level CDI in
   our files, so we can't fit it at all. *Need: TLO/FMW2013 CDI item
   data.*
3. **TotLot2 (TL2) has WG+WS CDI but no LENA in our files.** It could be
   a usable third IO cohort if its input exists somewhere; the
   `TL3TLO_LENA.csv` only covers TL3 + TLO. *Question: is TL2 one of
   AM2018/FMW2013 under a different label, or a separate cohort? Does it
   have LENA?*
4. **TLO LENA is single-timepoint (18 mo only); TL3 is 16+18 mo.** Fine
   for the IO model (input rate is treated age-invariant), but worth
   noting.

### Crosswalk still to confirm (needs lab knowledge)

| TotLot label | LENA `Study` | Peekbank LWL `dataset_name` | public name |
|---|---|---|---|
| totlot3 / TL3 | TL3 | `adams_marchman_2018` (711) | **AM2018** ✓ confirmed |
| ? | TLO | `fmw_2013` (178) | **FMW2013** ✓ confirmed |
| totlot2 / TL2 | (none) | `fernald_totlot` (229)? `fernald_marchman_2012` (679)? | ? |

The lab subject IDs (`11xxx`, `20xxx`) are the join key; `data/peekbank/README.md`
says TL2/TL3 IDs line up with `peekbankr::get_subjects()` for
`adams_marchman_2018` and `fmw_2013`.

### Bottom line — UPDATED 2026-05-28 (gaps resolved)

Mike added the missing CDI files. All four IO datasets now have the CDI
they need:
- **BabyView + SEEDLingS**: ready (β_react removal aside).
- **AM2018**: now WG+WS (TL3_compiled_WG.xlsx added) → full age range.
- **FMW2013**: now WG@18 + WS@24,30 (TLO_*.xlsx added).
- **FM2012**: WG+WS present but no LENA → reserved for the *processing*
  study, not IO.

Remaining precondition is engineering, not sourcing: extend
`parse_stanford_cdi.R` to ingest the new TL3-WG / TLO files (age from
filename, drop `misc` rows) into `stanford_cdi_items_long.csv`, then
build IO bundles with WG+WS for AM2018 and FMW2013.

## Concrete order of operations

1. Crosswalk Peekbank LWL names → studies; document in `data/README`.
2. Re-prep all four IO bundles with WG+WS (production), channels flagged.
3. Two-stage δ_j: extract from big EN longitudinal, wire into IO fits.
4. Re-fit IO (CDI+input) per dataset with shared δ_j; check β_react
   sensitivity.
5. γ analyses: reduced-form raw-rate (§3.1), then model-based γ (§3.2),
   then the fan-out plot (§3.3).
6. Processing analysis (AM2018 + FMW2013, +LWL) as a separate study.
