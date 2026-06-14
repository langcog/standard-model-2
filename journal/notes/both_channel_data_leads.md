# Both-channel (input + processing) data leads

The io-proc separation (input vs processing contribution to efficiency) is
identified only by kids with **input AND processing AND vocabulary** in the same
child. We currently have **97** (AM2018 66 + FMW2013 31). These are the leads for
growing that, from two background lit/data reviews (2026-06-13).

## Lead 1 — Mahr & Edwards (2018, Dev Sci 21:e12685) — OUT (too old)
- n=109, English (Wisconsin); LENA input + LWL eye-tracking + vocabulary, same kids.
- **Fully downloadable**: `github.com/tjmahr/mahr-edwards-2017`,
  `data/03_input_vocab_eyetracking.csv` = analysis-ready merge of all 3 measures,
  keyed by `ResearchID`; raw frame-level eye-tracking also present.
- **Why it's out (MCF):** ages **28–39 mo at T1** — past our 13–30 mo window.
- Other caveats if revisited: vocab is **PPVT/EVT standardized scores, not
  item-level CDI** (needs a sumscore-vocab bridge); data formally "belong to
  UW-Madison" (permission needed to publish). Tristan Mahr is responsive.

## Lead 2 — SEEDLings LWL (Bergelson) — PAPER IN HAND (2026-06-14)
**Zhu, Amatuni, Egan-Dailey, Garrison, Kalenkovich, Koorathota, Righter, Tor,
Bergelson — "Experience Shapes Early Noun Comprehension from 8–18 Months"**
(submitted; PDF `papers/zhu_etal_submitted.pdf`, 93pp w/ SI). Data+code+stimuli:
**osf.io/m2kdz** (API auth-gated while under review — need MCF's login / a
view-only link / a manual download into `data/`).
- **Study 1 (n=44) = our SEEDLings cohort** (already in the io-proc model as
  *input-only*: LENA + monthly WG-CDI). Hand-tailored LWL every 2 mo, **8–18 mo**
  (≤6 sessions × 16 words/session). They ALSO computed **per-child per-word
  frequency from each child's own home recordings** (the "audio-nouns" were
  hand-picked from that child's prior 2 mo of daylong audio).
- **Study 2 (n=247)** = cross-sectional controls on Study-1 kids' tailored
  stimuli; NOT our cohort, no home recordings → no input channel → not useful here.

### Verdict for io-proc
- **Friction 1 — measure:** infant outcome is *increase in target looking*
  (accuracy), NOT Fernald RT (our channel is log-RT). RT computed for **adults
  only** (Figs E1/E2). They explicitly punt RT to future work and invite reuse
  ("reaction time and target-looking accuracy are of course closely linked").
  → RT is *derivable* from the raw OSF timecourse, their adult-RT code a head
  start, but it is pipeline work (± lab coordination for raw frames).
- **Friction 2 — age:** SEEDLings LWL 8–18 mo vs our Fernald RT channel **13–30 mo**
  (AM2018 13–30, FMW2013 17–26, FM2012 18–30, totlot 15–25). Overlap only
  **13–18 mo** (14/16/18-mo sessions); 8–12 mo is below window + near floor.
- **Upside:** the 44 kids → both-channel would grow the separation bottleneck
  **97 → ~141 (+45%)** — but most new processing data is below window. Treat
  LWL-as-processing as a **stretch / robustness add, not a core dependency**.
- **The real prize = the per-child per-word frequency** — observed child-specific
  *word-level* input, far richer than our LENA total-rate, no measure mismatch,
  and exactly the quantity the model currently imputes. Worth grabbing regardless
  of the LWL question; feeds a future word-level input channel.

### UPDATE 2026-06-14 — RAW EYETRACKING IN HAND (`data/seedlings/raw_eyetracking_data/`)
Friction 1 is **largely resolved**: RT is derivable from data we now have, no lab
coordination needed.
- **HaT/** = the hand-tailored Study 1 (our cohort). EyeLink Data Viewer exports
  (tab-delim, `.xls`): `fix_rep_eighttoeighteenmonth` (163,797 fixation rows, 235
  cols) + `mes_rep_eighttoeighteenmonth` (trial timing). NOT in peekbank d_sub
  (it's kid-specific-words, un-standardizable) → must derive ourselves.
- **Session label = the crosswalk**: `RECORDING_SESSION_LABEL` is literally
  `subj_age` (`01_08`…`44_18`), subj 01–44 = our LENA/CDI ids. 264 sessions ≈
  44 kids × 6 ages (8/10/12/14/16/18 mo), ~27k fixations/age.
- **RT is computable, standard pipeline**: `mes_rep.IP_START_TIME` = target word
  onset; `fix_rep` has CURRENT_FIX_START + CURRENT_FIX_INTEREST_AREA_LABEL +
  TargetSide/TargetImage. Fernald RT = first target-AOI fixation start − onset,
  distractor-initial trials, 300–1800 ms window. (Same computation peekbank runs.)
- **In-window coverage**: ages **14/16/18 mo overlap our 13–30 RT window** → ~44
  kids × 3 sessions ≈ **132 kid-sessions in window** (8/10/12 below; could extend
  the window down). Moving these to both-channel: bottleneck **97 → ~141 (+45%)**.
- **DiSCo/** = `dsc_{fix,mes}rep_001_380.Rds` (405+170 MB, already R-processed) —
  almost certainly the Study-2 cross-sectional controls (not our cohort, no home
  recordings) → lower priority; heavy, don't load casually.
- **Revised verdict:** no longer a "stretch" — **feasible and high-value**, gated
  on (a) building+QC'ing the RT extractor and (b) the age overlap (half the
  sessions below window). De-risk with a 1–2-session proof-of-concept first.

### BUILT + QC'd 2026-06-14 — `model/scripts/prepare_seedlings_lwl_rt.R`
Extractor done; QC report card vs paper targets in
[`seedlings_lwl_qc_targets.md`](seedlings_lwl_qc_targets.md). Output
`data/seedlings/seedlings_lwl_rt.csv` (1,697 RTs, 44 kids, median 710 ms,
774→642 ms over 8→18 mo; 1,097 in our 14–18 mo window). Added to the RT spaghetti
plot. **Onset marker = `EL_BUTTON_CRIT_WORD` (experimenter button), not
`IP_START_TIME`** (which is just trial start) — corrects the line above.
**⚠ The "trivial `01`–`44`" crosswalk below is WRONG**: subjects are `01..46` with
gaps (05, 24), io_pooled uses a dense factor rank → naive `::N=0N` misaligns 40/44
kids. Bundle wiring is gated on a verified subject_id crosswalk (recommend
retaining `subject_id` through prepare_seedlings → io_pooled). See the QC note.

### Decision this forces for the measurement-model redesign (#16)
Let the processing channel accept **accuracy (target-looking) alongside log-RT**
(heterogeneous-processing measurement model) → SEEDLings plugs in WITHOUT deriving
RT (cost: a 2nd efficiency→accuracy link + the young-age floor). Same lever also
unlocks Mahr&Edwards-style data later.

- **Access paths:** processed looking + per-child freq → OSF m2kdz (auth-gated);
  ID crosswalk to our LENA/CDI is trivial (`subj` = `"01"`–`"44"`, drop `_06`
  suffix, drop `ns01–ns12` lab-only). raw audio = HomeBank pw; home video/raw
  EyeLink = Databrary authorized-investigator. Contact: elika.bergelson@harvard.edu.

## Takeaway for the model
Both leads have a **format mismatch** with the current model: M&E vocab is
PPVT/EVT (sumscore), SEEDLings processing is accuracy (not RT). Growing the
both-channel set therefore likely needs the model to accept **heterogeneous
vocab** (item-CDI + sumscore) and/or **heterogeneous processing** (RT + accuracy).
Worth a design conversation alongside the measurement-model spec — it would
unlock both at once. Dead ends (single-channel, documented in the lit-review
results): Casillas Tseltal, Pomper/Saffran, Newman/Ratner (segmentation not RT),
TALK study (MEG, not yet released), Houston/Wang.
