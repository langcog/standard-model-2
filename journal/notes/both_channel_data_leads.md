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

## Lead 2 — SEEDLings LWL (Bergelson) — PROMISING, MCF following up
- The SEEDLings home cohort (our n=44, LENA + monthly CDI) WAS run in
  looking-while-listening: **Bergelson & Aslin 2017 (PNAS 114:12916)**, repo
  `ebergelson/sixmonth_seedlings_paper` (companion: Bergelson & Swingley 2018,
  Child Dev 89:1567).
- **IDs crosswalk one-to-one**: the eye-tracking `subj` field is literally
  `"01"`–`"44"` — the same SEEDLings child-id scheme as our LENA/CDI data. Linkage
  is trivial (drop the `_06` month suffix; `ns01–ns12` are lab-only kids to drop).
- **The catch:** the published processing measure is **proportion-target-looking
  (accuracy), not Fernald RT/latency** — our model's processing channel is log-RT.
  Public data is mostly **6 mo** (too young for our window). Raw EyeLink output
  reportedly carries a latency column → RT *could* be derived at 12/18 mo, but
  that needs the lab.
- **Access:** processed looking data open on GitHub; raw audio = HomeBank
  password; home video/lab = Databrary authorized-investigator. Contact:
  elika.bergelson@harvard.edu. **MCF is following up with Bergelson.**

## Takeaway for the model
Both leads have a **format mismatch** with the current model: M&E vocab is
PPVT/EVT (sumscore), SEEDLings processing is accuracy (not RT). Growing the
both-channel set therefore likely needs the model to accept **heterogeneous
vocab** (item-CDI + sumscore) and/or **heterogeneous processing** (RT + accuracy).
Worth a design conversation alongside the measurement-model spec — it would
unlock both at once. Dead ends (single-channel, documented in the lit-review
results): Casillas Tseltal, Pomper/Saffran, Newman/Ratner (segmentation not RT),
TALK study (MEG, not yet released), Houston/Wang.
