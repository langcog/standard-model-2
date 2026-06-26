# SEEDLingS raw data

Three source streams:
1. **Per-recording LENA + per-child CDI totals** from the Egan-Dailey
   & Bergelson (2025) paper.
2. **Item-level CDI** from a separate Bergelson lab repository
   (Dong & Bergelson 2026 working materials).
3. **Looking-while-listening (LWL) eyetracking** from Zhu et al. (2026) —
   the source of the processing (reaction-time) channel.

## Provenance

### Per-recording LENA + per-child CDI totals
- `lena_data.csv`, `seedlings_data.csv`, `all_vocab.csv`
- Source: <https://github.com/ShannonDailey/EganDailey_DevPsy_2025>
  (specifically `data/`).
- Publication: Egan-Dailey, S. & Bergelson, E. (2025). Early child
  measures outpredict input measures of preschool language skills in
  U.S. English learners. *Developmental Psychology* (advance online
  publication). DOI: `10.1037/dev0002019`.
- The 44 SEEDLingS children were enrolled 2014-2015 in Rochester, NY
  and recorded monthly from 0;6 to 1;5 plus a 4;6 follow-up.

### Item-level CDI (`cdi_ht_raw_temp.csv`)
- Source: <https://github.com/BergelsonLab/WordExposure/blob/main/data/ht/cdi_ht_raw_temp.csv>
- Publication: from Dong & Bergelson (2026) working materials.
- Wide-format export with `Talk_<item>` (production) and
  `Understand_<item>` (comprehension) columns plus `subj`, `month`,
  `CDIcomp`, `CDIprod`, `Date_Completed`, and a
  `SeedlingsFinalSample` flag. Longitudinal: rows for every monthly
  admin from ~6 to ~18 mo per child.
- `model/scripts/parse_seedlings_cdi.R` pivots this into the long
  item-level `cdi_items_long.csv` (mapping the integer `subj` to the
  `01`..`44` convention used in `seedlings_data.csv`).

### Looking-while-listening (LWL) eyetracking — the processing channel
- `raw_eyetracking_data/HaT/{fix_rep,mes_rep}_eighttoeighteenmonth_7-26-2016.xls`
  — hand-tailored ("HaT") EyeLink fixation + message reports, 8-18 mo.
- Publication: **Zhu, L. Z., Amatuni, A., Egan-Dailey, S., Garrison, H.,
  Kalenkovich, E., Koorathota, S., Righter, L., Tor, S., & Bergelson, E.
  (2026). Experience Shapes Early Noun Comprehension from 8-18 Months: The
  Roles of Word Frequency and Referent Familiarity.** PsyArXiv preprint.
  <https://osf.io/preprints/psyarxiv/zchbj_v2>
- Study 1 (n=44, the SEEDLingS cohort) ran hand-tailored LWL eyetracking
  every two months from 8-18 mo. `model/scripts/prepare_seedlings_lwl_rt.R`
  derives a per-session Fernald-style reaction time from these reports
  (QC'd against the paper's trial-funnel + accuracy targets, see
  `journal/notes/seedlings_lwl_qc_targets.md`) → `seedlings_lwl_rt.csv`
  (`dataset_name = seedlings_zhu`), which enters the io-proc model as the
  child-level RT (processing) channel.
- `raw_eyetracking_data/DiSCo/` holds reports from a second eyetracking
  paradigm in the same study; not currently used.

## Files

| file | what it contains |
|---|---|
| `lena_data.csv` | per-recording LENA stats (560 rows, 44 kids × ~13 monthly audio recordings 0;6-1;5). Columns: `subj`, `month`, `duration_hrs`, `awc_perhr`, `cvc_perhr`, `ctc_perhr`, outlier flags |
| `seedlings_data.csv` | per-child summary used in Egan-Dailey 2025: aggregated input measures + CDI totals + later assessments |
| `all_vocab.csv` | per-child CDI totals at 8/12/18 months + later language scores |
| `cdi_ht_raw_temp.csv` | wide-format CDI item-level data, longitudinal (Dong & Bergelson 2026 working file) |
| `cdi_items_long.csv` | long item-level CDI — output of `parse_seedlings_cdi.R` |
| `raw_eyetracking_data/HaT/*.xls` | hand-tailored LWL EyeLink fixation + message reports (Zhu et al. 2026, Study 1) |
| `raw_eyetracking_data/DiSCo/*.Rds` | second-paradigm eyetracking reports (not used) |
| `seedlings_lwl_rt.csv` | per-session Fernald RT derived from the HaT reports — output of `prepare_seedlings_lwl_rt.R` |

## What's NOT here

- Manual noun-token annotations (the high-fidelity input measure from
  Bergelson et al.). Available in
  <https://github.com/BergelsonLab/seedlings-nouns> (~46 MB raw)
  if/when we want to add a manual-input sensitivity check.
- Raw audio / video recordings — Databrary-restricted
  (<https://databrary.org/party/61>).

## Pipeline

- `model/scripts/parse_seedlings_cdi.R` → `cdi_items_long.csv`
  (+ `cdi_seedlings_short_code_map.csv` audit).
- `model/scripts/prepare_seedlings_lwl_rt.R` → `seedlings_lwl_rt.csv`
  (the RT/processing channel).
- `model/scripts/prepare_seedlings.R` consumes `lena_data.csv` for
  per-recording `log_r_obs` and `cdi_items_long.csv` for the item-level
  vocabulary, building the SEEDLingS bundle for the io-model.
