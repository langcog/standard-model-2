# SEEDLingS raw data

Two source streams:
1. **Per-recording LENA + per-child CDI totals** from the Egan-Dailey
   & Bergelson (2025) paper.
2. **Item-level CDI** from a separate Bergelson lab repository
   (Dong & Bergelson 2026 working materials).

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
- **Caveat (still being verified):** the `subj` column appears to
  use sequential integer IDs that should match the "01".."44"
  convention in `seedlings_data.csv`. Mike is confirming this
  alignment before we trust the linkage.

## Files

| file | what it contains |
|---|---|
| `lena_data.csv` | per-recording LENA stats (560 rows, 44 kids × ~13 monthly audio recordings 0;6-1;5). Columns: `subj`, `month`, `duration_hrs`, `awc_perhr`, `cvc_perhr`, `ctc_perhr`, outlier flags |
| `seedlings_data.csv` | per-child summary used in Egan-Dailey 2025: aggregated input measures + CDI totals + later assessments |
| `all_vocab.csv` | per-child CDI totals at 8/12/18 months + later language scores |
| `cdi_ht_raw_temp.csv` | wide-format CDI item-level data, longitudinal (Dong & Bergelson 2026 working file) |

## What's NOT here

- Manual noun-token annotations (the high-fidelity input measure from
  Bergelson et al.). Available in
  <https://github.com/BergelsonLab/seedlings-nouns> (~46 MB raw)
  if/when we want to add a manual-input sensitivity check.
- Raw audio / video recordings — Databrary-restricted
  (<https://databrary.org/party/61>).

## Pipeline

`model/scripts/prepare_seedlings.R` consumes `lena_data.csv` for
per-recording `log_r_obs` and (once subject-ID alignment is verified)
will pivot `cdi_ht_raw_temp.csv` into the long item-level format
`cdi_items_long.csv` expected by the io-model bundle builder.
