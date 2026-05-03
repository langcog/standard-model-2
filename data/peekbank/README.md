# Peekbank raw data

Two distinct sets of inputs live here:

1. **Per-admin LWL summaries** lifted from the peekbank-development
   working repo for the developmental-relationships paper.
2. **Stanford TotLot CDI source files** (Marchman lab) for use as
   item-level CDI inputs to our longitudinal IRT pipeline. These are
   the same children whose looking-while-listening data appears in
   `adams_marchman_2018`, `fmw_2013`, `fernald_marchman_2012`, and
   `fernald_totlot` Peekbank datasets.

## Provenance

### Per-admin LWL summaries
- `1_d_sub.Rds`, `0_cdi_subjects.Rds`
- Source: <https://github.com/peekbank/peekbank-development>
  (working repo for Frank et al., 2026, *eLife*).
- Database version: peekbank `2026.1`.
- These are exact copies of `cached_intermediates/1_d_sub.Rds` and
  `cached_intermediates/0_cdi_subjects.Rds` from that repo. We do
  not modify the source repo; refresh by re-running the
  peekbank-development `0_get_data.qmd` and `1_tidy_data.qmd`
  notebooks and copying the resulting Rds files here.

### Stanford TotLot CDI source files
- `TL2_WG_compiled.xlsx`, `TL2_WS_compiled.xlsx`, `TL3_compiled_WS.csv`
- Source: provided directly by the Stanford Marchman lab (Mike Frank
  has copies). These are wide-format CDI exports compiled from
  paper forms collected over the lifetime of the TotLot 2 and
  TotLot 3 longitudinal cohorts.
- Lab subject IDs (`11xxx`, `20xxx`) line up with the
  `lab_subject_id` field exposed by `peekbankr::get_subjects()` for
  `adams_marchman_2018` and `fmw_2013`. The TotLot 2 study covers
  WG admins at ~16 mo and WS admins from ~18-30 mo.

## Files

| file | rows | what it contains |
|---|---|---|
| `1_d_sub.Rds` | 4 124 admins × 25 cols | per-(subject, admin) LWL summary: RT, accuracy, CDI totals (`prod`, `comp` as proportion of form length) |
| `0_cdi_subjects.Rds` | one row per CDI admin in the LWL set | admin-level CDI scores normalized to instrument length |
| `TL3_compiled_WS.csv` | 185 rows × 821 cols | TotLot 3 (Adams 2018 cohort): WS form, ages 20-32 mo, 65 kids |
| `TL2_WS_compiled.xlsx` | 347 rows × 822 cols | TotLot 2: WS form, ages 18-30 mo, 119 kids |
| `TL2_WG_compiled.xlsx` | 120 rows × 522 cols | TotLot 2: WG form, ages 15-19 mo, 97 kids |

## Derived files (built by `model/scripts/parse_stanford_cdi.R`)

- `cdi_short_code_map_ws.csv`, `cdi_short_code_map_wg.csv` — short-code
  → Wordbank `item_definition` mapping. All 1 076 entries (680 WS,
  396 WG) resolve via `auto_exact` (fingerprint match) or
  `manual_disambig` (deterministic Marchman-lab conventions). Hand
  review of the manual_disambig rows is logged as a loose end in
  `outputs/experiments.md`.
- `stanford_cdi_items_long.csv` — long format, one row per
  (lab_subject_id, age, form, item, produces). 183 subjects × 511
  admins × 681 items ≈ 321K rows. Placeholder rows (no parent
  return) are filtered out.

## What's NOT here

- Item-level Peekbank CDI responses — Peekbank only stores totals.
- Raw frame-level (`0_d_aoi.Rds`) and trial-level (`1_d_trial.Rds`)
  data from peekbank-development — too large to commit. Re-derive
  from Peekbank if needed.
