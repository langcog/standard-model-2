# Peekbank / Stanford-Marchman raw data

Source data for the io-proc (input + processing) model: item-level CDI, LWL
reaction time, and daylong LENA input for the four Stanford-Marchman LWL cohorts,
plus the peekbank-development LWL summaries they were drawn from.

## ⭐ Label correspondence (read this first)

The lab's internal cohort labels (`TL`, `TL2`, `TL3`, `TLO`) are opaque and were a
recurring source of confusion. The authoritative key — verified against the data
2026-06-14 (channel counts below are kids with that channel):

| lab label | model `dataset_name` | publication | CDI | RT (LWL) | **input (LENA)** |
|---|---|---|:--:|:--:|:--:|
| `TL`  / "TotLot" (orig.) | **`fernald_totlot`** | Fernald, Perfors & Marchman 2006 | ✓ | ✓ 63 | **✗ none** |
| `TL2` / "TotLot 2"       | **`fernald_marchman_2012`** | Fernald & Marchman 2012 | ✓ | ✓ 122 | **✗ none** |
| `TL3` / "TotLot 3"       | **`adams_marchman_2018`** | Adams, Marchman et al. 2018 | ✓ | ✓ 69 | **✓ 66** |
| `TLO` / "TotLot Outreach"| **`fmw_2013`** | Fernald, Marchman & Weisleder 2013 | ✓ | ✓ 80 | **✓ 51** |

**Notes / caveats (don't skip):**
- **Only TL3 and TLO have daylong LENA input.** TL2 (`fernald_marchman_2012`) and
  TL (`fernald_totlot`) are **CDI + RT only — no input** anywhere in the repo. If
  `fernald_marchman_2012` ever collected daylong audio, it is **missing** and must
  be sourced; do not assume it has an input channel.
- The single file `lena_am2018_fmw2013.csv` (formerly `TL3TLO_LENA.csv`) carries
  **both** TL3 and TLO input, split by its internal `Study` column (`TL3`/`TLO`).
- **TLO input is under-used.** All 51 TLO LENA kids load, but the proc/joint bundle
  gates input on RT+CDI, so only the ~31 TLO kids who also have RT enter. The ~20
  TLO kids with input+CDI but **no RT** could contribute as input-only kids (as
  BabyView/SEEDLingS do) but currently don't — an open salvage lever, not a bug.

## Directory layout

```
data/peekbank/
  adams_marchman_2018/      TL3_compiled_WS.csv, TL3_compiled_WG.xlsx     (raw CDI)
  fernald_marchman_2012/    TL2_WG_compiled.xlsx, TL2_WS_compiled.xlsx    (raw CDI)
  fmw_2013/                 TLO_18m_WG.xlsx, TLO_24_WS.xlsx, TLO_30m_WS.xlsx (raw CDI)
  fernald_totlot/           TL_18m_WS.xlsx, TL_21m_WS.xlsx, TL_25m_WS.xlsx (raw CDI)
  lena_am2018_fmw2013.csv   daylong LENA (AWC/CTC/CVC @16+18mo) for TL3 + TLO
  1_d_sub.Rds               per-admin LWL summary (RT + accuracy) — peekbank 2026.1
  0_cdi_subjects.Rds        admin-level CDI scores (peekbank-development)
  stanford_cdi_items_long.csv   long CDI, TL2/TL3/TLO   (built by parse_stanford_cdi.R)
  totlot_cdi_items_long.csv     long CDI, TL             (built by parse_totlot_cdi.R)
  cdi_short_code_map_{ws,wg}.csv  short-code -> Wordbank item_definition
  cdi_master_item_key.csv   human-checkable item merge key (build_cdi_master_key.R)
```

## Which script reads what

- `parse_stanford_cdi.R` → reads `adams_marchman_2018/`, `fernald_marchman_2012/`,
  `fmw_2013/` raw CDI → writes `stanford_cdi_items_long.csv` (+ short-code maps).
- `parse_totlot_cdi.R` → reads `fernald_totlot/` raw CDI → `totlot_cdi_items_long.csv`.
- `prepare_proc_dp_bundle.R` → reads `lena_am2018_fmw2013.csv` for the input channel
  (maps internal `Study` TL3→`adams_marchman_2018`, TLO→`fmw_2013`); RT from `1_d_sub.Rds`.
- Downstream bundle builders read the long intermediate CSVs, **not** the raw files,
  so the per-study raw dirs only matter when re-parsing.

## Provenance

- **LWL summaries** (`1_d_sub.Rds`, `0_cdi_subjects.Rds`): exact copies of
  `cached_intermediates/` from <https://github.com/peekbank/peekbank-development>
  (Frank et al. 2026, *eLife*), peekbank DB version `2026.1`. Refresh by re-running
  that repo's `0_get_data.qmd` + `1_tidy_data.qmd`.
- **Raw CDI** (`TL*`/`TLO*` files): provided directly by the Stanford Marchman lab
  (wide-format exports from paper forms). Lab subject IDs (`11xxx`, `20xxx`) line up
  with `peekbankr::get_subjects()` `lab_subject_id` for the matching datasets.
- **LENA** (`lena_am2018_fmw2013.csv`): Marchman-lab daylong recordings, AWC/CTC/CVC
  per hour at 16 and 18 mo (two replicate reads/child; 63 have 16M, 114 have 18M).

## What's NOT here
- Item-level Peekbank CDI responses (Peekbank stores only totals — hence the raw
  lab files above).
- Frame/trial-level LWL (`0_d_aoi.Rds`, `1_d_trial.Rds`) — too large; re-derive from
  peekbank-development if needed.
