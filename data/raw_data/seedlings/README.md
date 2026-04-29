# SEEDLingS data staged for the io-model pipeline

## Public files (already in this folder)

- `lena_data.csv` — per-recording LENA stats (560 rows, 44 kids × ~13 monthly
  audio recordings 0;6–1;5). Columns: `subj`, `month`, `duration_hrs`,
  `awc_perhr`, `cvc_perhr`, `ctc_perhr`, outlier flags. Source:
  https://github.com/ShannonDailey/EganDailey_DevPsy_2025
- `seedlings_data.csv` — per-child summary used in Egan-Dailey & Bergelson
  (2025): aggregated input measures + CDI totals + later assessments.
  Same source.
- `all_vocab.csv` — per-child CDI totals at 8/12/18 months + later
  language scores. Same source.

These are sufficient for a between-child correlational analysis but
**not** for fitting the input-observed model — that needs item-level
CDI responses, which the public data do not include.

## Missing file: `cdi_items_long.csv` (item-level CDI)

`prepare_seedlings.R` will refuse to build a bundle until this file
exists. It must be in long format with exactly these columns:

| column     | type    | description                                     |
|------------|---------|-------------------------------------------------|
| subject_id | char    | "01".."44", matching `subj` in `lena_data.csv`  |
| age        | numeric | age in months at the CDI admin (~12 or ~18)     |
| form       | char    | "WG" (only WG is collected at these ages)       |
| item       | char    | CDI item-definition string (matches Wordbank)   |
| produces   | 0/1     | 1 if parent endorsed "produces", else 0         |

One row per (subject, age, form, item) cell. Expected size:
44 children × 2 ages × ~396 WG items ≈ 35K rows.

## How to obtain it

The Bergelson lab maintains a private `cdi_spreadsheet` repository
(loaded by `blabr::get_cdi_spreadsheet()`). A one-shot export of
SEEDLingS WG admins at 1;0 and 1;6 in long item-level format is what
we need — equivalent to:

```r
library(blabr)
cdi <- get_cdi_spreadsheet(version = '<latest>', type = 'csv')
# ... or, from a raw WebCDI export csv:
seedlings_wg <- wrangle_web_cdi(filepath = "<seedlings.csv>",
                                form = "WG", table = "wordlevel")
```

The `wordlevel` table from `wrangle_web_cdi()` has columns
`study_name`, `subject_id`, `repeat_num`, `item`, `response`. We just
need rows where `study_name` is the SEEDLingS study and
`response %in% c("produces", "understands", "neither")`, recoded to
the 0/1 `produces` column above.
