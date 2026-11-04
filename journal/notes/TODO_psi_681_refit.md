# TODO — regenerate the full 681-item δ_j file (`long_no_freq_slopes_psi.csv`)

Opened 2026-11-04 (from the CLASSIC study, `studies/classic_chunking/`). Not started.

## Problem
`fits/summaries/long_no_freq_slopes_psi.csv` (Oct 26; 682 rows, columns `jj, delta_j_mean, delta_j_median`)
does not line up with the current bundle `fits/long_subset_data.rds` (681 items in `word_info`, 4.25M-row `df`).
Measured on 2026-11-03: joined by `jj`, its δ_j has correlation ≈ −0.01 with item production rate at
22–28 mo (the 197-item named fit `figs/longitudinal/psi_freq_regression_per_word.csv` has −0.97), and
is uncorrelated with that named fit's ψ for the same items. The file carries no item column, so any
consumer can only join positionally. The same failure is described in `journal/HANDOFF.md` ("Stale psi")
and in the comment above the join in `paper/build_cache.R`.

## What to do
1. Refit (or re-extract from the existing draws, if they match the 681-item bundle) δ_j for all 681 items
   and write the summary **with an `item` column** (and `jj`), so joins can be by name.
2. Assert `nrow(psi) == nrow(word_info)` and name-match wherever the file is read.
3. Sanity check before use: cor(δ_j, −logit(production rate, 22–28 mo)) should be ≈ 0.95+.

## Consumers to audit (all read this file)
- `paper/build_cache.R` (Fig 4) — has a name-join guard; confirm it still fires correctly.
- `model/scripts/prepare_io_pooled.R` and `model/scripts/add_io_anchor.R` — join to `word_info` by **`jj`**
  with no length check; these set the global δ_j anchor for the pooled input-observed fits (Fig 3).
  If the file is misaligned they anchor δ_j to the wrong words. Check whether the shipped fits were built
  from an aligned copy.
- `model/scripts/prepare_io_dataset.R`, `journal/experiments.md` (mentions only).

## Why it matters for CLASSIC
The CLASSIC comparison currently uses the 197-item named fit (195 items heard in the streams) plus empirical
difficulty. With a valid 681-item δ_j the test uses all items and v2's own measure instead of the proxy.
Re-run: update `00_export_v2_items.R` to read the new file (join by name) and `sh run_all.sh`.
