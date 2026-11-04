# TODO — refit δ_j on the current 681-item bundle (psi alignment done, refit open)

Opened 2026-11-04 (from the CLASSIC study, `studies/classic_chunking/`). Alignment fix on branch
`fix/psi-681-alignment`; the refit itself is **not done**.

## What was wrong
`fits/summaries/long_no_freq_slopes_psi.csv` (682 rows, `jj, delta_j_mean, delta_j_median`, no item column)
was fit on the **June-2026 bundle** (`446820b`, 682 items; psi committed in `4d514ed`, 06-10). The
**2026-08-15** rebuild (`c97256d`, "port the acceleration repo's clean extraction; wordbankr 2.0") renumbered
`jj`, renamed items (`a` → `id:a`, `in` → `ul:in`) and merged the old `in` and `inside/in` items into one
(`ul:in`), giving 681. Joined by `jj` to the current bundle, psi has cor ≈ 0 with production difficulty
(0.007); on the June bundle it is 0.990 (checked against every historical version of `long_subset_data.rds`).
The file had no item column, so no consumer could tell.

## Were shipped results affected? No.
The anchored IO bundles were built while psi and bundle still matched, and carry δ_j anchors that equal the
name-aligned psi exactly:

| Bundle | Anchored to psi vintage | Check |
|---|---|---|
| `io_fmw2013_subset_data`, `io_pooled_subset_data` (Fig 3) | 2026-06-10 (`4d514ed`) | re-derived with the new loader: max \|diff\| = 0.00000 on all shared anchored items (670 FMW2013, 681 pooled) |
| `io_am2018_subset_data`, `seedlings_subset_data` | **2026-05-23** (`86e7856`, 671 items) | exact match to the May psi mapped through its bundle (`018c0fb`) |

So the pooled/FMW2013 fits and the am2018/seedlings fits rest on **different δ_j vintages** (shifted by
~0.4–0.8 in anchor value); each is self-consistent. Whether to re-anchor am2018/seedlings is a modelling
decision, not made here.

The exposure was **re-running** the anchoring scripts: `add_io_anchor.R`, `prepare_io_pooled.R` and
`prepare_io_dataset.R` joined psi to the bundle by `jj` with no length check, then keyed on `item`
strings that no longer match the current bundle's names. `paper/build_cache.R` (Fig 4) already had a guard
and stops loudly. The local `fits/summaries/long_no_freq_slopes.draws.rds` that build_cache also needs is
not in the checkout, so Fig 4 was not rebuilt here.

## Done on `fix/psi-681-alignment`
- `model/scripts/export_psi_aligned.R`: re-keys the June fit to the current bundle by `item_definition`
  (recovers the June bundle from git); writes `long_no_freq_slopes_psi.csv` with `jj, item, item_definition,
  delta_j_mean, delta_j_median, source_jj` and `long_no_freq_slopes_psi.README.txt`. 681 / 681 items matched,
  cor with production difficulty (22–28 mo) 0.984, cor 0.969 with the independent 197-item named fit.
- `model/R/psi_anchor.R`: `load_en_psi_by_name()` (requires an `item` column and an item-set match with
  the bundle, else stops) and `check_anchor_coverage()` (stops on 0 anchored). The three consumers use them.
- Not changed: shipped bundles in `fits/` and everything fitted from them.

## Still open
1. **Refit** the EN longitudinal model on the current 4.25M-row bundle (cluster; GCP/Sherlock scripts in
   `cluster/`) and export psi with `item`; the aligned file carries June-fit values on the old data.
   Expect δ_j to move; compare to the aligned file first (should still correlate ~0.98 with it).
2. Decide whether to re-anchor `io_am2018` / `seedlings` (May vintage) to the current psi, and re-fit.
3. Rebuild Fig 4 (`paper/build_cache.R`) once the draws file is available; the name join now returns 681 rows.
4. Drop the stale `jj`-keyed copies in `fits/summaries/long_no_freq_slopes_{english_I200,norwegian}_psi.csv`
   from consumers if they are joined positionally anywhere (not audited here).

## Why it matters for CLASSIC
`studies/classic_chunking/` can now test on all 681 items with a δ_j from v2's fit instead of the 197-item
named fit plus empirical difficulty.
