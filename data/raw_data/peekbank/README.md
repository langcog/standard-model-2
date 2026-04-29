# Peekbank cached intermediates

These files are copied verbatim from the peekbank-development repo
(`cached_intermediates/`) and represent the per-admin and CDI-joined
data used in the Peekbank developmental-relationships paper.

## Files

- `1_d_sub.Rds` (~380 KB) — one row per (subject, administration).
  Built by `1_tidy_data.qmd` in peekbank-development. Columns:
  - `dataset_name`, `subject_id`, `administration_id`, `age`, `log_age`
  - Trial counts: `n_trials`, `n_trials_rt`
  - Accuracy: `short_window_accuracy`, `long_window_accuracy`,
    `*_window_acc_var`, `*_window_elogit`
  - Reaction time: `rt`, `log_rt`, `rt_var`, `log_rt_var`
    (only D-T shifts; RT < 367 ms set to NA)
  - CDI: `prod`, `comp` (rawscore / instrument_length, fuzzy-joined
    within ±1 month of the LWL admin age)
  - Coverage: `*_window_prop_data`

- `0_cdi_subjects.Rds` (~50 KB) — admin-level CDI scores from
  peekbankr's subject_aux_data, normalized to instrument_length.
  Already merged into `1_d_sub.Rds`; kept here for provenance and
  in case we want to re-do the fuzzy join with different parameters.

## Source data version

Pulled from peekbank database `2026.1`. To refresh:

1. Run `0_get_data.qmd` and `1_tidy_data.qmd` in
   `~/Projects/peekbank/peekbank-development/` against the desired
   `db_version`.
2. Copy the resulting `cached_intermediates/0_cdi_subjects.Rds` and
   `1_d_sub.Rds` into this folder.

We do not depend on peekbankr / wordbankr at runtime, so Sherlock
fits read this snapshot directly.

## What's NOT here

- Item-level CDI responses (Peekbank only stores totals).
- Raw frame-level (`0_d_aoi.Rds`) and trial-level (`1_d_trial.Rds`)
  data — too large to commit. Re-derive from Peekbank if needed.
