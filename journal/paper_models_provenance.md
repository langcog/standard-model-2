# Paper model-fit provenance

**Purpose.** One place to answer "are the fits in the manuscript current?" For every
model fit whose output reaches `paper/standard_model.qmd` (via `paper/cache/*.rds`),
this records: what the model is, where the fit code lives, which cluster trained it,
when, and where its outputs are stored.

**Last audited:** 2026-06-13 (dates below are file mtimes + the cited
`journal/experiments.md` entries; treat them as "as of this audit").

**Path note (2026-08-15):** the manuscript is now `paper/input_paper.qmd`
(successor to `paper/standard_model.qmd`, later `paper/input_paper/input_paper_text.qmd`);
retired papers moved to `paper/old/`. The cache chain below is unchanged.

**The dependency chain.** manuscript chunk → `paper/cache/<x>.rds` → a build script
(`paper/build_cache.R` or `paper/build_input_cache.R`) → one or more **fit outputs**
in `fits/` → trained by a **fit driver** on a **cluster**. A fit is "stale" when a
fit output is older than its inputs, or when a build script points at a fit that has
been superseded.

**To refresh after new fits land:** re-run `Rscript paper/build_cache.R` and
`Rscript paper/build_input_cache.R`, then re-render. Both read committed fit outputs
only (no cluster access needed).

---

## ⚠️ Currency flags (2026-08-15 audit, post-split outline: demographics → direct io-proc → imputed; no LM/ladder/exposures)

The io-proc thread's terminal state is 2026-06-29; nothing on it moved during the
July acceleration push, so "current" = the entry-39 finale.

1. **Direct joint fit CURRENT.** `joint_io_proc_lean_d2_enct_2k.summary.rds` (06-29,
   converged; entry 39) on `joint_io_proc_english_count_subset_data.rds` (06-28);
   `fig_io_imputed_proc.rds` (06-29) reads it. But `data_check_io_panels.R` (the
   observed-data figure) still reads the mm-era bundle (`joint_io_proc_mm_subset_data.rds`,
   06-16), not the english_count bundle the fit uses — repoint before submission.
2. **Separate INPUT model exists and is current-era but not wired.** The no-proc
   isolation fit (`log_irt_long_io_count.stan`, Sherlock job 31764182, 06-28, same
   bundle; input→accel 0.76 [.08,1.41], best-converged of the four) has no
   draws/summary in `fits/summaries/` — numbers live only in entry 39's table.
   Re-extract from Sherlock (or refit) and cache before it can go in main text.
3. **Separate PROC model STALE.** `proc_dp*_all` (06-10) predates the all-items
   switch (06-13), SEEDLings RT wiring (06-14), and its own rebuilt bundle (06-16).
   No proc-only fit exists in the lean/enct era; main-text use needs a refit
   (mirror of io_count: lean model with the input channel stripped).
4. **Old io-pooled STALE/superseded.** `io_pooled_widedelta` (06-02) vs bundle
   rebuilt 06-16; its `fig5_io_summary` headline is superseded by the joint +
   no-proc fits. Cut or refit.
5. **Demographics longitudinal arm: Marchman BLUPs are from the broken
   child_id extraction** (WG↔WS unlinked → WG arm lost; see
   `journal/notes/marchman_data_issue.md`). Clean re-extract exists
   (`01b_extract_marchman_clean.R`, 2060 kids) but **clean D_log was never fit**
   (needs Sherlock). With demographics now leading the results: refit on Sherlock
   and rebuild `blups_demographics.rds`, or drop Marchman from that arm.
   Cross-sectional arm (06-10 uncapped refit) is fine.
6. **Imputed panel: PRE-CLEANUP — bundles REBUILT clean 2026-08-15, refits pending.**
   The extraction port is done and the committed `long_subset_data[_nor].rds`
   now match the acceleration repo's clean bundles exactly (EN I=3105/N=4.25M
   = thal+smith+marchman summed; NO I=1630/N=4.44M; see
   `journal/notes/imputed_refit_plan.md` for validation + the six Sherlock
   launch lines). The six fits in `fits/summaries/long_no_freq_slopes*` are
   still the old-data ones until those refits land. Original defect record:
   (Upgraded 08-15 after comparing against the acceleration repo's
   `studies/bayes_long/00_prepare_bundles.R`, which corrects three defects this
   repo's `pull_longitudinal.R` → `prepare_longitudinal_{data,norwegian}.R`
   pipeline still has: (a) child key = Wordbank `child_id`, which fails to link
   WG↔WS — costs **Marchman's entire WG arm and ~178 Norwegian cross-form kids**,
   selectively the long-span, slope-informative trajectories; (b) no uni_lemma
   cross-form item harmonization; (c) no crater/jump QC filter.) The EN and NO
   D fits + all four σ_r anchors (6 GCP fits) sit on these bundles. Redo: port
   the clean extraction (key `dataset::study_internal_id` via
   `include_study_internal_id=TRUE`, uni_lemma option-a harmonization, clean_child
   QC), rebuild `long_subset_data[_nor].rds`, refit. QC alone moved acceleration's
   fits little (their SI sensitivity analysis), but the key fix is a
   sample-composition change — and the two companion papers cannot describe the
   same named datasets with different Ns and cleaning.
7. **Tables stale vs the enct bundle:** `table1_datasets.csv` (06-10) and
   `si_io_data_table.csv` (06-14, mm-era) both predate the english_count bundle;
   re-run their build scripts once the section list settles.
8. **Wordbank moved to Redivis (discovered 08-15).** The old MySQL connection
   endpoint (`wordbank.stanford.edu/db_args`) is gone — wordbankr ≤ 1.0.3
   cannot pull at all. Fix: `redivis` ≥ 0.12.12 (langcog r-universe) +
   `wordbankr` 2.0.0 (GitHub), which reads the versioned Redivis dataset
   `datapages.wordbank:627v` (current = **v2.0**; pass `version=` to pin) and
   needs `REDIVIS_API_TOKEN` (repo `.secrets`). ⚠️ Dataset v2.0 changed the
   instrument `value` column to raw responses ("yes"/"never"/…) — any script
   testing `value == "produces"` silently gets all-zeros; use the logical
   `produces` column. Fixed in `pull_longitudinal.R`; still latent in
   `precursor_variability_plot.R` / `precursor_irt_crosslang.R` (dormant), and
   `studies/glmer_ladder/01_extract_one.R`, `01b_extract_marchman_clean.R`,
   `studies/cross_sectional_demographics/00_build.R`, `paper/build_cache.R`
   pull Wordbank at runtime and need a 2.0 compatibility pass before their
   next rerun. Data-release drift caution: v2.0 contents ≠ the June MySQL
   snapshot the acceleration bundles were built from — per-dataset Ns may
   differ beyond the key-fix effect.

## ⚠️ Currency flags (open as of 2026-06-13)

1. **RESOLVED 2026-06-13 (partition): Fig 3 variance partition now sourced from the
   joint model.** `build_input_cache.R` reads `joint_io_proc_d3.draws.rds` for the
   `panel_partition` shares (in-model `share_input_xi`/`share_proc_xi` and `var_*_k`),
   replacing the old io-pooled + proc_dp stitch. `fig3_input.rds` rebuilt. **Still
   open:** the observed-input *fan* (panel C, `panelD`) and processing fan (panel D,
   `panelE`) still read `io_pooled_widedelta` / `proc_dp1_all` — swap to the joint
   model if full panel consistency is wanted.
2. **The wired joint fit is the DONE 5-study run; a 6-study version may be coming.**
   `joint_io_proc_d3.draws.rds` (06-12) is a 5-study / 292-kid fit (`tau_s[1..5]`).
   The on-disk bundle `joint_io_proc_subset_data.rds` has since grown to **S=6 / I=348**
   — so a 6-study fit is in progress/queued. The partition shares are self-contained in
   the draws, so the wired numbers are valid 5-study estimates; repoint when/if the
   6-study fit should drive Fig 3.
3. **`io-imputed D` inventory label says "local," but sampling was on GCP**
   (sm2-fit-01 EN / sm2-fit-02 NO; entries 32/34/35). "Local" refers to where the
   summaries/partition are extracted, not where MCMC ran.

---

## Table A — model fits behind the paper

| # | Model (conceptual) | Fit code + Stan/formula | Cluster | Trained | Raw fit outputs |
|---|---|---|---|---|---|
| **1** | **glmer ladder** — 7 nested mixed-effects accumulator models per dataset (A; B/C/D × linear/log age). Paper's M1–M4 = A/B/C/D_log; **M_best = D_log** (`(1+log_age \| child)`). Source of Fig 1B fits, AIC/BIC tables, and the longitudinal-arm BLUPs for demographics. | `studies/glmer_ladder/02_fit_one.R` (fit), `04a_simulate.R` (quantile fans → `sim_cache.rds`). `lme4::glmer`, `produces ~ 1 + log_age + (1+log_age\|child) + (1\|item)`. | Sherlock | 5 paper datasets (thal/smith/marchman/NO/JP) **2026-06-07**; `sim_cache` **06-10** | `fits/glmer_ladder/fit_<lang>_<model>.rds`, `summary_<lang>_<model>.csv`, `sim_cache.rds` |
| **2** | **io-imputed D** (Bayesian M_best, σ_r-pinned, *no* per-child input) — population share of efficiency variance attributable to input (π_α). Fig 3A + inline `io_partition`. | GCP family runner (`cluster/gcp/run_family.sh`, `wait_then_no_mbest.sh`), variant `long_no_freq_slopes[_norwegian]`. Stan `model/stan/log_irt_long.stan` (no-freq toggle). Entries 14, 24/25, 35. | GCP sm2-fit-01 (EN) / sm2-fit-02 (NO) | EN draws **06-09**; NO post-dedup refit collected **06-11** (entry 35) | `fits/summaries/long_no_freq_slopes[_norwegian].{draws,summary}.rds`, `_psi.csv`; bundles `long_subset_data[_nor].rds` |
| **3** | **io-pooled** (Bayesian, *observed* LENA/head-cam input, 4 datasets; wide-delta + γ slope channel) — observed input share of efficiency ≈ 2.8%. Fig 3 io panels + `fig5_io_summary`. | `model/scripts/fit_io_pooled_widedelta.R`; Stan `model/stan/log_irt_io_pooled.stan` (+ `_gamma_add`). Entry 30. | local | **2026-06-02/03** | `fits/io_pooled_widedelta.rds`, `io_pooled_gamma_widedelta_add.rds`, bundle `io_pooled_subset_data.rds` |
| **4** | **proc_dp D′0–D′3** (Bayesian processing+input regression ladder, 3 LWL datasets: AM2018/FM2012/FMW2013) — processing→efficiency; selected D′1. **Currently feeds Fig 3.** | `model/scripts/fit_proc_dp.R`; Stan `model/stan/log_irt_long_proc_dp.stan`. Entry 33. | Sherlock | **2026-06-09/10** | `fits/summaries/proc_dp{0,1,2,3}_all.{draws,summary,loo}.rds`, `proc_dp1_all_psi.csv`; bundle `proc_dp_all_subset_data.rds` |
| **5** | **joint io+proc D′0–D′3** (proc_dp + BabyView/SEEDLingS input-only; per-study σ_lena; σ_r = 0.44; 5 datasets / 292 kids) — **now drives the Fig 3 variance partition** (replaced the io-pooled + proc_dp stitch, 2026-06-13). | `model/scripts/fit_joint_io_proc.R`; Stan `model/stan/log_irt_long_proc_dp_joint.stan`; bundle `model/scripts/prepare_joint_io_proc_bundle.R`; `cluster/sherlock/joint_io_proc_fit.slurm`. Entry 37. | Sherlock | d3 (wired) **06-12**, 5-study; bundle since grown to 6-study (06-13, fit pending) | `fits/summaries/joint_io_proc_d<rung>.{draws,summary}.rds`; bundle `joint_io_proc_subset_data.rds` |
| **6** | **LLM acceleration** — GPT-2 (small) retrained on 24.5M-word CHILDES (training axis + distinct-input/development axis) + Chang & Bergen 2022 4-architecture fits; per-word sigmoid slopes vs children's κ. Fig `fig-llm-acceleration`. | `studies/llm/` (feng2024 protocol); see `journal/experiments_llm.md`. Children's κ from fit #2 draws. | Sherlock / Marlowe (GPU) | surprisal + sigmoids **06-10**; finer ladder **06-11** | `fits/llm/sigmoids/`, `fits/llm/ladder_bestval_finer.csv`, `fits/llm/surprisal_*.csv`, `data/chang_bergen_2022/` |
| **7** | **cross-sectional demographics** — per-language `glmer(produces ~ la*p + (1\|item) + (1\|child_id))` on 31 Wordbank languages; sex & maternal-ed → efficiency/acceleration. Fig `fig-demographics` (with the #1 longitudinal BLUPs as the longitudinal arm). | `studies/cross_sectional_demographics/00_build.R` (+ `composite_figure.R`). `lme4::glmer`, `nAGQ=0`. Entry 31. | local | uncapped refit **2026-06-11** | `studies/cross_sectional_demographics/cache/{fits,scatter}.rds` |

---

## Table B — cache layer (`paper/cache/` ← build script ← fit outputs → manuscript)

| `paper/cache/` file (setup var) | Built by | Reads (Table A #) | Manuscript consumer |
|---|---|---|---|
| `fig2_glmer_ladder.rds` (`ladder`) | `build_cache.R` §1 | #1 `sim_cache.rds` | `fig-schematic` panel B |
| `aic_summary.rds` (`aic`) | `build_cache.R` §1b | #1 `sim_cache.rds` | SI `tbl-aic`/`tbl-bic`; inline `aic_loglin_range` |
| `blups_demographics.rds` (`demo_blups`) | `build_cache.R` §3 | #1 `fit_*_D_log.rds` + Wordbank demographics | feeds #7 (longitudinal arm of `fig-demographics`) |
| `fig5_io_summary.rds` (`io_summary`) | `build_cache.R` §4 | #3 io-pooled wide-delta + γ | io-pooled headline (setup-loaded; verify live use) |
| `fig4_exposure.rds` (`exposure`) | `build_cache.R` §5 | #2 EN draws + psi + `long_subset_data.rds` | `fig-efficiency` |
| `fig6_llm_slopes.rds` (`llm_slopes`) | `build_cache.R` §6 | #6 LLM sigmoids/ladder + #2 EN/NO draws | `fig-llm-acceleration` |
| `io_partition.rds` (`io_part`) | `build_cache.R` §8 | #2 EN/NO summaries | inline "Population input-related variation" |
| `input_rate_table.rds` | `build_cache.R` §7 | CSVs (Sperry/HR + validation set) — *not a model fit* | SI `tbl-input-rates` |
| `fig3_input.rds` (`input3`) | `build_input_cache.R` | **partition: #5 joint** (06-13); fans: #3 io-pooled (panel C) + #4 proc_dp (panel D); #2 EN/NO (imputed rows + panels A/B) | `fig-io-partition` |
| `fig_io_imputed_proc.rds` | `build_input_cache.R` | #2 + #4 | SI panel (not in main setup load) |

---

## Notes / conventions

- **Bundles vs fits.** `*_subset_data.rds` files are *input bundles* (assembled CDI/RT/input
  data for Stan), not fits. A bundle newer than its fit summary ⇒ a re-fit is pending
  (see flag 2).
- **Sherlock fits** are sampled on the cluster, then summarized to `fits/summaries/<tag>.{draws,summary,loo}.rds`
  via `cluster/sherlock/extract_*.R`; the large per-chain CSVs are not committed (gc'd).
- **GCP fits** (#2) likewise land as extracted summaries; raw CSVs were garbage-collected
  on the VMs (entry 34).
- **The 5 "paper" glmer-ladder datasets** are thal/smith/marchman (English by-study) + NO + JP.
  Other `fit_*_D_log.rds` in `fits/glmer_ladder/` (english_american, finnish, french_*, spanish_mexican,
  dated 05-27) are extra languages **not** in the main ladder.
- **Authoritative narrative:** `journal/experiments.md` (numbered entries cited above) and its
  top "Analysis inventory" table; LLM detail in `journal/experiments_llm.md`.
