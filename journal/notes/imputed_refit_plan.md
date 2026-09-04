# Imputed-panel refits on clean bundles — runbook (2026-08-15)

**Why.** The EN/NO longitudinal bundles were rebuilt on the acceleration repo's
data standard (clean child key `dataset::study_internal_id`, uni_lemma
cross-form harmonization, crater/jump QC, wordbankr 2.0 / Redivis dataset
v2.0). All six imputed fits must be redone on them. Full context:
`journal/paper_models_provenance.md`, 08-15 currency flags.

**Compute: Sherlock only** (owners partition now available; GCP retired for
this. Local cmdstan is CrowdStrike-killed — never fit locally).

## The six fits — SUBMITTED 2026-08-15 (jobs 39218376–39218385, owners, 48h)

**σ_r anchor decision (MCF approved 08-15):** pins = **{0.35, 0.44, 0.58}**.
0.35 = the direct enct/io_count fits' own σ_r estimate (0.35 [0.32, 0.39],
study-centered → population lower bound); 0.44 = channel-matched population
re-anchor (main dashed pin); 0.58 = literature upper. **0.53 (Sperry CDS)
retired** — no baseline (un-suffixed) refits; `build_fig_io_cache.R` now
reads σ_xi from the 0p44 fits and the band is [0.32, 0.58], so it will
error (not silently go stale) until the new summaries land.

| tag | dataset | σ_r pin |
|---|---|---|
| `long_no_freq_slopes_sigmaR_0p35` | english | 0.35 |
| `long_no_freq_slopes_sigmaR_0p44` | english | 0.44 |
| `long_no_freq_slopes_sigmaR_0p58` | english | 0.58 |
| `long_no_freq_slopes_norwegian_sigmaR_0p35` | norwegian | 0.35 |
| `long_no_freq_slopes_norwegian_sigmaR_0p44` | norwegian | 0.44 |
| `long_no_freq_slopes_norwegian_sigmaR_0p58` | norwegian | 0.58 |

*Confirmed: the rebuilt bundles pin σ_r = 0.534 (Sperry-anchored baseline),
so the anchor set stays {0.44 pin, 0.53 baseline, 0.58 pin} and **all six
fits are needed** (build_fig_io_cache.R plots anchors at those three values
and skips missing summaries gracefully).

**Bundle validation (2026-08-15):** rebuilt EN = I 3105 / A 7329 / J 681 /
N 4,253,076 — exactly the sum of the acceleration repo's clean thal+smith+
marchman bundles (Marchman post-QC 2136 ✓, QC removed 90 admins ✓); rebuilt
NO = I 1630 / A 6446 / N 4,442,722 — exactly their norwegian bundle (QC 50 ✓).
The two repos are on identical data. EN a0 moved 19 → 18 (median admin age
with the WG arm restored); NO a0 = 25 unchanged.

## Launch pattern (already run 08-15; kept for resubmission)

2026-08-15 (later): everything is now committed and pushed
(`0562376` on master) and the Sherlock clone fast-forwarded to it with a
clean tree — **normal `git pull` is the sync path again**. Launch pattern:

```bash
for v in 0.35 0.44 0.58; do for d in english norwegian; do
  sbatch -p owners --requeue -t 48:00:00 \
    --export=ALL,STAN_SIGMA_R_OVERRIDE=$v \
    cluster/sherlock/long_fit.slurm no_freq_slopes $d
done; done
```

Be generous with `-t` (see feedback memory: short walltimes silently kill
fits; normal QOS caps at 48h, `--qos=long` beyond). The clean EN bundle is
~2× the old one (Marchman WG arm restored) — expect the EN fits to run
closer to the old Norwegian times than the old English times.

## ⚠️ 2026-08-16 — the OOM-at-the-finish-line problem, and the fix

Job **39218381 (NO σ_r 0.44)** sampled cleanly for 20h (all 4 chains
finished) and then was **OOM-killed at 48G in the post-sampling
`save_object()`** step, which reads every log_lik draw into memory. The
per-chain CSVs (4 × 28 GB) survived in `$SCRATCH/.../fits/csvs_<tag>/`.
Root cause: the Sherlock `long_fit.slurm` never set `STAN_SKIP_SAVE_OBJECT=1`
(the retired GCP runner did); with the clean bundles both languages are
"big" now. **The other five in-flight jobs will die the same way** — but
harmlessly, since the CSVs hold everything.

Fixes applied:
- `long_fit.slurm` now sets `STAN_SKIP_SAVE_OBJECT=1` and runs
  `recover_from_csvs.R <tag>` after the fit (streams `pi_alpha`, `sigma_xi`,
  etc. from the CSVs by column name — tiny RAM). Future launches are safe.
- NO 0.44: `recover_from_csvs.R` run manually on the login node
  (`logs/recover_no_0p44.log`).
- The five in-flight jobs: a Sherlock-side watcher
  (`cluster/sherlock/auto_recover_imputed_20260815.sh`, nohup on the login
  node, `logs/auto_recover.log`) polls squeue every 10 min and runs recovery
  for each job as it leaves the queue. Nothing depends on a local session.
- Recovery reads the CSV headers, so a partial/incomplete run would recover
  garbage — the watcher only fires after squeue drop, and `sacct` state is
  logged next to each recovery for a sanity check.

Expected outputs: `fits/summaries/<tag>.{summary,draws}.rds` per tag; pull
home with `scp sherlock:'$SCRATCH/standard_model_2/fits/summaries/long_no_freq_slopes*sigmaR*' fits/summaries/`.
Then delete the ~500 GB of CSVs on scratch (`csvs_no_freq_slopes*0p*`).

## 2026-08-17 — collection status + two findings

**Recovered so far** (pulled to local `fits/summaries/`; NB the tag lacks
the `long_` prefix — launched as variant `no_freq_slopes`, June's GCP runs
were `long_no_freq_slopes`, same Stan model; `build_fig_io_cache.R` expects
`long_`-prefixed tags → rename on the final collection):
- `no_freq_slopes_sigmaR_0p35` (EN 0.35): π_α 0.964, σ_ξ 1.85, δ 9.89,
  **σ_ζ 6.25** — but **σ_ξ/π_α ESS ≈ 22, r̂ 1.13: NOT converged.**
- `no_freq_slopes_norwegian_sigmaR_0p44` (NO 0.44): π_α 0.969, σ_ξ 2.49,
  δ 11.80, σ_ζ 7.48 — ESS 40, r̂ 1.10: marginal.
- NO 0.58: recovery v2 running on the 4 complete chains (the CSV dir held
  3 attempts' files; the stubs of preempted attempts broke the glob).
- EN 0.44 (39218380): chains at 60–90% sampling as of 08:20; expect the
  OOM-then-watcher-recovery path.
- EN 0.58 (39218382): **preempted three times on owners** (21 min, 13.6 h,
  then a 26 h attempt that OOM'd at the finish); now on attempt 4 from
  zero. NO 0.35 (39218379): OUT_OF_MEMORY at 19.75 h — but its 28 GB
  chains look complete; watcher should have recovered it (check log).

**Finding 1 — the clean data moved the numbers, in the expected
direction.** π_α EN 0.947 → 0.964 (input share of efficiency variance
5.3% → 3.6% at these pins); NO 0.974 → 0.969. **σ_ζ jumped a lot** (EN
4.41 → 6.25; NO 8.32 → 7.48): restoring the WG arm / long-span
cross-form kids adds real slope information — this is the sample-
composition effect predicted on 08-15. δ EN 10.67 → 9.89.

**Finding 2 — the imputed model mixes badly on the clean bundles.** ESS
≈ 20–40 on the headline scalars at 1000 sampling iterations. Options:
(a) accept as-is with honest r̂ reporting; (b) refit at 2000/2000 +
adapt_delta 0.97 (as was done for the enct joint fit); (c) treat these
as pilot values and rethink whether the imputed panel needs six full-size
fits at all (the analytic curve share = σ_r²/σ_ξ² only needs σ_ξ, which
IS well-determined in every arm — the pinned refits only add the anchor
CIs). Decision for MCF.

**Operational lesson: owners + 25 h chains + `--requeue` = restart
roulette.** cmdstan has no checkpointing; each preemption discards the
run. For any further imputed refits: `-p mcfrank` (owned node, no
preemption; walltimes generous per the cluster memory) or `-p normal
--qos=long`. Never owners for >~12 h fits.

## 2026-08-17 — DECISION (MCF): do the convergence refits. Plan + queue

All six imputed anchors get a **2000/1000, adapt_delta 0.97** refit
(`STAN_TAG_SFX=_2k`), warm-started via `STAN_INIT_FROM=<1000-iter tag>`
where a recovered summary exists. **All on `-p mcfrank`** (the owned node:
24 cores/192 GB, no preemption, no walltime cap) — `owners` preempts 25 h
chains into restart roulette, `normal` is hard-capped at 48 h and this
account has no `long` QOS (only `normal`, `high_p`). One 16-core fit at a
time → serial, ~40–50 h each → **allow ~1 week for the set**. No one needs
to babysit: `long_fit.slurm` now skips save_object and recovers the slim
summary from the CSVs inline.

| tag (`_2k`) | job | init |
|---|---|---|
| `no_freq_slopes_sigmaR_0p35_2k` | 39400957 (running 08-17 ~09:00) | warm from 0p35 |
| `no_freq_slopes_sigmaR_0p58_2k` | 39400960 | cold |
| `no_freq_slopes_norwegian_sigmaR_0p44_2k` | 39401017 | warm from NO 0p44 |
| `no_freq_slopes_norwegian_sigmaR_0p35_2k` | 39401019 | cold |
| `no_freq_slopes_norwegian_sigmaR_0p58_2k` | 39401025 | cold |
| `no_freq_slopes_sigmaR_0p44_2k` | auto-submitted by `cluster/sherlock/queue_en044_refit_20260817.sh` (nohup, login node) once the 1000-iter 39218380 finishes + recovers | warm if summary lands |

`fit_longitudinal.R` now honors `STAN_TAG_SFX` (it didn't; the first `_2k`
launch would have silently overwritten the 1000-iter CSV dir — caught and
relaunched). Two 1000-iter cold NO recoveries (0.35, 0.58) are also still
streaming on the login node — useful as sanity comparisons, not for the paper.

**Collection, when done:** pull `fits/summaries/no_freq_slopes*_2k.*`,
rename to the `long_`-prefixed tags `build_fig_io_cache.R` expects (or drop
the prefix there — simpler), rebuild `fig_io_imputed_proc.rds` + §8
`io_partition`, delete ~600 GB of `csvs_*` on scratch, re-render.

## 2026-09-03 — post-mortem: all six _2k refits SAMPLED, all six died after

Checked in after 2.5 weeks: sacct shows every `_2k` job OUT_OF_MEMORY. But
every `csvs_*_2k/` dir has 4 chains with cmdstan's "Elapsed Time" footer
and "All 4 chains finished successfully" in the logs (EN ~2d16h, NO
~1d14h each, serial on the owned node). **The sampling is done and safe
(1.2 TB on scratch).** They died in `fit_longitudinal.R`'s in-R
`summarize_fit()` — which loads every draw incl. log_lik — because
`STAN_SKIP_SAVE_OBJECT=1` skipped `save_object()` but not that step; then
`set -e` aborted the slurm before its inline recovery. Separately, the
08-17 login-node auto-recoveries of the 1000-iter fits failed because
cancelled mis-tagged relaunches had left empty partial CSVs in those dirs
("no lines available in input").

Fixed (commit on master, 09-03): flag also skips `summarize_fit`; slurm
recovers even when R exits non-zero; `recover_from_csvs.R` skips CSVs
without the footer and picks the latest complete run per dir; new
`recover_csvs.slurm` streams recovery on a compute node. **Recovery jobs
submitted 09-03 ~21:30 PT** (normal partition, 2 cores/16G/12h each):
41960236 EN .35 · 41960238 EN .44 · 41960240 EN .58 ·
41960299 NO .35 · 41960300 NO .44 · 41960305 NO .58. Expect
`fits/summaries/<tag>_2k.{summary,draws}.rds` within hours. Then: pull
home, reconcile the `long_` prefix in `build_fig_io_cache.R`, rebuild
caches, **delete ~1.9 TB of `csvs_*` on scratch**, re-render.

Lesson for the runbook: on these bundles nothing in R may touch the full
draws — every post-sampling path must stream from the CSVs by column.

Tag note: June's `long_no_freq_slopes*` and the new `no_freq_slopes*_2k`
are the same variant — `variant_hyperpriors()` strips the `long_` prefix
(`sub("^(long_proc_|long_|io_)", "", name)`); the June GCP runner just
passed the prefixed name. `build_fig_io_cache.R` now reads the `_2k` tags
directly (and asserts all six anchors are present).

## ✅ 2026-09-03 — DONE. All six `_2k` anchors recovered, caches rebuilt

Results and convergence assessment: `journal/experiments.md` #45 (🟢). σ_ξ EN 1.85 /
NO 2.49; band [0.32, 0.58] → input share of efficiency variance EN 3.0–9.9%,
NO 1.7–5.4%. `fig_io_imputed_proc.rds` + `io_partition.rds` rebuilt from the
`_2k` tags; manuscript re-rendered. Scratch: 1000-iter + proc_dp CSVs deleted
(653 GB); the six `_2k` CSV sets (1.2 TB) kept for one more pass — delete with
`rm -rf $SCRATCH/standard_model_2/fits/csvs_no_freq_slopes*_2k` once nobody
wants another column out of them. Do not refit these again: the residual
slow mixing is Norwegian's (α,ζ) ridge, not iteration count.

## After the fits land

1. Extract with the existing `cluster/sherlock/extract_summaries.R` flow →
   `fits/summaries/<tag>.{draws,summary}.rds`, pull back locally.
2. `Rscript paper/build_cache.R` (§8 io_partition; §7 input_rate_table
   unaffected) and `Rscript paper/build_fig_io_cache.R` (panel A curves +
   anchors; panel B unaffected — the joint enct fit is current).
3. Re-render `paper/input_paper.qmd`; the imputed section inline values
   (`io_share()`, π_α) update from the caches.

## Queue state (2026-08-15 evening)

Everything EXCEPT the six imputed refits above is now submitted or done —
the imputed refits are deliberately held pending the σ_r anchor-set
discussion (MCF: current anchors {0.44/0.53/0.58} may be wrong; note the
direct enct/io_count fits estimate σ_r ≈ 0.35, below the band's lower edge).

- Sherlock **39216615**: proc-only fit (`log_irt_long_proc_count.stan`, new,
  syntax-checked; mirror of io_count) → tag `joint_io_proc_lean_proc_d2_enct`.
- Sherlock **39216878**: glmer ladder, 35-cell array (5 clean datasets × 7
  models) on the rebased extraction — extracts match the acceleration
  bundles exactly (marchman 2136/QC90, smith 316, thal 653, norwegian
  1630/QC50; japanese 187 kids new-standard). Both on `-p owners --requeue`.
- io_count (no-proc) summary+draws recovered from Sherlock scratch →
  `fits/summaries/joint_io_proc_lean_io_d0_enct.*` (eff_input_k 0.762
  [0.082, 1.41], matches entry 39).
- `01_extract_one.R` REBASED onto `long_items.rds` (supersedes 01b);
  pull now includes Japanese; `data_check_io_panels.R` repointed to the
  english_count bundle; glmer slurm+submit paths fixed (were pre-reorg).

Remaining local chores after fits land: `paper/build_cache.R` §3 must join
demographics via the saved `ckey_map` in each `data_<slug>.rds` (ckey no
longer equals a Wordbank child_id) + a wordbankr-2.0 pass; wire an io_count/
proc_count cache once the separate-vs-joint figure design settles.

## Still queued behind this (separate work)

- **glmer ladder re-extraction** for the demographics BLUPs: recommended
  design = add Japanese to `pull_longitudinal.R` LANGUAGES and rebase
  `studies/glmer_ladder/01_extract_one.R` onto `long_items.rds` (one
  extraction, one standard; kills the child_id/QC divergence for good).
  Then Sherlock `glmer_ladder.slurm` refits ×5 → rebuild
  `blups_demographics.rds` (its build_cache.R §3 demo join must switch to
  ckey and get a wordbankr-2.0 pass).
- **Separate proc-only fit** (mirror of `log_irt_long_io_count.stan`) +
  re-extract the io_count no-proc summary from Sherlock (entry 39, job
  31764182) — unaffected by the Wordbank cleanup (io/proc datasets).
