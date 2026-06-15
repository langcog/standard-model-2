# io-proc measurement-model (_mm) GCP run — OPTIMIZED relaunch

**Launched 2026-06-15 04:07 UTC** on **sm2-fit-01** (us-central1-a, n2d-standard-128).
PID 3437; log `~/standard_model_2/mm_d3_v2.log`. Auto-stop watcher halts the VM
~3min after the fit exits (output saved first).

## What changed from the first (killed) run
- **Optimized model** `log_irt_long_proc_dp_joint_mm.stan`: reduce_sum over the CDI
  likelihood + precomputed admin_base[A]/item_offset[J] (collapses autodiff O(N)→O(A+J))
  + dropped unused log_lik GQ. Compiled `stan_threads=TRUE`.
- **Run config**: 4 chains × **30 threads = 120 of 128 cores**; grainsize=12809;
  1000 warmup + 1000 iter, adapt_delta 0.95. (First run used 4 cores, ~3-day ETA.)
- **New data** (bundle `joint_io_proc_mm_subset_data.rds`, I=377 J=681 N=768599
  V_obs=6498 N_lwl=2774): FMW2013 LENA now 18+24mo (61 kids, 47 w/ replicate).
  ⚠ ELENA CDI parsed but BLOCKED on id mismatch (pending Virginia crosswalk).

## GOTCHA / ops
- VM `.bashrc` CMDSTAN=2.36.0 (nonexistent) → `export CMDSTAN=/opt/cmdstan/cmdstan-2.38.0`.
- VM is NOT a git repo — sync by `gcloud compute scp` (model/bundle/driver).
- Output: `fits/summaries/joint_io_proc_mm_d3.{summary,draws}.rds` on the VM.

## Check / retrieve
```
gcloud compute ssh sm2-fit-01 --zone us-central1-a --command \
  'cd ~/standard_model_2; tail -15 mm_d3_v2.log; ls -la fits/summaries/joint_io_proc_mm_d3.*'
gcloud compute scp sm2-fit-01:'~/standard_model_2/fits/summaries/joint_io_proc_mm_d3.*' \
  fits/summaries/ --zone us-central1-a
```

## Look at: sigma_r posterior vs 0.44 pin; share_input_xi/share_proc_xi (now data-est);
## sigma_meas[1] head-cam (~0.70) vs [2] LENA (~0.28); sigma_rt0 ~0.143; divergences/rhat.
## Known: SEEDLings RT ~280ms short (button onset, absorbed by tau_s[6]); AM2018@30 RT
## bump is real but vanilla-filtered in peekbank -> absorbed.

---
## D'0-D'3 ladder on SHERLOCK (2026-06-15) — γ_in competition diagnostic
**Why:** the full glmer (`produces ~ input_z*log_age + (1+log_age|child) + (1|item)`,
406k obs, all 681 items) finds input → SLOPE significant: **input_z:log_age = 0.902
(p=0.002)**, ~4.7% of slope variance, positive in all 4 studies. The mm D'3 fit said
input→accel ~0.8% (γ_in=0.744, CI [-0.52,1.94]). The glmer-implied γ_in ≈ 2.5 is
OUTSIDE the D'3 CI → real structural disagreement, NOT mixing (r̂ 1.03 fine).
Hypothesis: γ_in competes with processing (β_k0,β_k1) + huge residual slope
(σ_ζ≈4.3) in D'3; D'0 pins the processing betas off, so γ_in is uncontested.

**Submitted:** `sbatch --array=0-3 cluster/sherlock/joint_io_proc_mm_fit.slurm`
→ array **29606743** (rungs D'0..D'3), 4 chains × 8 threads = 32 cores, 1000+1000.
Output: `$SCRATCH/standard_model_2/fits/summaries/joint_io_proc_mm_d{0,1,2,3}.{summary,draws}.rds`.

**Retrieve + compare when done:**
```
ssh sherlock 'squeue -u $USER --name=joint_io_proc_mm'   # check status
scp sherlock:'$SCRATCH/standard_model_2/fits/summaries/joint_io_proc_mm_d*.summary.rds' fits/summaries/
# then compare gamma_in (and gamma_in*sigma_r vs the glmer's 0.90) across D'0->D'3:
Rscript -e 'for(r in 0:3){s<-readRDS(sprintf("fits/summaries/joint_io_proc_mm_d%d.summary.rds",r));
  g<-s[s$variable=="gamma_in",]; sr<-s[s$variable=="sigma_r",];
  cat(sprintf("D%d: gamma_in=%.2f [%.2f,%.2f]  gamma_in*sigma_r=%.2f (glmer: 0.90)\n",
      r,g$mean,g$q5,g$q95,g$mean*sr$mean))}'
```
**Expectation test:** if γ_in (or γ_in·σ_r) jumps at D'0 vs D'3 → processing competition
confirmed. If γ_in stays ~0.7 at D'0 too → the residual slope (σ_ζ) is the absorber,
or the latent-vs-observed input gap is the issue (MCF's worry: processing isn't the culprit).
