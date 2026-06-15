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
