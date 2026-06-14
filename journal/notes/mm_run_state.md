# io-proc measurement-model (_mm) GCP run — in flight

**Launched 2026-06-14 19:14 UTC** on **sm2-fit-01** (us-central1-a, n2d-standard-128).

- **Model:** `model/stan/log_irt_long_proc_dp_joint_mm.stan` (input = measurement model,
  sigma_r ESTIMATED ~N(0.44,0.10), frank2026 RT priors). Rung **D′3** (all channels free).
- **Bundle:** `fits/joint_io_proc_mm_subset_data.rds` (I=368, J=681 all-items, N=761747,
  V_obs=6441 raw input recs, N_lwl=2774, both-channel 141). Built locally + scp'd
  (VM has NO git — sync by scp).
- **Driver:** `model/scripts/fit_joint_io_proc_mm.R 3`; 4 chains × 1000 warmup + 1000 iter,
  adapt_delta 0.95.
- **Log:** `~/standard_model_2/mm_d3.log` (PID 3882). **Output:**
  `fits/summaries/joint_io_proc_mm_d3.{summary,draws}.rds` on the VM.
- **GOTCHA:** the VM `.bashrc` sets `CMDSTAN=/opt/cmdstan/cmdstan-2.36.0` (does NOT exist).
  Must `export CMDSTAN=/opt/cmdstan/cmdstan-2.38.0` before any cmdstanr call. Model
  already compiled there.

## To check / retrieve
```
gcloud compute ssh sm2-fit-01 --zone us-central1-a --command \
  'cd ~/standard_model_2; tail -20 mm_d3.log; ls -la fits/summaries/joint_io_proc_mm_d3.*'
# when done:
gcloud compute scp sm2-fit-01:'~/standard_model_2/fits/summaries/joint_io_proc_mm_d3.*' \
  fits/summaries/ --zone us-central1-a
# then STOP the VM to stop billing:
gcloud compute instances stop sm2-fit-01 --zone us-central1-a
```

## Key things to look at in the result vs the current (pinned-sigma_r) D′3
- **sigma_r posterior** vs the 0.44 pin — does the data confirm/tighten it? CI width.
- **share_input_xi / share_proc_xi** — now both data-estimated (sigma_r free).
- **sigma_meas[1] (head-cam) vs [2] (LENA)** — instrument noise, should differ.
- **sigma_rt0** ~0.143 (frank prior should confirm, not move).
- divergences / rhat (sampler health, esp. the within-study centering + estimated sigma_r).
