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

---
## gamma_in prior-shrinkage test + lambda_bar (2026-06-15) — job 29719815
**Why:** the G0-G3 glmer morph FALSIFIED the "pinned efficiency suppresses acceleration"
hypothesis. Pinning input-efficiency to 0.358 (vs free 0.13, processing controlled)
did NOT drag input->accel (0.845 -> 0.964 if anything UP). So neither (A) latent-vs-observed
(cor 0.997, moot) nor (B) efficiency-pinning explains the SM2's low gamma_in. The cause is
Bayesian-specific: candidates = N(0,1) prior on gamma_in, and/or latent-kappa + large sigma_zeta
partitioning. See [[glmer_ladder_benchmark]].

**Submitted:** `sbatch --export=ALL,STAN_GAMMA_IN_PRIOR_SD=5 cluster/sherlock/joint_io_proc_mm_fit.slurm 3`
-> job **29719815** -> `fits/summaries/joint_io_proc_mm_d3_gp5.{summary,draws}.rds`.
Driver now also saves **sigma_lambda** -> lambda_bar = exp(sigma_lambda^2/2) for the theta->logit bridge.

**Read when done:**
- gamma_in (wide prior 5) vs the canonical 0.70 [N(0,1)]: if ~1.1 -> prior shrinkage was a factor;
  if ~0.7 -> structural (latent-kappa/sigma_zeta). Decomposition predicted likelihood center ~1.1.
- lambda_bar: SM2 input->accel on logit scale ~ gamma_in*sigma_r*lambda_bar; check vs glmer 0.85.
  NOTE variance SHARES (SM2 0.8% vs glmer ~4%) are scale-free and already disagree, so a real gap
  is expected regardless of lambda_bar.

### RESULTS (job 29719815, retrieved 2026-06-16)
| param | N(0,1) canonical | N(0,5) wide |
|---|---|---|
| gamma_in | 0.701 [-0.49, 1.89] | **1.688 [-0.33, 3.77]** (rhat 1.06) |
| input->accel (logit, lambda_bar=1) | 0.25 | **0.60** |
| input->accel SHARE | 0.8% | **3.0%** (glmer ~4%) |
| sigma_r | 0.358 | 0.355 |
| sigma_zeta | 4.302 | 4.277 |
| sigma_lambda | (not saved) | **0.046 -> lambda_bar=1.001** |

**Resolution:** (1) **lambda_bar=1.00** -> theta~logit, overlay scale is correct, no artifact.
(2) The **N(0,1) prior on gamma_in was the main culprit**: widening to N(0,5), gamma_in jumps
0.70->1.69 and the input->accel share 0.8%->3.0%, ~reconciling with the glmer benchmark (~4%).
So the "input barely touches acceleration" result was **mostly prior shrinkage, not a structural
failure** of the io model. No re-spec needed (G0-G3 already exonerated the coeff-1 efficiency term).
**Caveat:** gamma_in is weakly identified (CI [-0.33,3.77], rhat 1.06) because sigma_zeta~4.3
dominates the slope variance -- which is exactly why it's so prior-sensitive.

**Broader (MCF's sigma_zeta question):** sigma_zeta's prior is NOT the lever -- prior is
half-normal(0,1) but posterior is 4.30, data-dominated (robust). BUT all THREE slope-channel
coefs are weakly identified / prior-sensitive, not just gamma_in:
gamma_in [-0.49,1.89], beta_k0 [-1.78,1.47], beta_k1 [-1.98,1.30]. So prior sensitivity is a
property of the whole SLOPE channel (sigma_zeta soaks the variance). -> systematic sensitivity
sweep should target the slope-coef priors (+ sigma_r's lit prior). See [[glmer_ladder_benchmark]].

---
## RT measurement-model fix: mixed grain (SEEDLings trials vs Peekbank admins) — job 29864855
**Finding (MCF caught it):** the LWL RT channel mixes grains. 4 Peekbank `d_sub` datasets are
per-ADMIN (~1 row/run, within-run SD ~0.12-0.19); **SEEDLings (my Zhu derivation) is per-TRIAL**
(6.6 rows/run, within-run SD 0.43) and is 61% of all 2774 RT rows. So the SM2's single global
sigma_lwl was inflated by SEEDLings trial noise and mis-applied to Peekbank admins (under-IDing
their rt0), and SEEDLings kids got ~38 "independent" RT rows vs ~3-7 for Peekbank (effective-N
imbalance). This also broke my reliability calc (sub-obs treated as independent -> rho=0.16).

**Fix:** collapse SEEDLings trials -> per-(child,age) session means; keep Peekbank rows as-is.
2774 -> 1333 rows. Per-child proc reliability **0.16 -> 0.53**; sigma_child 0.063 -> 0.127.
Bundle: `fits/joint_io_proc_mm_runlvl_subset_data.rds` (built by /tmp/build_runlevel.R, post-proc).
Driver now takes `STAN_MM_BUNDLE` override; TAG gets `_runlvl` suffix.

**Submitted:** `sbatch --export=ALL,STAN_MM_BUNDLE=joint_io_proc_mm_runlvl_subset_data.rds
cluster/sherlock/joint_io_proc_mm_fit.slurm 3` -> job **29864855** ->
`joint_io_proc_mm_d3_runlvl.{summary,draws}.rds` (same N(0,1) priors as canonical D'3).

**Read when done:** compare run-level vs canonical D'3 -- does the processing channel sharpen?
beta_xi (proc->level), beta_k0/beta_k1 (proc->accel), sigma_rt0/sigma_rt1, and gamma_in.
Expect sigma_rt0 up (~0.13 vs 0.09), beta_xi better-identified. If so, the trial-noise
mis-spec was attenuating processing (compounding the N(0,1) prior shrinkage).

**FOLLOW-UP (circle back):** fix at source -- `prepare_seedlings_lwl_rt.R` should emit per-session
mean RT, not trials, so the canonical bundle is run-level. (Post-processed bundle is the test.)

### RESULTS (job 29864855, run-level D'3, retrieved 2026-06-16) — IT WORKED (level channel)
| param | canonical (trial LWL) | run-level LWL |
|---|---|---|
| sigma_lwl | 0.380 | **0.216** (trial noise removed) |
| sigma_rt0 | 0.090 | **0.133** (signal recovered; =glmer sigma_child 0.127) |
| beta_xi (proc->level) | -1.67 [-3.32, 0.09] | **-2.41 [-3.66, -1.11]** (CI now excludes 0) |
| beta_k0,beta_k1 (proc->accel) | -0.16,-0.33 | -0.09,-0.34 (unchanged) |
| gamma_in | 0.70 | 0.75 (unchanged) |
| sigma_zeta | 4.30 | 4.30 (unchanged) |
| proc->level per-SD | 0.15 | **0.32** (doubled, toward glmer ~0.58 / disatt ~0.80) |

**Verdict:** the SEEDLings trial-noise mis-spec was attenuating the LEVEL channel. Fixing the grain
recovers sigma_rt0 and sharpens beta_xi to a confident, ~2x larger effect. It does NOT touch the
SLOPE channel (gamma_in/beta_k0/beta_k1/sigma_zeta) -- that's the separate N(0,1)-prior +
sigma_zeta-domination problem. Two cleanly separated fixes: (1) measurement grain [DONE], (2)
standardized-per-SD slope-coef prior [TODO]. Residual proc->level gap (0.32 vs ~0.80) = beta_xi
still N(0,1)-shrunk.

---
## τ-SWEEP LAUNCHED (2026-06-16 night, YOLO) — reparam + new data
**Done this session:**
1. **Reparam** (commit 890f060): slope/level coefs are per-SD effects (a_in, b_xi, a_k0, a_k1)
   on standardized predictors; common-scale priors. STAN_TAU_SLOPE = sweep var. Compiles on Sherlock.
2. **SEEDLings RT grain fix** (caf9ab7): per-session means, not trials (prepare_seedlings_lwl_rt.R).
3. **BabyView refresh** (496ebfd): data_june_2026 -> 22 input kids, 116 admins (was 103).
4. **ELENA integrated** (5861ba9): ids reconciled (4943.. = LENA); parse_elena_cdi.R reads the WG
   xlsx; 26 kids, 16-18mo, median 44 words/kid, all 395 items match universe.
5. **Bundle rebuilt**: I=377->**403**, input 193->**219**, proc 326->**350**, both 142->**166**,
   N_lwl 1333->**1374**, J=681. Data-check figures refreshed (figs/data_checks/): input flat
   within-child (trait assumption holds); SEEDLings RT now session-level.

**Sweep:** `sbatch --export=ALL,STAN_TAU_SLOPE=<tau> joint_io_proc_mm_fit.slurm 3` for tau in
{0.5,1,2,4} -> jobs **29922069-72** -> `joint_io_proc_mm_d3{,_tau0.5,_tau2,_tau4}.summary.rds`
(tau=1 has NO suffix = canonical). All D'3, reparam, new bundle, N(0,1) priors superseded.

**Read when done (retrieve + compare a_in across tau vs glmer benchmark 0.85):**
```
scp sherlock:'$SCRATCH/standard_model_2/fits/summaries/joint_io_proc_mm_d3*.summary.rds' fits/summaries/
Rscript -e 'for(t in c("_tau0.5","","_tau2","_tau4")){f<-sprintf("fits/summaries/joint_io_proc_mm_d3%s.summary.rds",t)
  if(file.exists(f)){s<-readRDS(f); g<-function(v)s$mean[s$variable==v]
  cat(sprintf("tau%s: a_in=%.2f  a_k0=%.2f  a_k1=%.2f  b_xi=%.2f  share_in_k=%.1f%%\n", ifelse(t=="","1",t),
    g("a_in"),g("a_k0"),g("a_k1"),g("b_xi"), 100*g("var_input_k")/(g("var_input_k")+g("var_proc_k")+g("var_resid_k"))))}}'
```
**Expectation:** at tau=1 (weakly informative), a_in should land near the glmer's 0.85 -- i.e. the
honest common-scale prior ALONE recovers input->acceleration (old N(0,1) raw = N(0,0.36) per-SD
shrank it to 0.25). Sweep shows robustness. b_xi (proc->level) should hold ~its run-level value.
**Glmers NOT refit** (per MCF: data may change once more with next BabyView).

### τ-sweep D'3 RESULTS (jobs 29922069-72) — BROKE on a_k1 non-identifiability
| tau | a_in [q5,q95] | a_k1 | sigma_zeta | maxrhat | verdict |
|---|---|---|---|---|---|
| 0.5 | 0.39 [-0.13,0.92] | -0.76 | - | 1.42 | not converged |
| 1   | 0.47 [-0.15,1.07] | -2.21 | 3.64 | 1.07 | borderline (a_k1 contaminated) |
| 2   | 0.51 | -2.07 | - | 1.54 | not converged |
| 4   | 0.58 | -4.59 | 1.25 | 1.59 | GARBAGE (ess~7) |
**Diagnosis:** `a_k1` (rt1->accel) runs away and trades off with sigma_zeta -- per-child RT
SLOPES (rt1) are barely measured, so a_k1*z_rt1 is non-identifiable vs residual slope. The old
N(0,1)-on-raw-coef = N(0,0.16) per-SD MASKED this; the reparam's free prior exposed it & broke
the sampler. a_in (0.47 @ tau=1) is contaminated, not trustworthy.

### FIX: D'1 sweep (jobs 29980831-34) -- drop processing-SLOPE (a_k0,a_k1), null per glmer & unidentifiable
Keeps a_in (input->slope) + b_xi (proc->level, the real processing finding); processing->slope
pinned off. tau in {0.5,1,2,4} -> `joint_io_proc_mm_d1{,_tau0.5,_tau2,_tau4}.summary.rds`.
Read: a_in across tau vs glmer 0.85; b_xi should firm up (reparam unshrinks it too, ~disatt glmer 0.8).
