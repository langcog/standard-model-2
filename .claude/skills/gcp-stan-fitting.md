---
name: gcp-stan-fitting
description: How to launch, monitor, and recover large Stan model fits on the standard_model_2 GCP VMs (sm2-fit-01, sm2-fit-02). Use whenever the user wants to fit a model, start a family run, check on a running fit, or diagnose a fit that didn't finish cleanly. Covers VM configuration, threading, compile flags, disk requirements, env vars, and the run_fit.sh / run_family.sh entry points.
---

# GCP Stan model fitting (standard_model_2)

## VMs in use

Two persistent VMs in `us-central1-a`, project `hs-babyview`:

- `sm2-fit-01` — English fits
- `sm2-fit-02` — Norwegian fits

Both are **`n2d-standard-128`**: 64 physical AMD EPYC cores (128 vCPUs with SMT), 128 GB RAM, **250 GB boot disk**.

The repo lives at `~/standard_model_2/` on each VM. Run all commands relative to that path.

## Threading rule (read before launching anything)

**4 chains × 16 threads/chain = 64 physical cores.** Do NOT set threads higher.

The machine has 128 vCPUs but only 64 physical cores (SMT 2x). Stan's HMC is FPU-bound, and SMT siblings share an FPU, so oversubscribing via threads causes ~3× slowdown. The first English run pre-dating this rule ran 4 chains × 8 threads on a 32-vCPU machine and was 3× slower than necessary.

Launch with: `STAN_THREADS_PER_CHAIN=16`.

## Compile flags (read before recompiling Stan)

The cmdstan `make/local` should contain **only** `STAN_THREADS=true`. **Do not add** `-march=native`, `-DNDEBUG`, or any extra optimization flags. These trigger "double free or corruption (out)" segfaults inside Stan 2.38 on AMD EPYC at default -O3.

If a fit is segfaulting at startup or mid-sample, suspect contamination from a previous bad-flags compile. Recover by:
```
cd $CMDSTAN_PATH && make clean-all && make build
```
Then recompile each model from its `.stan` source.

## Disk requirements

Stan writes log_lik for every (chain, iter, obs) tuple to /tmp CSV during sampling. **Measured, May 2026**: ~30 MB/iter per chain at N=2.18M (Norwegian). 4 chains × 1000 sampling iters = ~120 GB just for log_lik output.

**For SM2 VMs:**
- `sm2-fit-01` (English, N=1.1M): **250 GB** is correct.
- `sm2-fit-02` (Norwegian, N=2.18M): **400 GB**. The original 250 GB filled at iter ~1800/2000 of demo_alpha_norwegian and lost 15h of sampling. Resized 2026-05-19.

Online-resize a running VM:
```
gcloud compute disks resize <disk> --size=400GB --zone=us-central1-a
# Then on the VM:
sudo growpart /dev/sda 1 && sudo resize2fs /dev/sda1
```

If a fit is approaching disk-full mid-sample, **resize immediately** — the disk-fill failure mode is silent (Stan chain just dies with `iostream error`) and unrecoverable. The CSVs get cleaned up before you can rescue them.

## Standard launch env

For both EN and NO at the current model scale (I=500, N≈1–2M):

```bash
STAN_ITER=2000           # 1000 warmup + 1000 sampling
STAN_WARMUP=1000
STAN_THREADS_PER_CHAIN=16
# Don't set STAN_ADAPT_DELTA or STAN_MAX_TREEDEPTH unless a fit fails
# and you've confirmed it's not a cold-start fix (see stan-mixing-failures skill).
```

## Entry points

- **`gcp/run_fit.sh <variant> <dataset> [init_from_tag]`** — single fit + extracts.
- **`gcp/run_family.sh [english|norwegian]`** — runs the 5-variant family sequentially. Reads/honors `STAN_*` env vars. Skips variants whose `fits/<tag>.rds` already exists (so it doubles as a resume driver). Handles the cold-start vs warm-start policy automatically (see stan-mixing-failures).

Standard launch pattern (detached, survives SSH disconnect):

```bash
ssh sm2-fit-01:
cd ~/standard_model_2
STAN_ITER=2000 STAN_WARMUP=1000 STAN_THREADS_PER_CHAIN=16 \
  nohup bash gcp/run_family.sh english \
  > gcp_family_english.log 2>&1 < /dev/null & disown
```

The `< /dev/null & disown` is required for true detachment; without it, `nohup ... &` alone can still die when the SSH session drops.

## Monitoring a running fit

```bash
# Are chains alive and saturating CPU?
ps -eo etime,pcpu,cmd --sort=-pcpu | grep log_irt | grep -v grep | head

# What iter? (Stan logs every refresh=200)
tail -20 gcp_family_<lang>.log | grep "Iteration"

# Are extracts being written?
ls -lat fits/summaries/ | head
```

A healthy chain on `n2d-standard-128` should run at ~5–6 iter/min steady-state. Significantly slower (~1 iter/min) usually means treedepth saturation — see `stan-mixing-failures`.

## File layout (post-fit)

- `fits/<tag>.rds` — the cmdstanr fit object (large, 10–30 GB compressed)
- `fits/summaries/<tag>.summary.rds` — parameter posteriors
- `fits/summaries/<tag>.draws.rds` — scalar draws (small)
- `fits/summaries/<tag>_psi.csv` — per-item difficulty
- `fits/summaries/<tag>.loo.rds` — LOO output

Rsync summaries to local for analysis:
```
gcloud compute scp --zone=us-central1-a --recurse \
  sm2-fit-01:standard_model_2/fits/summaries/ ./fits/
```

## Recompiling a Stan model after edits

If you `scp` a modified `.stan` file to a VM, also `rm` the compiled binary before relaunching:
```
rm -f model/stan/log_irt_long          # the binary
rm -f model/stan/log_irt_long.rds      # cmdstanr's cache
```
Otherwise cmdstanr uses the stale binary.

**Mac → Linux warning:** never rsync the whole `model/stan/` directory from a Mac to a Linux VM — the Mac arm64 binary will clobber the Linux x86_64 binary and the next fit will segfault. Use selective `gcloud compute scp` of just the `.stan` source.
