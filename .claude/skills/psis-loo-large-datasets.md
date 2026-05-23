---
name: psis-loo-large-datasets
description: Compute PSIS-LOO on Stan fits with very large log_lik matrices (N >~ 2M observations). Use whenever extract_loo_thinned.R segfaults inside psis_apply, or when the user is about to run LOO on a fit with N > 1.5M observations and wants to avoid the segfault. Explains the chunked-LOO recipe that ships in sherlock/extract_loo_thinned.R.
---

# PSIS-LOO for very large datasets

## When this matters

A standard `loo(ll, cores=1)` call on a log_lik matrix of `S` draws × `N` observations can segfault inside the `loo` package's `psis_apply` function when the matrix is too large. The threshold we hit:

- Norwegian bundle: N = 2.18M obs, 2000 draws (thinned from 4000)
- Resulting matrix: 35 GB
- Plus `lw_list` internal copies pushes peak memory past ~70 GB
- Segfault at `psis_apply(lw_list, "log_weights", fun_val = numeric(S))` with "memory not mapped"

This happened on `sm2-fit-02` (n2d-standard-128, 128 GB RAM). English (N=1.1M) is borderline OK, Norwegian is not.

## The fix is already in the repo

`sherlock/extract_loo_thinned.R` implements chunked PSIS-LOO. Just call it:
```bash
Rscript sherlock/extract_loo_thinned.R <tag>
```

The script:
1. Reads the `.rds`, extracts `log_lik` as draws-matrix
2. Thins to `N_THIN = 2000` draws if more
3. Computes `loo(...)` on column chunks of `CHUNK_OBS = 250000` observations each (~4 GB per chunk at 2000 draws)
4. Aggregates pointwise elpd / pareto_k / n_eff across chunks
5. Synthesizes a slim `loo` S3 object (skips `save_psis`, which is the largest field)

## Why chunking is exact, not approximate

PSIS-LOO computes for each observation `i`:
- `log_weights_i = -log_lik_i`
- Fits a generalized Pareto to the top tail (PSIS smoothing)
- `elpd_loo_i = log(mean(exp(log_lik_i + smoothed_log_weights_i)))`
- `pareto_k_i = GPD shape param`

These computations are **independent per observation**. Chunking by columns is therefore exact: total `elpd_loo = sum(elpd_loo_i)`, SE = `sqrt(N * var(elpd_loo_i))`. Validated locally on synthetic data: chunked output matches vanilla `loo()` to machine precision (max abs diff `0.000e+00` on elpd, pareto_k, mcse).

## Where the result lives

The script reads from `$STANDARD_MODEL_FITS_DIR` (default `$SCRATCH/standard_model_2/fits`) and writes to `$SCRATCH/standard_model_2/summaries/<tag>.loo.rds`. On GCP this is set up by `gcp/run_fit.sh` so the file lands in `fits/summaries/<tag>.loo.rds`.

## Symlink gotcha (run_fit.sh related)

The first version of `run_fit.sh` did `mkdir -p $SCRATCH/standard_model_2/summaries` before checking whether to make it a symlink to `fits/summaries`. This created a regular directory that swallowed all extract output, which never reached `fits/summaries/` for rsync.

`run_fit.sh` now detects this case explicitly: if `summaries/` exists as a real directory, it migrates contents into `fits/summaries/` and replaces it with a proper symlink. Read `run_fit.sh` if confused about where files landed.

## What success looks like

```
Reading /home/mcfrank/standard_model_2/fits/long_demo_pure_norwegian.rds ...
Extracting log_lik draws (memory-heavy step)...
  log_lik shape: 4000 draws x 2180667 obs (69.78 GB)
Thinned to 2000 draws (34.89 GB)
Computing PSIS-LOO in chunks of 250000 obs ...
  chunk 1/9: obs 1..250000 (250000 obs)
  ...
  chunk 9/9: obs 2000001..2180667 (180667 obs)
Wrote LOO (elpd = -654836.6 +- 770.1, pareto_k > 0.7: 97168 / 2180667 obs).
```

A high fraction of pareto_k > 0.7 (e.g., 4.5% of obs as above) is expected for an under-specified variant like `demo_pure`; it's a model-misspecification signal, not a numerical problem with the chunking.

## Don't try to "save_psis = TRUE"

The full `psis_object` for N=2M is ~100 GB; we don't save it. `chunked_loo` passes `save_psis = FALSE` to each chunk's `loo()` call, so the returned object is slim enough to fit in the `~/standard_model_2/fits/summaries/<tag>.loo.rds` file (typically <100 MB).

If a downstream consumer (e.g., `loo_compare`) needs the `psis_object`, it'll have to recompute it from the `log_lik` array.
