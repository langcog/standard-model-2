# Sherlock deployment

This directory contains the files needed to run the fits on Stanford's
Sherlock HPC cluster via SLURM.

## One-time setup (on Sherlock)

```bash
# Log in
ssh <SUNetID>@login.sherlock.stanford.edu

# Clone this repo into $HOME
cd $HOME
git clone <your git URL> standard_model_2
cd standard_model_2

# Load R and install packages (interactive dev session)
ml R
sh_dev -c 4              # 4-core dev session; needed because rstan is compile-heavy
Rscript sherlock/setup_R.R
exit                     # back to login node
```

`setup_R.R` installs `rstan`, `posterior`, `dplyr`, `tidyr`, `ggplot2`,
`tibble`, `patchwork`, `MASS`, `arrow`, `wordbankr`, `childesr`,
`remotes`.

## Transferring data

**Do the Wordbank / CHILDES pulls LOCALLY on your laptop**, not on Sherlock.
Two reasons:

1. `wordbankr` depends on `RMySQL`, which needs mysql headers not available on
   Sherlock.
2. `childesr` requires R ≥ 4.4; Sherlock's R module is 4.2.

The pulled intermediate files (`model/fits/long_ws_items.rds`,
`model/fits/norwegian_word_freq.rds`) are committed to the repo, so on
Sherlock you just do:

```bash
# Prepare the Stan-ready bundles (reads the committed intermediate files):
cd $HOME/standard_model_2
ml R
Rscript model/scripts/prepare_longitudinal_data.R 'English (American)' 500 1000
Rscript model/scripts/prepare_longitudinal_norwegian.R 500 1000
```

If you need to refresh the intermediate files (e.g., new Wordbank release),
re-run `pull_longitudinal.R` and `pull_norwegian_freq.R` on your laptop,
`git add model/fits/long_ws_items.rds model/fits/norwegian_word_freq.rds`,
commit, push, and `git pull` on Sherlock.

## Submitting a fit

One SLURM script covers all longitudinal variants; pass variant name as
a job argument:

```bash
cd $HOME/standard_model_2
# English longitudinal 2PL+slopes
sbatch sherlock/long_fit.slurm long_2pl_slopes

# Norwegian longitudinal 2PL+slopes
sbatch sherlock/long_fit.slurm long_2pl_slopes_nor
```

Output `.rds` lands in `$SCRATCH/standard_model_2/fits/`, logs in
`$SCRATCH/standard_model_2/logs/`.

For the cross-sectional sensitivity sweep, each σ_r value is its own
job; submit them in parallel:

```bash
for sr in 0.3 0.53 0.8 1.2; do
  sbatch sherlock/sensitivity_fit.slurm 2pl $sr
done
```

## Pulling results back

```bash
# Locally
rsync -avz <SUNetID>@login.sherlock.stanford.edu:"\$SCRATCH/standard_model_2/fits/" model/fits/
rsync -avz <SUNetID>@login.sherlock.stanford.edu:"\$SCRATCH/standard_model_2/logs/" sherlock/logs/
```

Then run analysis / plotting locally via the existing `model/scripts/`.

## Design notes

- The 4 chains run within one SLURM job on 4 cores (not 4 SLURM tasks).
  `options(mc.cores = 4)` is already set in `model/R/helpers.R`.
- Memory: 16 GB is plenty for the cross-sectional, 32 GB recommended
  for the longitudinal fits. Set in `.slurm` files.
- Wall time: set generously (see table in the parent project's README).
- Partition: `normal` for jobs < 48 hr; no need for `long`.
