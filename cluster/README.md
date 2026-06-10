# cluster/ — compute helpers (job scripts + launchers)

Helpers for running fits on remote clusters. The analysis code these jobs invoke
lives in `model/` and `studies/`; these are the submit/launch/sync wrappers.

- **`gcp/`** — Google Cloud VMs (`sm2-fit-01/02`). `setup_vm.sh`, `run_fit.sh`,
  `run_family.sh`, queue/wait helpers. Used for the io-imputed D/D′ longitudinal fits.
- **`sherlock/`** — Stanford Sherlock SLURM. `*.slurm` job files + `extract_*.R`
  draw/summary extractors. Used for the glmer ladder, io_pooled, proc_dp, and GPT-2 fits.
  Manifest-driven array jobs: `glmer_ladder_submit.sh` builds `glmer_ladder_manifest.csv`.

**Marlowe (Stanford GPU cluster):** the only Marlowe helper is
`model/scripts/feng_eval/marlowe/stage_marlowe.sh`, kept with the LLM study pipeline
rather than moved here (see `studies/llm/`). Add a `cluster/marlowe/` dir if standalone
Marlowe launchers accumulate.

Paths inside these scripts assume the cluster-side layout (`$HOME/standard_model_2/…`,
`$SCRATCH/…`), not this repo's, so the 2026-06 reorg did not rewrite them. See `/MOVES.md`.
