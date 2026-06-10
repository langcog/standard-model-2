# llm — large language model comparison  ·  Fig 5

**Provenance stub.** Code lives in `model/scripts/feng_eval/`, indexed here.

- **Pipeline:** `model/scripts/feng_eval/` (GPT-2 / CHILDES training + CDI surprisal scoring)
- **Marlowe staging:** `model/scripts/feng_eval/marlowe/stage_marlowe.sh` (Stanford Marlowe GPU cluster)
- **Cluster jobs:** `cluster/sherlock/feng_train_gpt2.slurm`, `feng_smoke.slurm`
- **Chang & Bergen comparison:** `model/scripts/chang_bergen_comparison.R`
- **Eval outputs:** `fits/feng_eval/` (surprisal CSVs, ladder summaries)
- **Slopes:** `fits/glmer_mbest_*` → `paper/cache/fig6_llm_slopes.rds`
- **Headline:** LLMs show standard accumulator dynamics — no acceleration, very limited variability.

Narrative: see [`/journal/experiments_llm.md`](../../journal/experiments_llm.md).
