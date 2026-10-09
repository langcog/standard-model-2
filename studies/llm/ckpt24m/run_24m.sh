#!/bin/bash
# 24M-rung retrain with trajectory checkpoints (SPEC.md): one seed per GPU, L6 config
# verbatim + --traj_save_steps, then offline per-word surprisal for every checkpoint.
# Writes runs/seed{S}/{RC,SCORED_RC} and ALL_DONE. HF upload is a separate, gated step.
ROOT=/data2/mcfrank/ladder; PY=$ROOT/condaenv/bin/python; H=$ROOT/ckpt24m
export HF_HOME=$ROOT/hf_cache TMPDIR=$ROOT/.tmp PYTHONDONTWRITEBYTECODE=1 HF_HUB_OFFLINE=1 HF_DATASETS_OFFLINE=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
PYTHONPATH=$H/code $PY $H/make_steps.py --steps_per_epoch 5622 --epochs 12 --n 64 > $H/steps_24M.txt
run(){ GPU=$1; S=$2; R=$H/runs/seed$S; mkdir -p $R
  CUDA_VISIBLE_DEVICES=$GPU $PY -B $H/code/train_gpt2_childes.py \
    --train_file $ROOT/chunks/ladder_seed$S/CHILDES_ladder_seed${S}_24M.txt \
    --val_file $ROOT/TinyDialogues/data/CHILDES_val_ordered.txt \
    --tokenizer_dir $ROOT/TinyDialogues/tokenizers/GPT2_CHILDES \
    --config_file $ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json \
    --output_dir $R --cdi_contexts_jsonl $ROOT/cdi_contexts.jsonl \
    --surprisal_csv $R/surprisal_inrun.csv \
    --num_train_epochs 40 --per_device_batch_size 8 --learning_rate 1e-4 --seed $S \
    --n_log_points 3 --early_stopping_patience 4 --save_total_limit 1 --max_eval_blocks 500 \
    --eval_callback_batch_size 128 --eval_max_ctx 128 --eval_max_per_word 50 \
    --traj_save_steps $H/steps_24M.txt > $R/run.out 2>&1
  echo $? > $R/RC
  CUDA_VISIBLE_DEVICES=$GPU $PY -B $H/score_ckpts.py --traj_dir $R/trajectory \
    --contexts $ROOT/cdi_contexts.jsonl --out $R/traj_surprisal.csv > $R/score.out 2>&1
  echo $? > $R/SCORED_RC; }
run 0 42 &
run 1 0 &
run 2 123 &
wait; touch $H/ALL_DONE
