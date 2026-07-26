#!/bin/bash
# One GPU worker for the archival re-training sweep. Atomically claims (seed,rung)
# items (mkdir-based, no double-booking), trains with the EXACT original config,
# exports best ckpt -> bf16 -> HF, marks UPLOADED. Safe to add more workers later.
# Oak backup is a separate final bulk step (avoids DTN contention w/ babyview).
set -u
GPU=$1
ROOT=/data2/mcfrank/ladder
PY=$ROOT/condaenv/bin/python
FENG=$ROOT/standard-model-2/model/scripts/feng_eval
export HF_HOME=$ROOT/hf_cache TMPDIR=$ROOT/.tmp PYTHONDONTWRITEBYTECODE=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
REPO=mcxfrank/childes-gpt2-ladder
OUT=$ROOT/retrain_dev; LOCK=$OUT/.upload.lock
mkdir -p "$OUT/claims"
SEEDS=(0 1 2 3 4 5 7 42 99 123)
RUNGS=(0.5M 0.75M 1M 1.25M 1.5M 1.75M 2M 2.5M 3M 3.5M 4M 5M 6M 7M 8M 12M 16M 24M)
log(){ echo "[$(date +%T) GPU$GPU] $*" >> $OUT/dispatch.log; }
for seed in "${SEEDS[@]}"; do for tag in "${RUNGS[@]}"; do
  rd=$OUT/seed${seed}_${tag}
  [ -f "$rd/UPLOADED" ] && continue
  mkdir "$OUT/claims/seed${seed}_${tag}" 2>/dev/null || continue   # atomic claim
  tf=$ROOT/chunks/ladder_seed${seed}/CHILDES_ladder_seed${seed}_${tag}.txt
  [ -f "$tf" ] || { log "MISSING_CHUNK $tf"; continue; }
  mkdir -p "$rd"; rm -rf "$rd"/checkpoint-* 2>/dev/null
  log "START seed$seed $tag"
  CUDA_VISIBLE_DEVICES=$GPU $PY -B $FENG/train_gpt2_childes.py \
    --train_file "$tf" --val_file $ROOT/TinyDialogues/data/CHILDES_val_ordered.txt \
    --tokenizer_dir $ROOT/TinyDialogues/tokenizers/GPT2_CHILDES \
    --config_file $ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json \
    --output_dir "$rd" --cdi_contexts_jsonl $ROOT/cdi_contexts.jsonl \
    --surprisal_csv "$OUT/surprisal_seed${seed}_${tag}.csv" \
    --num_train_epochs 40 --per_device_batch_size 8 --learning_rate 1e-4 --seed $seed \
    --n_log_points 3 --early_stopping_patience 4 --save_total_limit 1 --max_eval_blocks 500 \
    --eval_callback_batch_size 128 --eval_max_ctx 128 --eval_max_per_word 50 > "$rd/run.out" 2>&1
  rc=$?
  [ $rc -ne 0 ] && { log "TRAIN_FAIL seed$seed $tag rc=$rc"; touch "$rd/FAILED_train_$rc"; rmdir "$OUT/claims/seed${seed}_${tag}" 2>/dev/null; continue; }
  ckpt=$(ls -d "$rd"/checkpoint-* 2>/dev/null | head -1)
  [ -z "$ckpt" ] && { log "NOCKPT seed$seed $tag"; touch "$rd/FAILED_nockpt"; rmdir "$OUT/claims/seed${seed}_${tag}" 2>/dev/null; continue; }
  ( flock 9
    CUDA_VISIBLE_DEVICES= $PY -B $ROOT/export_and_upload.py --ckpt "$ckpt" \
      --tokenizer_dir $ROOT/TinyDialogues/tokenizers/GPT2_CHILDES \
      --repo $REPO --path_in_repo development/seed${seed}/rung${tag} \
      --export_dir $OUT/export_${seed}_${tag} \
      --meta "{\"condition\":\"development\",\"seed\":$seed,\"budget_words\":\"$tag\",\"epochs_max\":40,\"early_stop_patience\":4,\"dtype\":\"bfloat16\"}" > "$rd/upload.out" 2>&1
  ) 9>"$LOCK"
  erc=$?
  [ $erc -ne 0 ] && { log "UPLOAD_FAIL seed$seed $tag rc=$erc"; touch "$rd/FAILED_upload_$erc"; rmdir "$OUT/claims/seed${seed}_${tag}" 2>/dev/null; continue; }
  rm -rf "$rd"/checkpoint-*; touch "$rd/UPLOADED"
  log "DONE seed$seed $tag"
done; done
log "worker $GPU exit"
