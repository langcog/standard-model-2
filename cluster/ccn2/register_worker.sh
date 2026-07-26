#!/bin/bash
# Composition-control worker for one GPU. Atomic-claim over (dataset,seed,rung);
# EXACT L6 config; export bf16 -> HF; mark UPLOADED. Self-heal-friendly:
# writes gpu$GPU.pid, releases claim on transient failure, GAVEUP after 3 fails.
set -u
GPU=$1
ROOT=/data2/mcfrank/ladder; PY=$ROOT/condaenv/bin/python
FENG=$ROOT/standard-model-2/model/scripts/feng_eval
export HF_HOME=$ROOT/hf_cache TMPDIR=$ROOT/.tmp PYTHONDONTWRITEBYTECODE=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
REPO=mcxfrank/gpt2-composition-control
OUT=$ROOT/register_dev; LOCK=$OUT/.upload.lock
mkdir -p "$OUT/claims"; echo $$ > "$OUT/gpu${GPU}.pid"
DATASETS=(babylm_S1 babylm_S2 babylm_S3 climbmix_C1 climbmix_C2 climbmix_C3)
SEEDS=(0 1 2 3 4 5 7 42)
RUNGS=(0.5M 0.75M 1M 1.5M 2M 3M 4M 6M 8M 12M 16M 24M)
log(){ echo "[$(date +%m-%dT%T) GPU$GPU] $*" >> $OUT/dispatch.log; }
for ds in "${DATASETS[@]}"; do for seed in "${SEEDS[@]}"; do for tag in "${RUNGS[@]}"; do
  item="${ds}_seed${seed}_${tag}"; rd=$OUT/$item
  { [ -f "$rd/UPLOADED" ] || [ -f "$rd/GAVEUP" ]; } && continue
  mkdir "$OUT/claims/$item" 2>/dev/null || continue
  tf=$ROOT/register_chunks/${ds}_seed${seed}/${item}.txt
  [ -f "$tf" ] || { log "MISSING_CHUNK $tf"; touch "$rd/GAVEUP"; continue; }
  mkdir -p "$rd"; rm -rf "$rd"/checkpoint-* 2>/dev/null
  release(){ n=$(ls "$rd"/FAILED_* 2>/dev/null|wc -l); if [ "$n" -ge 3 ]; then log "GAVEUP $item ($n fails)"; touch "$rd/GAVEUP"; else rmdir "$OUT/claims/$item" 2>/dev/null; fi; }
  log "START $item"
  CUDA_VISIBLE_DEVICES=$GPU $PY -B $FENG/train_gpt2_childes.py \
    --train_file "$tf" --val_file $ROOT/TinyDialogues/data/CHILDES_val_ordered.txt \
    --tokenizer_dir $ROOT/TinyDialogues/tokenizers/GPT2_CHILDES \
    --config_file $ROOT/TinyDialogues/tokenizers/GPT2-small_config/config.json \
    --output_dir "$rd" --cdi_contexts_jsonl $ROOT/cdi_contexts.jsonl \
    --surprisal_csv "$OUT/surprisal_${item}.csv" \
    --num_train_epochs 40 --per_device_batch_size 8 --learning_rate 1e-4 --seed $seed \
    --n_log_points 3 --early_stopping_patience 4 --save_total_limit 1 --max_eval_blocks 500 \
    --eval_callback_batch_size 128 --eval_max_ctx 128 --eval_max_per_word 50 > "$rd/run.out" 2>&1
  rc=$?
  [ $rc -ne 0 ] && { log "TRAIN_FAIL $item rc=$rc"; touch "$rd/FAILED_train_$rc"; release; continue; }
  ckpt=$(ls -d "$rd"/checkpoint-* 2>/dev/null | head -1)
  [ -z "$ckpt" ] && { log "NOCKPT $item"; touch "$rd/FAILED_nockpt"; release; continue; }
  corpus=${ds%%_*}
  ( flock 9
    CUDA_VISIBLE_DEVICES= $PY -B $ROOT/export_and_upload.py --ckpt "$ckpt" \
      --tokenizer_dir $ROOT/TinyDialogues/tokenizers/GPT2_CHILDES --repo $REPO \
      --path_in_repo composition-control/${ds}/seed${seed}/rung${tag} --export_dir $OUT/export_${item} \
      --meta "{\"experiment\":\"composition-control\",\"corpus\":\"$corpus\",\"subset\":\"$ds\",\"seed\":$seed,\"budget_words\":\"$tag\",\"epochs_max\":40,\"early_stop_patience\":4,\"dtype\":\"bfloat16\"}" > "$rd/upload.out" 2>&1
  ) 9>"$LOCK"
  erc=$?
  [ $erc -ne 0 ] && { log "UPLOAD_FAIL $item rc=$erc"; touch "$rd/FAILED_upload_$erc"; release; continue; }
  rm -rf "$rd"/checkpoint-*; touch "$rd/UPLOADED"; rmdir "$OUT/claims/$item" 2>/dev/null
  log "DONE $item"
done; done; done
rm -f "$OUT/gpu${GPU}.pid"; log "worker $GPU exit"
