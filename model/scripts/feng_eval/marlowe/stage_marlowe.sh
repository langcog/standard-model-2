#!/usr/bin/env bash
# Stage the data-variance pilot on Marlowe (run on a LOGIN node, not in a job).
#
# Produces, under $PILOT_ROOT:
#   TinyDialogues/                     <- Feng et al. repo (CHILDES corpus + tokenizer)
#   TinyDialogues/data/CHILDES/CHILDES_{train,val}_ordered.txt
#   TinyDialogues/tokenizers/GPT2_CHILDES/, GPT2-small_config/config.json
#   standard-model-2/                  <- this repo (feng_eval pipeline)
#   cdi_contexts.jsonl                 <- eval set (val + 611 CDI words)
#   chunks/CHILDES_chunkA.txt, chunkB.txt, split_manifest.json
#   venv/                              <- python env
#
# EDIT THESE TWO before running:
PILOT_ROOT="${PILOT_ROOT:?set PILOT_ROOT to your Marlowe scratch path, e.g. /scratch/$USER/llm_var_pilot}"
# Python provisioning: this script uses a venv on top of a cluster python module.
# If Marlowe prefers conda/apptainer, swap the "Python env" block accordingly.

set -euo pipefail
mkdir -p "$PILOT_ROOT"
cd "$PILOT_ROOT"

echo "=== 1. Clone Feng TinyDialogues (CHILDES corpus + tokenizer) ==="
# CHILDES_data.zip is a Git-LFS object — install LFS first or it downloads as a pointer.
if [ ! -d TinyDialogues ]; then
  git lfs install
  git clone https://github.com/styfeng/TinyDialogues.git
  ( cd TinyDialogues && git lfs pull )
fi
# Unzip the corpus -> data/CHILDES/CHILDES_{train,val}_ordered.txt
if [ ! -f TinyDialogues/data/CHILDES/CHILDES_train_ordered.txt ]; then
  ( cd TinyDialogues/data && unzip -o CHILDES_data.zip )
fi
ls -la TinyDialogues/data/CHILDES/CHILDES_train_ordered.txt \
       TinyDialogues/data/CHILDES/CHILDES_val_ordered.txt \
       TinyDialogues/tokenizers/GPT2_CHILDES/tokenizer.json \
       TinyDialogues/tokenizers/GPT2-small_config/config.json

echo "=== 2. Clone standard-model-2 (feng_eval pipeline + cdi_words.txt) ==="
if [ ! -d standard-model-2 ]; then
  git clone https://github.com/langcog/standard-model-2.git
fi
FENG="$PILOT_ROOT/standard-model-2/model/scripts/feng_eval"
ls "$FENG/train_gpt2_childes.py" "$FENG/surprisal_callback.py" \
   "$FENG/extract_cdi_contexts.py" "$FENG/make_disjoint_chunks.py" \
   "$FENG/cdi_words.txt"

echo "=== 3. Python env (venv-on-module; adapt module names to Marlowe) ==="
# Marlowe modules per docs: slurm, nvhpc, cudnn/cuda12. Find the python module:
#   module avail python   (then load it below)
module load slurm 2>/dev/null || true
# EDIT: replace with Marlowe's actual python module name if not python/3.12
module load python/3.12 2>/dev/null || module load python 2>/dev/null || true
if [ ! -d venv ]; then
  python3 -m venv "$PILOT_ROOT/venv"
fi
source "$PILOT_ROOT/venv/bin/activate"
export PIP_CACHE_DIR="$PILOT_ROOT/.pipcache"; mkdir -p "$PIP_CACHE_DIR"
export TMPDIR="$PILOT_ROOT/.tmp"; mkdir -p "$TMPDIR"
pip install "pip<25"
# torch matched to Marlowe CUDA runtime (H100/cu12x). cu124 wheels work on H100.
pip install "torch==2.4.*" --index-url https://download.pytorch.org/whl/cu124
pip install --only-binary=:all: \
    "numpy>=2.0,<2.3" "scipy>=1.13" "pandas>=2.2,<3" "pyarrow>=14" \
    "transformers>=4.40,<4.50" "datasets>=2.20,<4" \
    "accelerate>=0.30,<2" "tokenizers>=0.15,<0.22"

echo "=== 4. HF caches to scratch ==="
export HF_HOME="$PILOT_ROOT/hf_cache"; mkdir -p "$HF_HOME"

echo "=== 5. Build the CDI-context eval set (shared by BOTH chunks) ==="
# Eval contexts come from the held-out VAL set, so they are independent of
# whichever training chunk a model sees -> no leakage asymmetry between chunks.
python "$FENG/extract_cdi_contexts.py" \
    --tokenizer "$PILOT_ROOT/TinyDialogues/tokenizers/GPT2_CHILDES" \
    --val_file "$PILOT_ROOT/TinyDialogues/data/CHILDES/CHILDES_val_ordered.txt" \
    --cdi_words "$FENG/cdi_words.txt" \
    --out_jsonl "$PILOT_ROOT/cdi_contexts.jsonl" \
    --out_coverage "$PILOT_ROOT/cdi_coverage.csv" \
    --max_per_word 200 --max_ctx 1024

echo "=== 6. Split CHILDES into two DISJOINT chunks ==="
mkdir -p "$PILOT_ROOT/chunks"
python "$FENG/make_disjoint_chunks.py" \
    --train_file "$PILOT_ROOT/TinyDialogues/data/CHILDES/CHILDES_train_ordered.txt" \
    --tokenizer_dir "$PILOT_ROOT/TinyDialogues/tokenizers/GPT2_CHILDES" \
    --out_a "$PILOT_ROOT/chunks/CHILDES_chunkA.txt" \
    --out_b "$PILOT_ROOT/chunks/CHILDES_chunkB.txt" \
    --target_tokens 19000000 \
    --mode random \
    --shuffle_seed 2026 \
    --manifest "$PILOT_ROOT/chunks/split_manifest.json"

echo "=== STAGED. Next: sbatch the training array (see train_gpt2_childes.slurm) ==="
