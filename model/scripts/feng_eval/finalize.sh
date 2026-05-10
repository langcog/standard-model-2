#!/usr/bin/env bash
# End-to-end finalization after the 3-seed Sherlock training jobs complete.
#
# Pulls surprisal CSVs from Sherlock, fits per-word 4-PL sigmoids, runs the
# Feng vs. C&B comparison plot+CSV, and updates the report's headline table.
#
# Run from repo root:
#   bash model/scripts/feng_eval/finalize.sh
#
# Optional env:
#   SHERLOCK_HOST   ssh alias (default: sherlock)
#   FENG_REMOTE_DIR remote /scratch dir (default: /scratch/users/mcfrank/feng_eval)
#   SEEDS           space-separated seeds (default: "42 0 123")
#   PYTHON          python interpreter with numpy+scipy installed (default: tries
#                   /tmp/feng_local_venv/bin/python then `python3`)

set -euo pipefail

SHERLOCK_HOST="${SHERLOCK_HOST:-sherlock}"
FENG_REMOTE_DIR="${FENG_REMOTE_DIR:-/scratch/users/mcfrank/feng_eval}"
SEEDS="${SEEDS:-42 0 123}"

if [ -n "${PYTHON:-}" ]; then
  PY="$PYTHON"
elif [ -x /tmp/feng_local_venv/bin/python ]; then
  PY=/tmp/feng_local_venv/bin/python
else
  PY=$(command -v python3)
fi
echo "Python: $PY"
$PY -c "import numpy, scipy" || {
  echo "ERROR: numpy/scipy not available in $PY."
  echo "Either set PYTHON= to an interp with them installed, or create a venv:"
  echo "  python3 -m venv /tmp/feng_local_venv && /tmp/feng_local_venv/bin/pip install numpy scipy matplotlib"
  exit 2
}

mkdir -p data/feng_2026 outputs/feng_eval outputs/figs/longitudinal

# 1. Rsync surprisal CSVs from Sherlock
echo
echo "=== 1. Rsync surprisal CSVs from $SHERLOCK_HOST ==="
for seed in $SEEDS; do
  src="$SHERLOCK_HOST:$FENG_REMOTE_DIR/runs/surprisal_gpt2_childes_seed${seed}.csv"
  dst="outputs/feng_eval/surprisal_gpt2_childes_seed${seed}.csv"
  rsync -av "$src" "$dst" | tail -1 || {
    echo "WARN: could not pull seed $seed surprisal CSV; skipping that seed."
    continue
  }
  echo "  seed $seed: $(wc -l < "$dst") lines, $(cut -d, -f1 "$dst" | sort -un | wc -l) unique steps"
done

# 2. Fit per-word 4-PL sigmoids
echo
echo "=== 2. Fit per-word sigmoids ==="
for seed in $SEEDS; do
  csv="outputs/feng_eval/surprisal_gpt2_childes_seed${seed}.csv"
  tsv="data/feng_2026/gpt2_childes_seed${seed}_sigmoids.txt"
  [ -f "$csv" ] || { echo "  skip seed $seed (no CSV)"; continue; }
  $PY model/scripts/feng_eval/fit_per_word_sigmoid.py \
    --surprisal_csv "$csv" --out_tsv "$tsv"
done

# 3. Quick slope summary per seed
echo
echo "=== 3. Per-seed slope summary ==="
$PY - <<'PYEOF'
import csv, glob
for path in sorted(glob.glob("data/feng_2026/gpt2_childes_seed*_sigmoids.txt")):
    slopes = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                ps = float(r["ParamScale"])
                pu = float(r["ParamUpper"])
                pl = float(r["ParamLower"])
            except (KeyError, ValueError):
                continue
            if 0.01 < ps < 10 and (pu - pl) > 1.0:
                slopes.append(0.434 / ps)
    slopes.sort()
    n = len(slopes)
    if n == 0:
        print(f"  {path}: no kept words")
        continue
    print(f"  {path}: N={n} median={slopes[n//2]:.3f} "
          f"IQR=[{slopes[n//4]:.3f}, {slopes[3*n//4]:.3f}]")
PYEOF

# 4. Snapshot density plot (matplotlib, no R required)
echo
echo "=== 4. Snapshot density plot ==="
tsvs=()
labels=()
for seed in $SEEDS; do
  tsv="data/feng_2026/gpt2_childes_seed${seed}_sigmoids.txt"
  [ -f "$tsv" ] || continue
  tsvs+=("$tsv")
  labels+=("GPT-2 CHILDES seed $seed")
done
if [ ${#tsvs[@]} -gt 0 ]; then
  $PY model/scripts/feng_eval/make_partial_plot.py \
    --feng_tsvs "${tsvs[@]}" --feng_labels "${labels[@]}" \
    --out_png outputs/figs/longitudinal/feng_slope_comparison.png \
    --title_suffix ""
fi

# 5. Full R comparison script (uses real kid posterior draws if present)
echo
echo "=== 5. R comparison script (requires fits/summaries/long_no_freq_slopes*.draws.rds) ==="
if [ -f fits/summaries/long_no_freq_slopes.draws.rds ] && command -v Rscript >/dev/null 2>&1; then
  Rscript model/scripts/feng_eval/feng_chang_bergen_comparison.R || true
else
  echo "  Skipped: kid posterior draws not present locally, or Rscript not in PATH."
  echo "  Run on a machine with the Stan fits: Rscript model/scripts/feng_eval/feng_chang_bergen_comparison.R"
fi

echo
echo "=== Done. Artifacts ==="
ls -la data/feng_2026/*_sigmoids.txt 2>/dev/null | head
ls -la outputs/figs/longitudinal/feng_*comparison*.png 2>/dev/null | head
ls -la outputs/feng_eval/*.csv 2>/dev/null | head
echo
echo "Next: fill in the headline table in outputs/feng_evaluation_report.md using the"
echo "summary above, then commit + push as appropriate."
