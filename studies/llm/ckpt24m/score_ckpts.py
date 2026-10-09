"""Offline per-CDI-word surprisal for saved trajectory checkpoints, using the same
scoring code and settings as the in-training callback (max_ctx 128, 50 contexts/word).
    python score_ckpts.py --traj_dir RUN/trajectory --contexts cdi_contexts.jsonl --out traj.csv
Rows: step, epoch, word, n_occurrences, mean_nll, sum_nll (same schema as surprisal_*.csv).
"""
import argparse
import glob
import json
import os
import sys

import torch
from transformers import AutoModelForCausalLM

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "code"))  # deployed layout
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from surprisal_callback import WordSurprisalCallback

ap = argparse.ArgumentParser()
ap.add_argument("--traj_dir", required=True)
ap.add_argument("--contexts", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--steps", default="", help="comma-separated subset of steps (default all)")
a = ap.parse_args()

scorer = WordSurprisalCallback(a.contexts, a.out, eval_batch_size=128, max_ctx=128,
                               max_per_word=50)
scorer._init_csv()
want = {int(s) for s in a.steps.split(",")} if a.steps else None
for d in sorted(glob.glob(os.path.join(a.traj_dir, "step-*"))):
    st = json.load(open(os.path.join(d, "state.json")))
    if want is not None and st["step"] not in want:
        continue
    # bf16-stored weights, fp32 compute (the in-training callback scores in fp32)
    model = AutoModelForCausalLM.from_pretrained(d, torch_dtype=torch.bfloat16).float().cuda()
    scorer._evaluate(model, "cuda", step=st["step"], epoch=st["epoch"])
    print("scored", os.path.basename(d), flush=True)
