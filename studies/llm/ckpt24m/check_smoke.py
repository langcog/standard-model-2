"""Smoke-test gates for the trajectory-checkpoint run (SPEC.md). Run on pulled files:
    python check_smoke.py SMOKE_DIR PUBLISHED_CSV
SMOKE_DIR has surprisal_{A,B}.csv, traj_best.csv, A_saved.json, steps_smoke.txt.
"""
import json
import sys

import numpy as np
import pandas as pd

d, pub_csv = sys.argv[1], sys.argv[2]
SEED, RUNG, SPE = 42, "0.5M", 119  # steps per epoch at the 0.5M rung

def final(path):  # rows at the run's last logged step = the loaded best-val model
    x = pd.read_csv(path)
    return x[x.step == x.step.max()].set_index("word").mean_nll

def cmp(a, b, label):
    j = pd.concat([a, b], axis=1, join="inner").dropna()
    r, dm = np.corrcoef(j.iloc[:, 0], j.iloc[:, 1])[0, 1], (j.iloc[:, 0] - j.iloc[:, 1]).mean()
    print(f"{label:38s} n={len(j):4d}  r={r:.4f}  mean diff={dm:+.3f} nats")
    return r, abs(dm)

pub = pd.read_csv(pub_csv)
pub = pub[(pub.seed == SEED) & (pub.rung == RUNG)].set_index("word").surprisal
A, B = final(f"{d}/surprisal_A.csv"), final(f"{d}/surprisal_B.csv")
saved = json.load(open(f"{d}/A_saved.json"))
stop = max(h["step"] for h in saved["log_history"] if "step" in h)
best = int(saved["best_model_checkpoint"].rsplit("-", 1)[1])
grid = {int(s) for s in open(f"{d}/steps_smoke.txt").read().split()}
expect = {s for s in grid if s <= stop} | set(range(SPE, stop + 1, SPE))

ok = []
print(f"stop step {stop}, best step {best} (epoch {best / SPE:g})")
miss, extra = expect - set(saved["saved_steps"]), set(saved["saved_steps"]) - expect
print(f"gate 1  saved steps: {len(saved['saved_steps'])} saved, missing {sorted(miss)}, extra {sorted(extra)}")
ok.append(not miss and not extra)
r, dm = cmp(A, pub, "gate 2  A (with saving) vs published")
ok.append(r >= 0.98 and dm <= 0.05)
cmp(B, pub, "        B (control) vs published")
cmp(A, B, "        A vs B")
r, dm = cmp(final(f"{d}/traj_best.csv"), A, "gate 3  trajectory@best vs A in-run eval")
ok.append(r >= 0.999 and dm <= 0.01)
print("PASS" if all(ok) else f"FAIL {ok}")
