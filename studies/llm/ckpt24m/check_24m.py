"""Per-seed gates for the 24M trajectory runs (SPEC.md). Run on ccn2 from ckpt24m/:
    python check_24m.py ladder_bestval_finer.csv
"""
import json
import sys

import numpy as np
import pandas as pd

pub = pd.read_csv(sys.argv[1])
grid = {int(s) for s in open("steps_24M.txt").read().split()}

def cmp(a, b):
    j = pd.concat([a, b], axis=1, join="inner").dropna()
    return len(j), np.corrcoef(j.iloc[:, 0], j.iloc[:, 1])[0, 1], (j.iloc[:, 0] - j.iloc[:, 1]).mean()

allok = True
for s in (42, 0, 123):
    R = f"runs/seed{s}"
    h = json.load(open(f"{R}/trajectory/log_history.json"))
    stop = max(x["step"] for x in h["log_history"] if "step" in x)
    spe = int(round(stop / max(x["epoch"] for x in h["log_history"] if "epoch" in x)))
    best = int(h["best_model_checkpoint"].rsplit("-", 1)[1])
    saved = set(h["saved_steps"])
    expect = {g for g in grid if g <= stop} | set(range(spe, stop + 1, spe))
    inrun = pd.read_csv(f"{R}/surprisal_inrun.csv")
    inrun = inrun[inrun.step == inrun.step.max()].set_index("word").mean_nll  # loaded best model
    traj = pd.read_csv(f"{R}/traj_surprisal.csv")
    tbest = traj[traj.step == best].set_index("word").mean_nll
    p = pub[(pub.seed == s) & (pub.rung == "24M")].set_index("word").surprisal
    n2, r2, d2 = cmp(inrun, p)
    n3, r3, d3 = cmp(tbest, inrun)
    ok = [expect <= saved, r2 >= 0.98 and abs(d2) <= 0.05, r3 >= 0.999 and abs(d3) <= 0.01]
    allok &= all(ok)
    print(f"seed{s}: stop {stop}, best {best} (epoch {best / spe:g}), {len(saved)} ckpts, "
          f"missing {sorted(expect - saved)}\n"
          f"  best vs published   n={n2} r={r2:.4f} diff={d2:+.3f}\n"
          f"  traj@best vs in-run n={n3} r={r3:.4f} diff={d3:+.4f}\n  {'PASS' if all(ok) else 'FAIL ' + str(ok)}")
print("ALL PASS" if allok else "SOME FAIL")
