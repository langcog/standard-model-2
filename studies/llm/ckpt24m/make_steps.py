"""Step grid for --traj_save_steps: step 0 + n log-spaced steps in [1, epochs*steps_per_epoch].
(Epoch ends are saved by the callback itself.) See SPEC.md.
    python make_steps.py --steps_per_epoch 5622 --epochs 12 --n 64 > steps_24M.txt
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from surprisal_callback import log_spaced_steps

ap = argparse.ArgumentParser()
ap.add_argument("--steps_per_epoch", type=int, required=True)
ap.add_argument("--epochs", type=int, default=12)
ap.add_argument("--n", type=int, default=64)
a = ap.parse_args()
print("\n".join(str(s) for s in [0] + log_spaced_steps(a.steps_per_epoch * a.epochs, a.n)))
