#!/usr/bin/env python3
"""Bank variant: constant learning rate (flat LR). Removes the schedule-induced
late-sharpening confound, giving the honest returns-to-experience curve. Applied
on top of the base (A40-adapted + probe + single-pass) train.py. Run from a variant
run dir (edits ./train.py)."""
import sys
P = "train.py"
s = open(P).read()
if "[constLR]" in s:
    print("already constLR"); sys.exit(0)
a = "lrm = get_lr_multiplier(progress)"
assert a in s, "lrm anchor not found"
s = s.replace(a, "lrm = 1.0  # [constLR] flat LR (no schedule-driven late sharpening)")
open(P, "w").write(s)
print("patched: constLR")
