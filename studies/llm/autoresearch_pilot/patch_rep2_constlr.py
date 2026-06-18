#!/usr/bin/env python3
"""Bank variant: 2 passes over the corpus (repetition) at constant LR. Tests the
research-menu lead hypothesis that single-pass-over-fresh data is WHY returns
diminish -- repetition pushes toward the grokking / Critical-Data-Size regime,
which may flatten or flip the returns curve. Flat LR isolates repetition from the
schedule. Stop at epoch>=3 = two full passes (dispatcher sets CDI_SP_CAP=2500 so
the step-cap won't fire first). Run from a variant run dir (edits ./train.py)."""
import sys
P = "train.py"
s = open(P).read()
if "[rep2]" in s:
    print("already rep2"); sys.exit(0)
a = "lrm = get_lr_multiplier(progress)"
assert a in s, "lrm anchor not found"
s = s.replace(a, "lrm = 1.0  # [rep2] flat LR + two passes")
b = "epoch >= 2 or step"
assert b in s, "single-pass stop anchor not found"
s = s.replace(b, "epoch >= 3 or step")  # two full passes
open(P, "w").write(s)
print("patched: rep2_constlr (flat LR, 2 passes)")
