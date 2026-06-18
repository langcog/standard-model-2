#!/usr/bin/env python3
"""Bank variant: learning-rate MAGNITUDE — flat LR at CDI_LR_MULT x the base LR
(scales all param-group LRs via the loop multiplier). CDI_LR_MULT=1.0 == constlr.
Tests whether LR magnitude changes the returns curve. Run from a variant run dir."""
import sys
P = "train.py"
s = open(P).read()
if "CDI_LR_MULT" in s:
    print("already lr_mult"); sys.exit(0)
a = "lrm = get_lr_multiplier(progress)"
assert a in s, "lrm anchor not found"
s = s.replace(a, "lrm = float(os.environ.get('CDI_LR_MULT', '1.0'))  # [lrmult] flat LR x mult", 1)
open(P, "w").write(s)
print("patched: lr_mult")
