#!/usr/bin/env python3
"""Bank variant: initialization scale for the hidden matrices, via CDI_INIT_SCALE.
Scales s = sqrt(3)/sqrt(n_embd) (the uniform init of c_q/k/v, c_fc, value_embeds).
CDI_INIT_SCALE<1 = small init -> rich/feature-learning regime (saddle-to-saddle
stagewise learning, the top genuine increasing-returns candidate); >1 = lazy/NTK.
Run from a variant run dir (edits ./train.py)."""
import sys
P = "train.py"
s = open(P).read()
if "CDI_INIT_SCALE" in s:
    print("already initscale"); sys.exit(0)
a = "s = 3**0.5 * n_embd**-0.5"
assert a in s, "init-scale anchor not found"
s = s.replace(a, "s = 3**0.5 * n_embd**-0.5 * float(os.environ.get('CDI_INIT_SCALE', '1.0'))", 1)
open(P, "w").write(s)
print("patched: initscale (CDI_INIT_SCALE)")
