#!/usr/bin/env python3
"""Bank variant: model SIZE via DEPTH from env CDI_DEPTH (model_dim = depth*64,
rounded to 128; num_heads = model_dim/128). Tests whether scale changes the returns
curvature / the timing-sharpness of capability transitions. Base config otherwise
(decay, 1 pass). Run from a variant run dir (edits ./train.py)."""
import sys
P = "train.py"
s = open(P).read()
if "CDI_DEPTH" in s:
    print("already arch"); sys.exit(0)
a = "DEPTH = 8"
assert a in s, "DEPTH anchor not found"
s = s.replace(a, "DEPTH = int(os.environ.get('CDI_DEPTH', '8'))", 1)
open(P, "w").write(s)
print("patched: arch (CDI_DEPTH)")
