#!/usr/bin/env python3
"""Generic bank variant: override module-level train.py constant(s) from env
CDI_CONST = "NAME=value" or "NAME1=v1;NAME2=v2". Lets the overnight controller
sweep ANY hyperparameter (WARMUP_RATIO, WARMDOWN_RATIO, FINAL_LR_FRAC, WEIGHT_DECAY,
ADAM_BETAS, TOTAL_BATCH_SIZE, MATRIX_LR, EMBEDDING_LR, WINDOW_PATTERN, ASPECT_RATIO,
...) with no new patch code. Replaces the first module-level `NAME = <oldval>`.
Run from a variant run dir (edits ./train.py)."""
import os, re, sys
P = "train.py"
s = open(P).read()
spec = os.environ.get("CDI_CONST", "")
if not spec:
    print("no CDI_CONST set"); sys.exit(0)
for kv in spec.split(";"):
    kv = kv.strip()
    if not kv:
        continue
    name, val = kv.split("=", 1)
    name, val = name.strip(), val.strip()
    pat = re.compile(r"^(%s\s*=\s*)([^#\n]*)" % re.escape(name), re.M)
    if not pat.search(s):
        raise SystemExit("const not found: " + name)
    s = pat.sub(lambda m: m.group(1) + val + " ", s, count=1)
open(P, "w").write(s)
print("patched const:", spec)
