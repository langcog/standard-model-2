#!/usr/bin/env python3
"""Switch autoresearch train.py from the 5-min wall-clock budget to a single pass
over the data. The dataloader's `epoch` is 1 during the first pass and increments
to 2 on wrap, so one full pass = stop when epoch >= 2 (plus a hard step-cap
backstop). Also rescale the LR/weight-decay `progress` from time-fraction to
step-fraction (step / CDI_SP_STEPS) so the schedule spans the whole pass.

Env: CDI_SP_STEPS = expected single-pass step count (~tokens/TOTAL_BATCH);
     CDI_SP_CAP   = hard safety ceiling on steps.
Idempotent. Run from the autoresearch repo dir, AFTER the probe patch."""
import sys

P = "train.py"
src = open(P).read()
if "CDI_SP_STEPS" in src:
    print("train.py already single-pass")
    sys.exit(0)

a = "    progress = min(total_training_time / TIME_BUDGET, 1.0)"
b = ('    progress = min(step / max(1, int(os.environ.get("CDI_SP_STEPS", "700"))), 1.0)'
     '  # [single-pass] step-fraction schedule')
assert a in src, "progress anchor not found"
src = src.replace(a, b)

c = "    if step > 10 and total_training_time >= TIME_BUDGET:"
d = ('    if step > 10 and (epoch >= 2 or step >= int(os.environ.get("CDI_SP_CAP", "1200"))):'
     '  # [single-pass] one full pass (epoch->2), with step-cap backstop')
assert c in src, "stop-condition anchor not found"
src = src.replace(c, d)

open(P, "w").write(src)
print("patched train.py for single-pass (epoch>=2 stop + step cap, step-fraction schedule)")
