#!/usr/bin/env python3
"""Bank variant: N passes (env CDI_PASSES) over the corpus with the TUNED decaying
schedule kept (no flat-LR override). Pair with patch_rep_constlr.py to test the
schedule x repetition interaction: does the decay schedule (which can manufacture
late sharpening) help repetition produce non-diminishing returns, vs flat LR?
Set CDI_SP_STEPS to the run length (~passes * 580) so the schedule spans the run.
Run from a variant run dir (edits ./train.py)."""
import sys
P = "train.py"
s = open(P).read()
if "CDI_PASSES" in s:
    print("already rep_decay"); sys.exit(0)
b = "epoch >= 2 or step"
assert b in s, "single-pass stop anchor not found"
s = s.replace(b, "epoch >= int(os.environ.get('CDI_PASSES', '1')) + 1 or step")
for k in ("torch.manual_seed(42)", "torch.cuda.manual_seed(42)"):
    s = s.replace(k, k.replace("42", "int(os.environ.get('CDI_SEED', '42'))"))
open(P, "w").write(s)
print("patched: rep_decay (decay schedule, CDI_PASSES passes, CDI_SEED)")
