#!/usr/bin/env python3
"""Bank variant: flat LR + N passes over the corpus, N from env CDI_PASSES.
Maps the repetition-depth axis at constant LR (the lead returns-hypothesis: does
repetition push the learner out of the smooth diminishing-returns regime toward
the grokking / Critical-Data-Size regime?). Stop at epoch >= CDI_PASSES+1.
Run from a variant run dir (edits ./train.py)."""
import sys
P = "train.py"
s = open(P).read()
if "[repconst]" in s:
    print("already rep_constlr"); sys.exit(0)
a = "lrm = get_lr_multiplier(progress)"
assert a in s, "lrm anchor not found"
s = s.replace(a, "lrm = 1.0  # [repconst] flat LR")
b = "epoch >= 2 or step"
assert b in s, "single-pass stop anchor not found"
s = s.replace(b, "epoch >= int(os.environ.get('CDI_PASSES', '1')) + 1 or step")
for k in ("torch.manual_seed(42)", "torch.cuda.manual_seed(42)"):
    s = s.replace(k, k.replace("42", "int(os.environ.get('CDI_SEED', '42'))"))
open(P, "w").write(s)
print("patched: rep_constlr (flat LR, CDI_PASSES passes, CDI_SEED)")
