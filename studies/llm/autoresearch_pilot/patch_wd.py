#!/usr/bin/env python3
"""Bank variant: weight-decay multiplier (the #1 grokking lever per the research
menu). Scales the muon weight-decay schedule by env CDI_WD_MULT. Base config (decay,
1 pass) so it's comparable to baseline. Run from a variant run dir (edits ./train.py)."""
import sys
P = "train.py"
s = open(P).read()
if "[wdmult]" in s:
    print("already wd_mult"); sys.exit(0)
a = "muon_weight_decay = get_weight_decay(progress)"
assert a in s, "weight-decay anchor not found"
s = s.replace(a, "muon_weight_decay = get_weight_decay(progress) * float(os.environ.get('CDI_WD_MULT', '1.0'))  # [wdmult]")
open(P, "w").write(s)
print("patched: wd_mult")
