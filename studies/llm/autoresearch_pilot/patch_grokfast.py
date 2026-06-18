#!/usr/bin/env python3
"""Bank variant: Grokfast (Lee et al. 2024) -- EMA-amplify the slow (generalizing)
gradient component to accelerate grokking-style transitions. Tests whether pushing
on the generalizing direction changes the returns curve / acquisition dynamics.
lamb from env CDI_GF_LAMB (default 2.0). Applied on the base config (decay, 1 pass)
so it's directly comparable to baseline. Run from a variant run dir (edits ./train.py)."""
import sys
P = "train.py"
s = open(P).read()
if "[grokfast]" in s:
    print("already grokfast"); sys.exit(0)

DEF = ('''def _gf_ema(model, grads, alpha=0.98, lamb=2.0):
    if grads is None:
        return {n: p.grad.detach().clone() for n, p in model.named_parameters() if p.grad is not None}
    for n, p in model.named_parameters():
        if p.grad is not None and n in grads:
            grads[n].mul_(alpha).add_(p.grad.detach(), alpha=1 - alpha)
            p.grad.add_(grads[n], alpha=lamb)
    return grads
_gf_grads = None
_gf_lamb = float(os.environ.get("CDI_GF_LAMB", "2.0"))
_gf_alpha = float(os.environ.get("CDI_GF_ALPHA", "0.98"))  # [grokfast]

''')
a = "t_start_training = time.time()"
assert a in s, "loop-start anchor not found"
s = s.replace(a, DEF + a, 1)
b = "    optimizer.step()"
assert b in s, "optimizer.step anchor not found"
s = s.replace(b, "    _gf_grads = _gf_ema(model, _gf_grads, alpha=_gf_alpha, lamb=_gf_lamb)  # [grokfast]\n" + b, 1)
open(P, "w").write(s)
print("patched: grokfast")
