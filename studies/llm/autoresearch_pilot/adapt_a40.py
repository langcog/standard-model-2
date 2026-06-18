#!/usr/bin/env python3
"""Adapt autoresearch train.py to run on Ampere (A40).

Stock uses FlashAttention-3 (Hopper sm_90 only) via HF `kernels`. Replace with
PyTorch scaled_dot_product_attention. Because the whole model is torch.compiled,
we force the non-materializing SDP backends GLOBALLY (disable the math backend,
which inductor would otherwise lower to a full B*T*T attention matrix and OOM the
A40). Architecture/math identical; per-layer sliding window_size ignored (full
causal). Also shrink the micro-batch (TOTAL_BATCH_SIZE unchanged -> training
identical, just more grad-accum) for headroom on the 45 GB card. Idempotent.
"""
import sys

P = "train.py"
lines = open(P).read().split("\n")
if any("enable_math_sdp" in l for l in lines):
    print("train.py already adapted for A40")
    sys.exit(0)

out, swapped, flagged = [], False, False
for l in lines:
    s = l.strip()
    if (s.startswith("from kernels import")
            or (s.startswith("repo =") and "flash-attn" in s)
            or s.startswith("fa3 = get_kernel")):
        out.append("# [A40 adapt] removed (Hopper-only): " + s)
        continue
    if "fa3.flash_attn_func" in l:
        ind = l[:len(l) - len(l.lstrip())]
        out += [
            ind + "# [A40 adapt] SDPA replaces FA3; contiguous, full causal (window_size ignored).",
            ind + "q = q.transpose(1, 2).contiguous()",
            ind + "k = k.transpose(1, 2).contiguous()",
            ind + "v = v.transpose(1, 2).contiguous()",
            ind + "if self.n_kv_head != self.n_head:",
            ind + "    _rep = self.n_head // self.n_kv_head",
            ind + "    k = k.repeat_interleave(_rep, dim=1)",
            ind + "    v = v.repeat_interleave(_rep, dim=1)",
            ind + "y = F.scaled_dot_product_attention(q, k, v, is_causal=True)",
            ind + "y = y.transpose(1, 2)",
        ]
        swapped = True
        continue
    out.append(l)
    if not flagged and s.startswith("import torch.nn.functional"):
        out += [
            "# [A40 adapt] force non-materializing SDP backends (FA3 was Hopper-only)",
            "torch.backends.cuda.enable_flash_sdp(True)",
            "torch.backends.cuda.enable_mem_efficient_sdp(True)",
            "torch.backends.cuda.enable_math_sdp(False)",
        ]
        flagged = True

src = "\n".join(out)
src = src.replace("DEVICE_BATCH_SIZE = 128",
                  "DEVICE_BATCH_SIZE = 64  # [A40 adapt] was 128 (45GB OOM)")
assert swapped, "did not find fa3.flash_attn_func call"
assert flagged, "did not find functional import anchor"
assert "DEVICE_BATCH_SIZE = 64" in src, "batch reduction not applied"
open(P, "w").write(src)
print("adapted train.py for A40 (SDPA no-math backend + batch 64)")
