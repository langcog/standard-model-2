#!/usr/bin/env python3
"""Idempotently insert the CDI probe hooks into autoresearch train.py.
Three anchored inserts; probe is inert unless env CDI_PROBE is set, so the
baseline run (CDI_PROBE unset) is behaviorally identical to upstream."""
import sys

P = "train.py"
src = open(P).read()
if "CDI_PROBE" in src:
    print("train.py already patched")
    sys.exit(0)

INSERT_A = '''
# --- CDI per-word surprisal probe (active only if env CDI_PROBE is set) ---
probe = None
if os.environ.get("CDI_PROBE"):
    from cdi_probe import CdiProbe
    _here = os.path.dirname(os.path.abspath(__file__))
    probe = CdiProbe(
        os.path.join(_here, os.environ.get("CDI_PROBE_CONTEXTS", "cdi_contexts.jsonl")),
        os.environ.get("CDI_PROBE_CSV", "word_surprisal.csv"),
        total_steps=int(os.environ.get("CDI_PROBE_STEPS", "1000")),
        n_points=int(os.environ.get("CDI_PROBE_NPOINTS", "25")),
    )
'''

INSERT_B = '''    if probe is not None:
        probe.maybe(model, step, epoch)
'''

INSERT_C = '''
if probe is not None:
    probe.final(model, step, epoch)
'''


def insert_after(s, anchor, block):
    idx = s.find(anchor)
    if idx == -1:
        raise SystemExit(f"ANCHOR NOT FOUND: {anchor!r}")
    end = s.find("\n", idx) + 1
    return s[:end] + block + s[end:]


# ensure `import os`
if "import os" not in src:
    lines = src.splitlines(keepends=True)
    for i, l in enumerate(lines):
        if l.startswith("import ") or l.startswith("from "):
            lines.insert(i, "import os\n")
            break
    src = "".join(lines)

src = insert_after(src, "x, y, epoch = next(train_loader)  # prefetch first batch", INSERT_A)
src = insert_after(src, "    step += 1\n", INSERT_B)
src = insert_after(src, "    val_bpb = evaluate_bpb(model, tokenizer, DEVICE_BATCH_SIZE)", INSERT_C)

open(P, "w").write(src)
print("patched train.py (3 inserts)")
