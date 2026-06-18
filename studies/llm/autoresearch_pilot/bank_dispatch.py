#!/usr/bin/env python3
"""Fan a bank of train.py variants across FREE A40 GPUs (one run per GPU), each a
single-pass(+repetition) + CDI-probe run, for the acceleration-experiment hill-climb.

A variant = an (idempotent) patch script that edits ./train.py on top of the base
(already A40-adapted + probe + single-pass), plus optional per-variant env. Manifest:
    label                              # baseline (no patch)
    label : path/to/patch.py           # apply patch.py
    label : path/to/patch.py : K=v K=v # + per-variant env overrides
    label : : K=v K=v                  # env only, no patch
'#' comments and blank lines ignored.

Per variant -> runs_bank/<label>/{train.py, run.log, word_surprisal.csv}.
Only launches on GPUs with < FREE_MB used (co-exists politely; never double-books).
kappa/returns fitting is done afterward by pulling each word_surprisal.csv.

Usage (from the autoresearch repo dir):  python bank_dispatch.py bank_manifest.txt
"""
import os, sys, time, subprocess, shutil

REPO = "/data2/mcfrank/autoresearch"
PY = "/data2/mcfrank/ladder/condaenv/bin/python"
BANK = os.path.join(REPO, "runs_bank")
FREE_MB = 8000
STAGGER_S = 8
DEPS = ["prepare.py", "cdi_probe.py", "cdi_contexts.jsonl", "cdi_words.txt"]
ENV_BASE = dict(CDI_PROBE="1", CDI_PROBE_STEPS="5000", CDI_SP_STEPS="700",
                CDI_PROBE_NPOINTS="36", CDI_SP_CAP="6000",
                CDI_PROBE_CSV="word_surprisal.csv")


def free_gpus():
    out = subprocess.check_output(
        ["nvidia-smi", "--query-gpu=index,memory.used", "--format=csv,noheader,nounits"]).decode()
    g = []
    for line in out.strip().splitlines():
        idx, used = [x.strip() for x in line.split(",")]
        if int(used) < FREE_MB:
            g.append(int(idx))
    return g


def prep(label, patch):
    rdir = os.path.join(BANK, label)
    os.makedirs(rdir, exist_ok=True)
    shutil.copy(os.path.join(REPO, "train.py"), os.path.join(rdir, "train.py"))
    for f in DEPS:
        link = os.path.join(rdir, f)
        if not os.path.lexists(link):
            os.symlink(os.path.join(REPO, f), link)
    if patch:
        for one in patch.split(","):
            one = one.strip()
            if not one:
                continue
            p = one if os.path.isabs(one) else os.path.join(REPO, one)
            subprocess.check_call([PY, os.path.abspath(p)], cwd=rdir)
    return rdir


def launch(rdir, gpu, extra_env=None):
    env = {**os.environ, **ENV_BASE, **(extra_env or {}), "CUDA_VISIBLE_DEVICES": str(gpu)}
    log = open(os.path.join(rdir, "run.log"), "w")
    return subprocess.Popen([PY, "train.py"], cwd=rdir, env=env,
                            stdout=log, stderr=subprocess.STDOUT)


def parse_manifest(path):
    q = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.split(":")]
        lab = parts[0]
        patch = parts[1] if len(parts) > 1 and parts[1] else None
        env = {}
        if len(parts) > 2 and parts[2]:
            for kv in parts[2].split():
                k, v = kv.split("=", 1)
                env[k] = v
        q.append((lab, patch, env))
    return q


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: bank_dispatch.py <manifest>")
    queue = parse_manifest(sys.argv[1])
    os.makedirs(BANK, exist_ok=True)
    print(f"[bank] {len(queue)} variants: {[v[0] for v in queue]}", flush=True)

    running, done = {}, []
    while queue or running:
        for g in free_gpus():
            if g in running or not queue:
                continue
            lab, patch, venv = queue.pop(0)
            try:
                rdir = prep(lab, patch)
            except Exception as e:
                print(f"[skip] {lab}: prep failed: {e}", flush=True)
                done.append((lab, "prep_fail"))
                continue
            running[g] = (lab, launch(rdir, g, venv))
            print(f"[launch] {lab} -> GPU{g} (pid {running[g][1].pid}) env={venv}", flush=True)
            time.sleep(STAGGER_S)
        for g, (lab, p) in list(running.items()):
            if p.poll() is not None:
                print(f"[done] {lab} rc={p.returncode}", flush=True)
                done.append((lab, p.returncode))
                del running[g]
        time.sleep(20)
    print("BANK_COMPLETE", done, flush=True)


if __name__ == "__main__":
    main()
