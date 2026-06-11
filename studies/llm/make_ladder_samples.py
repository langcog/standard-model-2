"""
Build a NESTED developmental ladder of CHILDES training samples for one seed.

Each seed is treated as an "individual": shuffle the conversations with the
seed, then emit cumulative prefixes that first reach each word-count target.
Larger budgets therefore CONTAIN smaller ones (a child accumulating input),
and "seed" bundles initialization + that individual's particular input stream
together — the total between-individual factor matched to children's σ_κ.

Budget is counted in WHITESPACE WORDS (the developmental currency, ~the CDS
word counts), so no tokenizer is needed here; the model still trains on the
fixed GPT2_CHILDES BPE (~1.9 tokens/word).

Usage:
  python make_ladder_samples.py \
      --train_file TinyDialogues/data/CHILDES_train_ordered.txt \
      --out_dir chunks/ladder_seed42 \
      --seed 42 \
      --targets_millions 1 2 4 8 16 24
"""

import argparse
import json
import os
import random


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_file", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--targets_millions", type=float, nargs="+",
                    default=[1, 2, 4, 8, 16, 24])
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.train_file) as f:
        lines = f.readlines()
    word_counts = [len(l.split()) for l in lines]
    total_words = sum(word_counts)

    order = list(range(len(lines)))
    random.Random(args.seed).shuffle(order)

    targets = sorted(int(t * 1_000_000) for t in args.targets_millions)
    manifest = {"seed": args.seed, "total_words": total_words,
                "total_conversations": len(lines), "rungs": []}

    cum_words = 0
    cum_idx = []           # conversation indices accumulated so far (shuffled order)
    ti = 0                 # next target index
    for pos, conv in enumerate(order):
        cum_idx.append(conv)
        cum_words += word_counts[conv]
        # emit every target that this step has now reached (handles big jumps)
        while ti < len(targets) and (cum_words >= targets[ti] or pos == len(order) - 1):
            tgt = targets[ti]
            mtag = f"{targets[ti] / 1_000_000:g}M"   # 0.5M, 1M, 1.5M, 24M
            out_path = os.path.join(args.out_dir, f"CHILDES_ladder_seed{args.seed}_{mtag}.txt")
            with open(out_path, "w") as wf:
                for i in cum_idx:
                    wf.write(lines[i])
            manifest["rungs"].append({
                "target_words": tgt, "tag": mtag,
                "actual_words": cum_words, "conversations": len(cum_idx),
                "file": out_path,
            })
            print(f"[ladder] {mtag}: {cum_words:,} words, {len(cum_idx):,} conversations -> {out_path}")
            ti += 1
        if ti >= len(targets):
            break

    # If the largest target exceeds the corpus, the final rung is the full set.
    with open(os.path.join(args.out_dir, f"ladder_manifest_seed{args.seed}.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[ladder] nested: each rung contains all smaller rungs. "
          f"corpus total = {total_words:,} words.")


if __name__ == "__main__":
    main()
