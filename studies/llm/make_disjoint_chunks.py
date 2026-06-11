"""
Split CHILDES_train_ordered.txt into two DISJOINT training chunks for the
data-variance pilot.

The data-variance question is: does the *identity* of the training data move
the per-word Chang & Bergen sigmoid slope, holding architecture, seed, tokenizer,
and eval set fixed? To answer it we need two training corpora that share no
conversations. This script produces them.

Unit of measurement is GPT2_CHILDES BPE tokens (so the two chunks have an
~equal number of gradient steps per epoch — "training amount" held constant).
We assign whole lines (= whole CHILDES conversations) so we never split a
conversation across chunks.

Modes:
  contiguous (default): chunk A = the first `target` tokens' worth of lines,
      chunk B = the LAST `target` tokens' worth of lines. The middle is left
      out, guaranteeing a clean gap so no boundary conversation is shared.
      Per-chunk line ORDER is preserved (matches Feng et al.'s ordered
      no-shuffle curriculum). Interpretation: "two different developmental
      slices of CHILDES" — the closest corpus analog to two children hearing
      systematically different input.
  random: shuffle conversations with a fixed seed, then take two disjoint
      groups of `target` tokens each. Removes the early-vs-late content
      confound; estimates pure subsample variance. Within-chunk order is the
      shuffled order (departs from the ordered curriculum).

Usage:
  python make_disjoint_chunks.py \
      --train_file TinyDialogues/data/CHILDES/CHILDES_train_ordered.txt \
      --tokenizer_dir TinyDialogues/tokenizers/GPT2_CHILDES \
      --out_a CHILDES_chunkA.txt \
      --out_b CHILDES_chunkB.txt \
      --target_tokens 19000000 \
      --mode contiguous

Default target_tokens=19_000_000 BPE ≈ 10M words ≈ the CDI-completion anchor
(GPT2_CHILDES fertility is ~1.9 BPE tokens/word). For two equal disjoint
HALVES of the corpus instead, pass --target_tokens 0 (uses all data, splits
at the token midpoint).
"""

import argparse
import json
import random
import sys

from transformers import AutoTokenizer


def count_line_tokens(tok, lines):
    """Return a list of BPE token counts, one per line. Batched for speed."""
    counts = []
    B = 1000
    for i in range(0, len(lines), B):
        batch = [l.rstrip("\n") for l in lines[i:i + B]]
        enc = tok(batch, add_special_tokens=False)["input_ids"]
        counts.extend(len(ids) for ids in enc)
        if (i // B) % 20 == 0:
            print(f"  tokenized {i + len(batch)}/{len(lines)} lines", file=sys.stderr)
    return counts


def take_prefix(lines, counts, target):
    """Smallest prefix of lines whose token sum first reaches `target`.
    Returns (indices, token_sum)."""
    s = 0
    out = []
    for i, c in enumerate(lines):
        out.append(i)
        s += counts[i]
        if s >= target:
            break
    return out, s


def take_suffix(lines, counts, target):
    """Smallest suffix of lines whose token sum first reaches `target`.
    Returns (indices in original order, token_sum)."""
    s = 0
    out = []
    for i in range(len(lines) - 1, -1, -1):
        out.append(i)
        s += counts[i]
        if s >= target:
            break
    out.sort()
    return out, s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_file", required=True)
    ap.add_argument("--tokenizer_dir", required=True)
    ap.add_argument("--out_a", required=True)
    ap.add_argument("--out_b", required=True)
    ap.add_argument("--target_tokens", type=int, default=19_000_000,
                    help="BPE tokens per chunk. 0 => two equal disjoint halves "
                         "(use all data, split at midpoint).")
    ap.add_argument("--mode", choices=["contiguous", "random"], default="contiguous")
    ap.add_argument("--shuffle_seed", type=int, default=2026,
                    help="Used only in --mode random.")
    ap.add_argument("--manifest", default=None,
                    help="Optional JSON file recording the split for provenance.")
    args = ap.parse_args()

    print(f"[split] loading tokenizer {args.tokenizer_dir}", file=sys.stderr)
    tok = AutoTokenizer.from_pretrained(args.tokenizer_dir)

    with open(args.train_file) as f:
        lines = f.readlines()
    print(f"[split] {len(lines)} lines (conversations)", file=sys.stderr)

    counts = count_line_tokens(tok, lines)
    total = sum(counts)
    print(f"[split] total BPE tokens: {total:,}", file=sys.stderr)

    if args.mode == "contiguous":
        if args.target_tokens == 0:
            # two equal disjoint halves at the token midpoint
            half = total // 2
            idx_a, sa = take_prefix(lines, counts, half)
            used_a = set(idx_a)
            idx_b = [i for i in range(len(lines)) if i not in used_a]
            sb = sum(counts[i] for i in idx_b)
        else:
            idx_a, sa = take_prefix(lines, counts, args.target_tokens)
            idx_b, sb = take_suffix(lines, counts, args.target_tokens)
            overlap = set(idx_a) & set(idx_b)
            if overlap:
                sys.exit(f"ERROR: chunks overlap ({len(overlap)} lines). "
                         f"target_tokens={args.target_tokens} too large for "
                         f"corpus of {total:,} tokens; max per disjoint chunk "
                         f"is ~{total // 2:,}.")
    else:  # random
        order = list(range(len(lines)))
        random.Random(args.shuffle_seed).shuffle(order)
        target = (total // 2) if args.target_tokens == 0 else args.target_tokens
        idx_a, sa = take_prefix([lines[i] for i in order],
                                [counts[i] for i in order], target)
        idx_a = sorted(order[i] for i in idx_a)
        remaining = [i for i in order if i not in set(idx_a)]
        idx_b, sb = take_prefix([lines[i] for i in remaining],
                                [counts[i] for i in remaining], target)
        idx_b = sorted(remaining[i] for i in idx_b)

    def write_chunk(path, idx):
        with open(path, "w") as f:
            for i in idx:
                f.write(lines[i])

    write_chunk(args.out_a, idx_a)
    write_chunk(args.out_b, idx_b)

    words_a = sum(len(lines[i].split()) for i in idx_a)
    words_b = sum(len(lines[i].split()) for i in idx_b)
    print(f"[split] chunk A: {len(idx_a):,} lines, {sa:,} BPE tokens, "
          f"~{words_a:,} words -> {args.out_a}", file=sys.stderr)
    print(f"[split] chunk B: {len(idx_b):,} lines, {sb:,} BPE tokens, "
          f"~{words_b:,} words -> {args.out_b}", file=sys.stderr)
    print(f"[split] disjoint: {len(set(idx_a) & set(idx_b)) == 0}", file=sys.stderr)

    if args.manifest:
        with open(args.manifest, "w") as f:
            json.dump({
                "mode": args.mode,
                "target_tokens": args.target_tokens,
                "shuffle_seed": args.shuffle_seed,
                "total_bpe_tokens": total,
                "chunk_a": {"lines": len(idx_a), "bpe_tokens": sa, "words": words_a,
                            "line_indices": idx_a},
                "chunk_b": {"lines": len(idx_b), "bpe_tokens": sb, "words": words_b,
                            "line_indices": idx_b},
            }, f)
        print(f"[split] wrote manifest {args.manifest}", file=sys.stderr)


if __name__ == "__main__":
    main()
