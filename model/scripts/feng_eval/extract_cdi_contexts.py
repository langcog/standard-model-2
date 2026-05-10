"""
Pre-extract CDI-word contexts from CHILDES validation set.

For each of the 611 CDI words from Chang & Bergen (2022) that is a SINGLE TOKEN
in the GPT-2 CHILDES tokenizer (Feng et al. 2024/2026), find up to N_PER_WORD
occurrences in CHILDES_val_ordered.txt. For each occurrence, record:
    word, token_id, left_context_ids (truncated to last (max_ctx-1) tokens), position_in_left_context_plus_1

The resulting JSONL is the eval set used by the training-time callback.

Usage:
    python extract_cdi_contexts.py \
        --tokenizer <tokenizer_dir> \
        --val_file <CHILDES_val_ordered.txt> \
        --cdi_words <cdi_words.txt> \
        --out_jsonl <out.jsonl> \
        --max_per_word 200 --max_ctx 1024

Output JSONL row schema:
    {"word": "dog", "token_id": 1234, "ctx": [..int..], "pos": 17}
where ctx[pos] == token_id and (ctx[:pos]) is the left context.
"""

import argparse
import json
import os
import sys
from collections import defaultdict

from transformers import AutoTokenizer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tokenizer", required=True)
    ap.add_argument("--val_file", required=True)
    ap.add_argument("--cdi_words", required=True)
    ap.add_argument("--out_jsonl", required=True)
    ap.add_argument("--out_coverage", default=None,
                    help="optional CSV: word,is_single_token,token_id,n_occurrences")
    ap.add_argument("--max_per_word", type=int, default=200)
    ap.add_argument("--max_ctx", type=int, default=1024)
    args = ap.parse_args()

    tok = AutoTokenizer.from_pretrained(args.tokenizer)

    with open(args.cdi_words) as f:
        words = [w.strip() for w in f if w.strip()]
    print(f"Loaded {len(words)} CDI words", file=sys.stderr)

    # Determine single-token coverage. GPT-2 BPE: words preceded by a space
    # tokenize differently than at sequence start. We use the leading-space form
    # since most occurrences are mid-sentence. Record both for diagnostics.
    coverage = []
    word_to_id = {}
    for w in words:
        ids_with_space = tok.encode(" " + w, add_special_tokens=False)
        ids_no_space = tok.encode(w, add_special_tokens=False)
        single = (len(ids_with_space) == 1)
        token_id = ids_with_space[0] if single else None
        coverage.append({
            "word": w,
            "single_token": single,
            "token_id_with_space": ids_with_space[0] if len(ids_with_space) == 1 else None,
            "token_id_no_space": ids_no_space[0] if len(ids_no_space) == 1 else None,
            "n_pieces_with_space": len(ids_with_space),
            "n_pieces_no_space": len(ids_no_space),
        })
        if single:
            word_to_id[w] = token_id

    n_single = sum(1 for c in coverage if c["single_token"])
    print(f"Single-token (with leading space): {n_single}/{len(words)}", file=sys.stderr)

    id_to_word = {tid: w for w, tid in word_to_id.items()}
    counts = defaultdict(int)

    out_f = open(args.out_jsonl, "w")
    total_written = 0
    total_lines = 0

    with open(args.val_file) as f:
        for line in f:
            total_lines += 1
            if total_lines % 200 == 0:
                print(f"  line {total_lines}, written {total_written}", file=sys.stderr)
            line = line.rstrip("\n")
            if not line:
                continue
            ids = tok.encode(line, add_special_tokens=False)
            for pos, tid in enumerate(ids):
                if tid not in id_to_word:
                    continue
                w = id_to_word[tid]
                if counts[w] >= args.max_per_word:
                    continue
                # left context: tokens [0, pos]; we keep up to max_ctx-1 to leave
                # room for the target token. The target is at the end of ctx.
                lo = max(0, pos + 1 - args.max_ctx)
                ctx = ids[lo:pos + 1]
                tgt_pos = len(ctx) - 1  # index of the target within ctx
                rec = {"word": w, "token_id": tid, "ctx": ctx, "pos": tgt_pos}
                out_f.write(json.dumps(rec) + "\n")
                counts[w] += 1
                total_written += 1
            # quick early-exit if every word is at the cap (won't generally trigger)
            if total_lines % 500 == 0:
                if all(counts[w] >= args.max_per_word for w in word_to_id):
                    print("All words at cap, stopping early.", file=sys.stderr)
                    break

    out_f.close()

    # Coverage report
    for c in coverage:
        c["n_occurrences"] = counts.get(c["word"], 0)
    if args.out_coverage:
        import csv
        with open(args.out_coverage, "w", newline="") as cf:
            w = csv.DictWriter(cf, fieldnames=list(coverage[0].keys()))
            w.writeheader()
            for row in coverage:
                w.writerow(row)
        print(f"Wrote coverage to {args.out_coverage}", file=sys.stderr)

    # Summary
    occs = [counts.get(w, 0) for w in word_to_id]
    have_any = sum(1 for n in occs if n > 0)
    have_50 = sum(1 for n in occs if n >= 50)
    have_full = sum(1 for n in occs if n >= args.max_per_word)
    print(
        f"Wrote {total_written} contexts to {args.out_jsonl}. "
        f"single-token CDI words: {len(word_to_id)}, with ≥1 occ: {have_any}, "
        f"with ≥50: {have_50}, with full ({args.max_per_word}): {have_full}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
