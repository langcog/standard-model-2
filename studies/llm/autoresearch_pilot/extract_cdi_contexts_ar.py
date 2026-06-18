#!/usr/bin/env python3
"""
Extract fixed-length CDI-word contexts from the autoresearch ClimbMix validation
shard, for the per-word surprisal probe (cdi_probe.py).

Restricts to CDI words that are a SINGLE token (space-prefixed " word") under the
autoresearch 8192-BPE tokenizer -- the only words with a well-defined single
target token to score (cf. Chang & Bergen 2022). For each such word we collect up
to --max-per-word occurrences in running validation text, each stored as a fixed
window of the K tokens ENDING in the word token.

Run inside the autoresearch repo (needs `from prepare import ...`):
    uv run extract_cdi_contexts_ar.py --cdi-words cdi_words.txt --out cdi_contexts.jsonl
"""
import argparse, json, os, sys
import pyarrow.parquet as pq
from prepare import Tokenizer, DATA_DIR, VAL_FILENAME


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cdi-words", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ctx-len", type=int, default=48,
                    help="fixed context length K in tokens (including the target word)")
    ap.add_argument("--max-per-word", type=int, default=80)
    ap.add_argument("--max-docs", type=int, default=300000)
    args = ap.parse_args()

    tok = Tokenizer.from_directory()
    enc = tok.enc

    with open(args.cdi_words) as f:
        words = [w.strip() for w in f if w.strip()]

    # single-token (space-prefixed) subset -> target token id
    word_to_id, id_to_word = {}, {}
    for w in words:
        ids = enc.encode_ordinary(" " + w)
        if len(ids) == 1:
            word_to_id[w] = ids[0]
            id_to_word[ids[0]] = w
    print(f"vocab={enc.n_vocab}  CDI words={len(words)}  single-token={len(word_to_id)}",
          file=sys.stderr)

    K = args.ctx_len
    target_ids = set(id_to_word)
    counts = {w: 0 for w in word_to_id}
    val_path = os.path.join(DATA_DIR, VAL_FILENAME)
    pf = pq.ParquetFile(val_path)

    n_ctx, docs_seen, done = 0, 0, False
    with open(args.out, "w") as out:
        for batch in pf.iter_batches(columns=["text"], batch_size=1000):
            texts = batch.column("text").to_pylist()
            id_lists = enc.encode_ordinary_batch(texts, num_threads=8)
            for ids in id_lists:
                docs_seen += 1
                for i, t in enumerate(ids):
                    if t in target_ids and i >= K - 1:
                        w = id_to_word[t]
                        if counts[w] >= args.max_per_word:
                            continue
                        out.write(json.dumps({"word": w, "ids": ids[i - (K - 1): i + 1]}) + "\n")
                        counts[w] += 1
                        n_ctx += 1
                if docs_seen >= args.max_docs:
                    done = True
                    break
            if done or all(c >= args.max_per_word for c in counts.values()):
                break

    have = sum(1 for c in counts.values() if c > 0)
    full = sum(1 for c in counts.values() if c >= args.max_per_word)
    print(f"contexts={n_ctx}  single-token words with >=1 occ={have}  at-cap={full}  "
          f"docs_seen={docs_seen}", file=sys.stderr)


if __name__ == "__main__":
    main()
