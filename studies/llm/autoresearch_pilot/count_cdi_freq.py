#!/usr/bin/env python3
"""Count ClimbMix training-corpus frequency of each single-token CDI word
(uncapped), for use as the frequency control in the onset->kappa analysis.
Run from the autoresearch repo dir: python count_cdi_freq.py cdi_words.txt freq.csv"""
import csv, glob, os, sys
import pyarrow.parquet as pq
from prepare import Tokenizer, DATA_DIR, VAL_FILENAME

tok = Tokenizer.from_directory(); enc = tok.enc
words = [w.strip() for w in open(sys.argv[1]) if w.strip()]
wid = {}
for w in words:
    ids = enc.encode_ordinary(" " + w)
    if len(ids) == 1:
        wid[ids[0]] = w
counts = {w: 0 for w in wid.values()}
total = 0
for f in sorted(glob.glob(os.path.join(DATA_DIR, "*.parquet"))):
    if f.endswith(VAL_FILENAME):
        continue  # train shards only (the model's actual exposure)
    for b in pq.ParquetFile(f).iter_batches(columns=["text"], batch_size=1000):
        for ids in enc.encode_ordinary_batch(b.column("text").to_pylist(), num_threads=8):
            total += len(ids)
            for t in ids:
                w = wid.get(t)
                if w is not None:
                    counts[w] += 1
with open(sys.argv[2], "w", newline="") as o:
    wr = csv.writer(o); wr.writerow(["word", "count", "total_tokens"])
    for w, c in counts.items():
        wr.writerow([w, c, total])
print(f"counted {len(counts)} words over {total} train tokens", file=sys.stderr)
