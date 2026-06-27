"""Build a novel-word continuation corpus by NONCE-BY-SUBSTITUTION.

Take held-out CHILDES carrier UTTERANCES (split on the literal '\\n\\n' separator used in the
corpus; each utterance is '**SPK**: words ...'); for each (anchor, nonce) pair replace the real
anchor word (e.g. "dog") with a novel string (e.g. "blicket") in a controlled number of
utterances. The nonce then occupies the anchor's natural contexts at a frequency we set, and the
saved models have never seen it. Nonces are grouped into COHORTS introduced at STAGGERED points
in the (no-shuffle) stream, supporting both:
  (A) across warm-started regimes -- who acquires the nonces faster?
  (B) within a run -- are LATER cohorts acquired faster? (kappa / increasing-returns test)

The stream is re-joined into the corpus's native line format (utterances joined by ' \\n\\n ').
Outputs in --out_dir: continuation_train.txt, nonce_probe.jsonl, nonce_manifest.csv.
"""
import argparse, csv, json, os, random, re
from transformers import AutoTokenizer

# pool of common CHILDES nouns; we auto-pick the n_pairs most frequent that have enough utterances
ANCHOR_POOL = ["dog","ball","cup","shoe","hat","car","book","cat","milk","sock","duck","bear",
               "apple","chair","truck","bird","fish","box","cookie","baby","cow","hand","door","spoon",
               "bottle","banana","cracker","blanket","towel","spoon","table","window","flower","horse",
               "balloon","bubble","kitty","doggie","blocks","crayon","diaper","bike","boat","plane"]
NONCES = ["dax","blicket","wug","fep","gorp","zud","kiv","nale","sprock","glim","vask","murn",
          "bink","lerp","tron","quib","frell","snod","plon","draz","vimp","churg","yomp","tisp"]
SEP = " \\n\\n "   # literal separator as it appears in the corpus


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tokenizer_dir", required=True)
    ap.add_argument("--carrier_file", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--cohorts", type=int, default=4)
    ap.add_argument("--per_cohort", type=int, default=6)
    ap.add_argument("--train_freq", type=int, default=100)
    ap.add_argument("--probe_freq", type=int, default=20)
    ap.add_argument("--utts_per_line", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    os.makedirs(args.out_dir, exist_ok=True)

    tok = AutoTokenizer.from_pretrained(args.tokenizer_dir)
    assert tok.is_fast, "need a fast tokenizer for offset mapping"

    raw = open(args.carrier_file).read()
    utts = [u.strip() for u in raw.split("\\n\\n") if u.strip()]
    rng.shuffle(utts)
    print(f"[corpus] carrier utterances: {len(utts)}", flush=True)

    n_pairs = args.cohorts * args.per_cohort
    need = args.train_freq + args.probe_freq
    pats = {a: re.compile(r"\b" + re.escape(a) + r"\b", re.I) for a in dict.fromkeys(ANCHOR_POOL)}
    by_anchor = {a: [] for a in pats}
    used = set()
    for i, u in enumerate(utts):
        for a in pats:
            if pats[a].search(u):
                by_anchor[a].append(i)
    # pick the n_pairs commonest anchors with >= need utterances
    ok = sorted([a for a in pats if len(by_anchor[a]) >= need], key=lambda a: -len(by_anchor[a]))
    if len(ok) < n_pairs:
        raise SystemExit(f"only {len(ok)} anchors have >= {need} utts; lower train/probe_freq or n_pairs")
    anchors = ok[:n_pairs]
    print(f"[corpus] anchors: " + ", ".join(f"{a}({len(by_anchor[a])})" for a in anchors), flush=True)

    # claim exposure utterances (disjoint across anchors), mark used
    claim = {}
    for a in anchors:
        avail = [i for i in by_anchor[a] if i not in used][:need]
        for i in avail:
            used.add(i)
        claim[a] = avail

    probe_rows, manifest = [], []
    cohort_train = {c: [] for c in range(args.cohorts)}
    for k, (a, nz) in enumerate(zip(anchors, NONCES[:n_pairs])):
        c = k // args.per_cohort
        intro_frac = c / args.cohorts
        sub = [pats[a].sub(nz, utts[i]) for i in claim[a]]
        train_s, probe_s = sub[:args.train_freq], sub[args.train_freq:need]
        cohort_train[c].extend(train_s)
        npr = 0
        for s in probe_s:
            m = re.search(r"\b" + re.escape(nz) + r"\b", s)
            if not m:
                continue
            enc = tok(s, return_offsets_mapping=True)
            ids, offs = enc["input_ids"], enc["offset_mapping"]
            span = [j for j, (st, en) in enumerate(offs) if st < m.end() and en > m.start() and en > st]
            if not span:
                continue
            probe_rows.append({"nonce": nz, "cohort": c, "intro_frac": intro_frac,
                               "train_freq": args.train_freq, "ids": ids, "span": [span[0], span[-1] + 1]})
            npr += 1
        manifest.append((nz, a, c, intro_frac, args.train_freq, npr))

    # filler = unused utterances; assemble stream cohort-by-cohort (staggered intro)
    filler = [utts[i] for i in range(len(utts)) if i not in used]
    fil_per = len(filler) // args.cohorts
    stream = []
    for c in range(args.cohorts):
        block = list(cohort_train[c]) + filler[c * fil_per:(c + 1) * fil_per]
        rng.shuffle(block)
        stream.extend(block)
    # re-join into native line format (utterances joined by literal ' \n\n ')
    lines = [SEP.join(stream[i:i + args.utts_per_line]) for i in range(0, len(stream), args.utts_per_line)]
    with open(os.path.join(args.out_dir, "continuation_train.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(args.out_dir, "nonce_probe.jsonl"), "w") as f:
        for r in probe_rows:
            f.write(json.dumps(r) + "\n")
    with open(os.path.join(args.out_dir, "nonce_manifest.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["nonce","anchor","cohort","intro_frac","train_freq","probe_n"])
        w.writerows(manifest)
    toks = sum(len(s.split()) for s in stream) * 13 // 10
    print(f"[corpus] {n_pairs} nonces x {args.cohorts} cohorts; stream utts={len(stream)} "
          f"({len(lines)} lines); probe ctx={len(probe_rows)}; ~{toks} tokens "
          f"(~{toks//16384} steps @ bs16)", flush=True)


if __name__ == "__main__":
    main()
