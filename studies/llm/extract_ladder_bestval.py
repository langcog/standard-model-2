# Build ladder_bestval_finer.csv: per (seed,budget), each word's surprisal at the
# converged (max-step = best-val, via load_best_model_at_end) checkpoint.
# Combines runs_grid (original-seed grid budgets) + runs_finer (inserted rungs + new seeds).
import glob, os, csv, re

def parse_budget(label):  # "0.75M" -> 750000
    return int(round(float(label[:-1]) * 1_000_000))

rows, words_of = {}, {}
for d in ["runs_grid", "runs_finer"]:        # runs_finer wins on any collision
    for f in sorted(glob.glob(os.path.join(d, "surprisal_seed*_*.csv"))):
        m = re.match(r"surprisal_seed(\d+)_(.+)\.csv$", os.path.basename(f))
        if not m:
            continue
        seed, blabel = int(m.group(1)), m.group(2)
        data = list(csv.DictReader(open(f)))
        if not data:
            continue
        mx = max(int(x["step"]) for x in data)
        rows[(seed, blabel)] = {x["word"]: float(x["mean_nll"]) for x in data if int(x["step"]) == mx}
        words_of[(seed, blabel)] = parse_budget(blabel)

with open("ladder_bestval_finer.csv", "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["seed", "rung", "words", "word", "surprisal"])
    for k in sorted(rows, key=lambda k: (k[0], words_of[k])):
        for word in sorted(rows[k]):
            w.writerow([k[0], k[1], words_of[k], word, "%.6f" % rows[k][word]])

seeds = sorted(set(s for s, _ in rows))
print("seeds:", seeds)
print("(seed,budget) pairs:", len(rows), "(expect 10x18=180)")
print("budgets/seed:", dict((s, sum(1 for (ss, _) in rows if ss == s)) for s in seeds))
print("total word-rows:", sum(len(v) for v in rows.values()))
