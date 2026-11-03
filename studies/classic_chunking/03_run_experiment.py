"""Run CLASSIC on age-ordered CHILDES streams; record, per CDI item and run, the
utterance at which the item first becomes a lexical/multiword chunk and the number
of tokens of the item heard before then.

usage: python run_experiment.py WORKDIR N_UTTS N_RUNS AGE_LO AGE_HI TAG [LR]
"""
import os
import sys
import pickle
import random
from bisect import bisect_right
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classic_fast import run  # noqa: E402

W, N, R, LO, HI, TAG = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), \
    float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
LR = float(sys.argv[7]) if len(sys.argv) > 7 else 0.10

D = pickle.load(open(f'{W}/prep.pkl', 'rb'))
ages = np.array(D['ages'])
pool = np.where((ages >= LO) & (ages < HI))[0]
targets = D['targets']
var2item = {}
for item, vs in targets.items():
    for v in vs:
        var2item.setdefault(v, []).append(item)
MAXSPAN = 4


def occurrences(stream):
    """item -> sorted list of utterance indices (1-based) containing it as a
    contiguous run of whole words (up to MAXSPAN words)."""
    occ = {i: [] for i in targets}
    for u, words in enumerate(stream, start=1):
        seen = set()
        nW = len(words)
        for a in range(nW):
            for b in range(a + 1, min(nW, a + MAXSPAN) + 1):
                s = ' '.join(' '.join(w) for w in words[a:b])
                if s in var2item:
                    for it in var2item[s]:
                        if it not in seen:
                            seen.add(it)
                            occ[it].append(u)
    return occ


def one(seed):
    rng = random.Random(10_000 + seed)
    idx = sorted(rng.sample(list(pool), N), key=lambda i: (ages[i], rng.random()))
    stream = [D['utts'][i] for i in idx]
    res = run(stream, targets, seed, D['phonemes'], lr=LR)
    occ = occurrences(stream)
    rows = []
    for item in targets:
        u = res['learned_at'].get(item)
        o = occ[item]
        tot = len(o)
        if u is None:
            rows.append((seed, item, 0, N, tot, tot))
        else:
            rows.append((seed, item, 1, u, bisect_right(o, u), tot))
    growth = [(seed,) + g for g in res['growth']]
    return rows, growth


if __name__ == '__main__':
    with Pool(min(8, os.cpu_count() or 2)) as p:
        out = p.map(one, range(R))
    rows = [r for o in out for r in o[0]]
    gr = [g for o in out for g in o[1]]
    pd.DataFrame(rows, columns=['seed', 'jj', 'learned', 'learn_utt', 'tokens_before',
                                'tokens_total']).to_csv(f'{W}/runs_{TAG}.csv', index=False)
    pd.DataFrame(gr, columns=['seed', 'utt', 'cdi_known', 'n_lex', 'n_multi', 'n_sub']
                 ).to_csv(f'{W}/growth_{TAG}.csv', index=False)
    print('done', TAG, 'pool', len(pool), 'runs', R)
