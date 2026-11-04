"""What does CLASSIC's learning time carry beyond frequency?  Compute per-item distributional
features of the CDI words in the same CHILDES pool (utterance position, utterance length, context
diversity) and ask how much of log_t's increment over frequency survives once they are controlled.

usage: python 06_features.py WORKDIR [TAG] [AGE_LO AGE_HI]     (reads work/prep.pkl, work/runs_<TAG>.csv)
"""
import sys
import io
import contextlib
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

import analyze as A
from cvlib import cv_r2

W = sys.argv[1]
TAG = sys.argv[2] if len(sys.argv) > 2 else 'main'
LO = float(sys.argv[3]) if len(sys.argv) > 3 else 10
HI = float(sys.argv[4]) if len(sys.argv) > 4 else 36
D = A.D
MAXSPAN = 4

var2item = defaultdict(list)
for item, vs in D['targets'].items():
    for v in vs:
        var2item[v].append(item)

# ---- occurrences of each item in the age pool (contiguous whole-word spans, as in 03) ----
wstr = {}                                   # word tuple -> space-joined phonemes
n_tok = Counter()
n_final, n_alone, n_init = Counter(), Counter(), Counter()
len_sum = Counter()
left, right = defaultdict(Counter), defaultdict(Counter)
ages = D['ages']
for u, words in enumerate(D['utts']):
    if not (LO <= ages[u] < HI):
        continue
    nW = len(words)
    ss = [wstr.setdefault(w, ' '.join(w)) for w in words]
    for a in range(nW):
        s = ''
        for b in range(a, min(nW, a + MAXSPAN)):
            s = ss[b] if b == a else s + ' ' + ss[b]
            for it in var2item.get(s, ()):
                n_tok[it] += 1
                len_sum[it] += nW
                if b == nW - 1:
                    n_final[it] += 1
                    if a == 0:
                        n_alone[it] += 1
                if a == 0:
                    n_init[it] += 1
                left[it][words[a - 1] if a else 'BOS'] += 1
                right[it][words[b + 1] if b + 1 < nW else 'EOS'] += 1


def entropy(c):
    p = np.array(list(c.values()), float)
    p /= p.sum()
    return float(-(p * np.log(p)).sum())


rows = []
for it in D['targets']:
    n = n_tok[it]
    if n < 5:
        continue
    rows.append(dict(jj=it, n_tok=n, log_cnt=np.log(n),
                     p_final=n_final[it] / n, p_alone=n_alone[it] / n, p_init=n_init[it] / n,
                     log_len=np.log(len_sum[it] / n),
                     ctx_ent=(entropy(left[it]) + entropy(right[it])) / 2))
F = pd.DataFrame(rows)
F.to_csv(f'{W}/item_features_{TAG}.csv', index=False)

FEATS = ['p_final', 'p_alone', 'log_len', 'ctx_ent']
for ycol in ('delta198', 'delta681', 'emp_diff'):
    with contextlib.redirect_stdout(io.StringIO()):
        d = A.table(TAG, ycol)
    d = d.merge(F, on='jj').copy()
    d['lp2'] = d.log_p ** 2
    d['lc2'] = d.log_cnt ** 2
    y = d.delta_j_median.values
    cls = pd.get_dummies(d.cc, drop_first=True).astype(float).values.T.tolist()
    print(f'\n== {TAG}, y = {ycol}, n = {len(d)} (items with >= 5 tokens in ages {LO:g}-{HI:g} mo)')

    print('\nitem features vs frequency / difficulty (Spearman):')
    cols = ['log_p', 'log_cnt'] + FEATS + ['n_ph', 'log_t']
    print(pd.DataFrame({c: [A.corr(d[c], d.delta_j_median)] + [A.corr(d[c], d[k]) for k in ['log_p', 'log_t']]
                        for c in cols}, index=['delta', 'log_p', 'log_t']).T.round(2).to_string())

    print('\nwhat predicts log_t? partial r(feature, log_t | class, log_p):')
    print('  ' + '  '.join(f'{f} {A.pcorr(d[f], d.log_t, cls + [d.log_p]):+.2f}' for f in FEATS + ['n_ph']))
    r_base = A.r2(d.log_t.values, cls + [d.log_p, d.lp2])
    r_all = A.r2(d.log_t.values, cls + [d.log_p, d.lp2] + [d[f] for f in FEATS + ['n_ph']])
    print(f'  R2 of log_t on class+log_p+log_p^2 {r_base:.3f}; + features + n_ph {r_all:.3f} '
          f'(unexplained share {1 - r_all:.3f})')

    models = {
        'B   class x log_p + log_p^2':                           ['class*log_p', 'lp2'],
        'B + features':                                         ['class*log_p', 'lp2'] + FEATS,
        'B + features + n_ph':                                  ['class*log_p', 'lp2'] + FEATS + ['n_ph'],
        'B + log_t':                                            ['class*log_p', 'lp2', 'log_t'],
        'B + features + n_ph + log_t':                          ['class*log_p', 'lp2'] + FEATS + ['n_ph', 'log_t'],
        'Bc  class x log_cnt (CHILDES) + sq':                   ['class*log_cnt', 'lc2'],
        'Bc + features + n_ph':                                 ['class*log_cnt', 'lc2'] + FEATS + ['n_ph'],
        'Bc + features + n_ph + log_t':                         ['class*log_cnt', 'lc2'] + FEATS + ['n_ph', 'log_t'],
        'class x log_t + features + n_ph (no frequency)':       ['class*log_t'] + FEATS + ['n_ph'],
    }
    print('\n10-fold CV R^2 (30 repeats):')
    res = {}
    for name, terms in models.items():
        res[name] = cv_r2(d, terms)
        print(f'  {name:50s} {res[name]:6.3f}')
    inc = res['B + features + n_ph + log_t'] - res['B + features + n_ph']
    inc0 = res['B + log_t'] - res['B   class x log_p + log_p^2']
    print(f'  increment of log_t over B: {inc0:+.3f};  over B + features + n_ph: {inc:+.3f}')

    full = ['class*log_p', 'lp2'] + FEATS + ['n_ph', 'log_t']
    print('  leave-one-feature-out (from B + features + n_ph + log_t), change in CV R^2:')
    for f in FEATS + ['n_ph']:
        t = [x for x in full if x != f]
        print(f'    drop {f:8s} {cv_r2(d, t) - res["B + features + n_ph + log_t"]:+.3f}')

    ctrl = cls + [d.log_p] + [d[f] for f in FEATS + ['n_ph']]
    rng = np.random.default_rng(0)
    bs = []
    for _ in range(1000):
        s = d.iloc[rng.integers(0, len(d), len(d))]
        b = pd.get_dummies(s.cc, drop_first=True).astype(float).values.T.tolist()
        if len(b) != len(cls):
            continue
        bs.append(A.pcorr(s.log_t, s.delta_j_median.values, b + [s.log_p] + [s[f] for f in FEATS + ['n_ph']]))
    print(f'  partial r(log_t, difficulty | class, log_p, features, n_ph) = '
          f'{A.pcorr(d.log_t, y, ctrl):+.3f}  95% boot CI {np.percentile(bs, [2.5, 97.5]).round(3)}')
