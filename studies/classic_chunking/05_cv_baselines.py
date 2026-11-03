"""Fairer baselines: class-specific frequency slopes (v2's M4-style), cross-validated R2.
usage: python analyze2.py WORKDIR TAG"""
import sys
import io
import contextlib
import numpy as np
import pandas as pd
import analyze as A

tag = sys.argv[2] if len(sys.argv) > 2 else 'main'


def design(d, terms):
    cols = [np.ones(len(d))]
    cls = pd.get_dummies(d.cc).astype(float)
    for t in terms:
        if t == 'class':
            cols += [cls[c].values for c in cls.columns[1:]]
        elif t.startswith('class*'):
            v = d[t.split('*')[1]].values
            cols += [cls[c].values * v for c in cls.columns]          # class-specific slopes
            if 'class' not in terms:
                cols += [cls[c].values for c in cls.columns[1:]]
        else:
            cols.append(d[t].values)
    return np.column_stack(cols)


def cv_r2(d, terms, k=10, reps=30, seed=0):
    y = d.delta_j_median.values
    X = design(d, terms)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(reps):
        perm = rng.permutation(len(y))
        pred = np.empty(len(y))
        for f in np.array_split(perm, k):
            tr = np.setdiff1d(perm, f)
            b, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=1e-10)
            pred[f] = X[f] @ b
        vals.append(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    return np.mean(vals)


models = {
    'class only':                          ['class'],
    'class + log_p (common slope)':        ['class', 'log_p'],
    'class x log_p (class slopes)':        ['class*log_p'],
    'class x log_p + n_ph':                ['class*log_p', 'n_ph'],
    'class x log_p + log_p^2':             ['class*log_p', 'lp2'],
    'class x log_p + log_p^2 + n_ph':      ['class*log_p', 'lp2', 'n_ph'],
    'class x log_p + log_p^2 + log_t':     ['class*log_p', 'lp2', 'log_t'],
    'class x log_t (CLASSIC only)':        ['class*log_t'],
    'class + log_t':                       ['class', 'log_t'],
    'class x log_p + log_t':               ['class*log_p', 'log_t'],
    'class x log_p + class x log_t':       ['class*log_p', 'class*log_t'],
    'class x log_p + log_t + log_tok':     ['class*log_p', 'log_t', 'log_tok'],
}
for ycol in ('delta198', 'emp_diff'):
    with contextlib.redirect_stdout(io.StringIO()):
        d = A.table(tag, ycol)
    d = d.copy()
    d['lp2'] = d.log_p ** 2
    print(f'\n== {tag}, y = {ycol}, n = {len(d)}  (10-fold CV R^2, 30 repeats)')
    for name, terms in models.items():
        print(f'  {name:36s} {cv_r2(d, terms):6.3f}')
