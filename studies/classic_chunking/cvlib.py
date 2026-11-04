"""Shared design-matrix / cross-validation helpers (class-specific slopes via `class*x`)."""
import numpy as np
import pandas as pd


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
