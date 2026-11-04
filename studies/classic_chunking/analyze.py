"""Compare CLASSIC per-word learning difficulty with v2's delta_j."""
import sys
import pickle
import numpy as np
import pandas as pd

W = sys.argv[1]
TAGS = sys.argv[2:] or ['main']
D = pickle.load(open(f'{W}/prep.pkl', 'rb'))
items = pd.read_csv(f'{W}/v2_items_delta.csv')
items = items[items.jj.isin(D['targets'])].copy()
items['n_ph'] = [len(D['targets'][j][0].split(' ')) for j in items.jj]
items['log_p'] = np.log(items.prob.clip(lower=1e-9))


def rank(x):
    return pd.Series(x).rank().to_numpy()


def corr(a, b, method='spearman'):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if method == 'spearman':
        a, b = rank(a), rank(b)
    return np.corrcoef(a, b)[0, 1]


def resid(y, X):
    X = np.column_stack([np.ones(len(y))] + [np.asarray(x, float) for x in X])
    beta, *_ = np.linalg.lstsq(X, np.asarray(y, float), rcond=None)
    return np.asarray(y, float) - X @ beta


def pcorr(a, b, ctrl):
    return corr(resid(a, ctrl), resid(b, ctrl), 'pearson')


def r2(y, X):
    e = resid(y, X)
    y = np.asarray(y, float)
    return 1 - (e ** 2).sum() / ((y - y.mean()) ** 2).sum()


def summarize(tag):
    r = pd.read_csv(f'{W}/runs_{tag}.csv')
    N = r.learn_utt.max()
    r['log_t'] = np.log(r.learn_utt)
    r['log_tok'] = np.log(r.tokens_before.clip(lower=1))
    g = r.groupby('jj').agg(frac=('learned', 'mean'), log_t=('log_t', 'mean'),
                            log_tok=('log_tok', 'mean'),
                            tok_total=('tokens_total', 'mean')).reset_index()
    return g, N


def table(tag, ycol='delta198'):
    g, N = summarize(tag)
    d = items.merge(g, on='jj')
    d = d[(d.tok_total > 0) & d[ycol].notna()].copy()   # items never heard say nothing
    print(f'\n=== {tag}  y={ycol} (N={N} utts; {len(d)} items; mean frac learned {d.frac.mean():.2f})')
    d['delta_j_median'] = d[ycol]
    y = d.delta_j_median.values
    rows = []
    for name in ['log_p', 'n_ph', 'log_t', 'log_tok', 'frac']:
        rows.append((name, corr(d[name], y), corr(d[name], y, 'pearson')))
    print(pd.DataFrame(rows, columns=['predictor', 'spearman_w_delta', 'pearson_w_delta']).round(3).to_string(index=False))
    print('partials vs delta_j (Pearson of residuals):')
    for name in ['log_t', 'log_tok']:
        print(f'  {name} | log_p: {pcorr(d[name], y, [d.log_p]):+.3f}   '
              f'| log_p, n_ph: {pcorr(d[name], y, [d.log_p, d.n_ph]):+.3f}')
    base = r2(y, [d.log_p])
    base2 = r2(y, [d.log_p, d.n_ph])
    print(f'R2 delta ~ log_p {base:.3f}; + n_ph {base2:.3f}; + log_t {r2(y, [d.log_p, d.log_t]):.3f}; '
          f'+ log_tok {r2(y, [d.log_p, d.log_tok]):.3f}; + n_ph + log_t {r2(y, [d.log_p, d.n_ph, d.log_t]):.3f}; '
          f'log_t alone {r2(y, [d.log_t]):.3f}')
    print('by lexical class (spearman delta vs log_t | n):')
    for c, s in d.groupby('cc'):
        if len(s) > 15:
            print(f'  cc={c}: {corr(s.log_t, s.delta_j_median):+.3f} (n={len(s)}); log_p {corr(s.log_p, s.delta_j_median):+.3f}')
    return d


def fe(d, nboot=2000, seed=0):
    """Class fixed effects (v2 has a lexical-class hierarchy on delta_j): does CLASSIC
    learning time add to log_p within class?"""
    y = d.delta_j_median.values
    cls = pd.get_dummies(d.cc, drop_first=True).astype(float).values.T.tolist()
    base = cls
    out = {
        'class': r2(y, base), '+log_p': r2(y, base + [d.log_p]),
        '+log_p+n_ph': r2(y, base + [d.log_p, d.n_ph]),
        '+log_t': r2(y, base + [d.log_t]),
        '+log_p+log_t': r2(y, base + [d.log_p, d.log_t]),
        '+log_p+n_ph+log_t': r2(y, base + [d.log_p, d.n_ph, d.log_t]),
    }
    print('class-FE R2:', {k: round(v, 3) for k, v in out.items()})
    rng = np.random.default_rng(seed)
    pc, dr = [], []
    n = len(d)
    for _ in range(nboot):
        i = rng.integers(0, n, n)
        s = d.iloc[i]
        yy = s.delta_j_median.values
        b = pd.get_dummies(s.cc, drop_first=True).astype(float).values.T.tolist()
        if len(b) != len(cls):
            continue
        pc.append(pcorr(s.log_t, yy, b + [s.log_p]))
        dr.append(r2(yy, b + [s.log_p, s.log_t]) - r2(yy, b + [s.log_p]))
    q = lambda v: np.percentile(v, [2.5, 97.5]).round(3)
    print(f'partial r(log_t, delta | class, log_p) = {pcorr(d.log_t, y, cls + [d.log_p]):+.3f} '
          f'95% boot CI {q(pc)};  dR2 over class+log_p = {out["+log_p+log_t"] - out["+log_p"]:.3f} CI {q(dr)}')


if __name__ == '__main__':
    for t in TAGS:
        for ycol in ('delta198', 'delta681', 'emp_diff'):
            d = table(t, ycol)
            fe(d)
            if t == 'main' and ycol == 'emp_diff':
                d.to_csv(f'{W}/item_table_main.csv', index=False)
