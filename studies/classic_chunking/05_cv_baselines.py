"""Fairer baselines: class-specific frequency slopes (v2's M4-style), cross-validated R2.
usage: python 05_cv_baselines.py WORKDIR TAG"""
import sys
import io
import contextlib
import numpy as np
import pandas as pd
import analyze as A
from cvlib import design, cv_r2  # noqa: F401

tag = sys.argv[2] if len(sys.argv) > 2 else 'main'


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
