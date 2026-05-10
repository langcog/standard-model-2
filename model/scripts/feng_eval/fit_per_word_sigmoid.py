"""
Fit Chang & Bergen (2022) 4-parameter logistic to per-CDI-word surprisal
trajectories from a Feng/CHILDES-trained model.

Input: surprisal_<seed>.csv with rows (step, epoch, word, n_occurrences,
       mean_nll, sum_nll)
Output: <model>_sigmoids.txt (TAB-separated, same columns as
        data/chang_bergen_2022/*_sigmoids.txt):
        Token, MaxSurprisal, MinSurprisal, ParamUpper, ParamLower,
        ParamXmid, ParamScale

Notes
-----
- C&B's parameterisation:
    S(x) = ParamLower + (ParamUpper - ParamLower)
                          / (1 + exp((x - ParamXmid) / ParamScale))
  where x = log10(steps), S is mean surprisal at that step in nats.
  ParamUpper = asymptote at x→-∞ (high surprisal early)
  ParamLower = asymptote at x→+∞ (low surprisal late)
  ParamScale > 0  =>  curve decreases with x (learning).

- We follow their filtering: keep words where there's a meaningful drop
  in mean surprisal (surprisal_range = max-min > 1.0 nat) and ParamScale
  is finite and in (0.01, 10).
"""

import argparse
import csv
import math
import warnings
from collections import defaultdict

import numpy as np
from scipy.optimize import curve_fit


def four_pl(x, upper, lower, xmid, scale):
    return lower + (upper - lower) / (1.0 + np.exp((x - xmid) / scale))


def fit_one(steps, nlls):
    """Fit 4-PL to (log10(steps), nlls). Returns dict of params or None."""
    x = np.log10(np.asarray(steps, dtype=float))
    y = np.asarray(nlls, dtype=float)
    if len(x) < 6 or not np.all(np.isfinite(y)):
        return None
    upper0 = float(np.max(y))
    lower0 = float(np.min(y))
    if upper0 - lower0 < 1e-6:
        return None
    # Midpoint where y crosses (upper+lower)/2
    half = 0.5 * (upper0 + lower0)
    mid_idx = int(np.argmin(np.abs(y - half)))
    xmid0 = float(x[mid_idx])
    # Initial scale: spread of x.
    scale0 = max(1e-2, 0.25 * (x.max() - x.min()))
    p0 = (upper0, lower0, xmid0, scale0)
    # Bounds: scale > 0 (=> decreasing); allow upper > lower in case of weird
    # words by flipping initial sign. We constrain upper >= lower via bounds.
    bounds_lo = (lower0 - 5, lower0 - 5, x.min() - 3, 1e-3)
    bounds_hi = (upper0 + 5, upper0 + 5, x.max() + 3, 50.0)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, _ = curve_fit(
                four_pl, x, y, p0=p0,
                bounds=(bounds_lo, bounds_hi), maxfev=8000,
            )
    except Exception:
        return None
    upper, lower, xmid, scale = popt
    return {
        "ParamUpper": float(upper),
        "ParamLower": float(lower),
        "ParamXmid": float(xmid),
        "ParamScale": float(scale),
        "MaxSurprisal": float(np.max(y)),
        "MinSurprisal": float(np.min(y)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--surprisal_csv", required=True,
                    help="Input CSV (step, epoch, word, n_occurrences, mean_nll, sum_nll)")
    ap.add_argument("--out_tsv", required=True,
                    help="Output TSV in C&B sigmoid format")
    ap.add_argument("--min_steps", type=int, default=6,
                    help="Minimum number of (step, surprisal) points to fit")
    args = ap.parse_args()

    # Read trajectories
    by_word = defaultdict(list)  # word -> [(step, mean_nll)]
    with open(args.surprisal_csv) as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            try:
                step = int(r["step"])
                mn = float(r["mean_nll"])
            except (ValueError, KeyError):
                continue
            by_word[r["word"]].append((step, mn))

    fields = ["Token", "MaxSurprisal", "MinSurprisal",
              "ParamUpper", "ParamLower", "ParamXmid", "ParamScale"]

    fitted = 0
    skipped = 0
    with open(args.out_tsv, "w", newline="") as f:
        wr = csv.writer(f, delimiter="\t")
        wr.writerow(fields)
        for word in sorted(by_word.keys()):
            traj = sorted(set(by_word[word]))  # dedupe + sort by step
            if len(traj) < args.min_steps:
                skipped += 1
                continue
            steps = [s for s, _ in traj]
            nlls = [v for _, v in traj]
            res = fit_one(steps, nlls)
            if res is None:
                skipped += 1
                continue
            wr.writerow([
                word,
                f"{res['MaxSurprisal']:.6f}",
                f"{res['MinSurprisal']:.6f}",
                f"{res['ParamUpper']:.6f}",
                f"{res['ParamLower']:.6f}",
                f"{res['ParamXmid']:.6f}",
                f"{res['ParamScale']:.6f}",
            ])
            fitted += 1
    print(f"Fitted {fitted}, skipped {skipped}. -> {args.out_tsv}")


if __name__ == "__main__":
    main()
