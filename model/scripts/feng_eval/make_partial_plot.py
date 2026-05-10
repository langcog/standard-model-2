"""Quick density plot comparing partial Feng GPT-2-CHILDES slopes
to Chang & Bergen BookCorpus LM slopes (kids overlay drawn as a
synthetic Gaussian based on the published kappa posterior summary;
final plot uses the real posterior draws).

Run locally; matplotlib only.
"""
import argparse
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def slopes_from_sigmoid_tsv(path, scale_lo=0.01, scale_hi=10.0, range_min=1.0):
    slopes = []
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for r in rdr:
            try:
                ps = float(r["ParamScale"])
                pu = float(r["ParamUpper"])
                pl = float(r["ParamLower"])
            except (KeyError, ValueError):
                continue
            if not (scale_lo < ps < scale_hi):
                continue
            if (pu - pl) <= range_min:
                continue
            slopes.append(0.434 / ps)
    return np.array(slopes)


def slopes_from_cb(cb_dir, name):
    return slopes_from_sigmoid_tsv(os.path.join(cb_dir, f"{name}_sigmoids.txt"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--feng_tsvs", nargs="+", required=True,
                    help="Sigmoid TSVs from fit_per_word_sigmoid.py (one per seed).")
    ap.add_argument("--feng_labels", nargs="+", required=True,
                    help="Display labels for each Feng TSV.")
    ap.add_argument("--cb_dir", default="data/chang_bergen_2022")
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--kid_kappa_pop", type=float, default=10.3,
                    help="Published English M_best kappa_pop = 1+delta")
    ap.add_argument("--kid_sigma_zeta", type=float, default=3.47,
                    help="Published English M_best sigma_zeta")
    ap.add_argument("--title_suffix", default="(partial training)")
    args = ap.parse_args()
    assert len(args.feng_tsvs) == len(args.feng_labels)

    # ---- Kid kappa distribution (synthetic from published summary) ----
    rng = np.random.default_rng(2026)
    kid_slopes = args.kid_kappa_pop + rng.normal(0, args.kid_sigma_zeta, size=5000)

    cb_slopes = {
        "BERT (BookCorpus)":   slopes_from_cb(args.cb_dir, "bert"),
        "GPT-2 (BookCorpus)":  slopes_from_cb(args.cb_dir, "gpt2"),
        "BiLSTM (BookCorpus)": slopes_from_cb(args.cb_dir, "bilstm"),
        "LSTM (BookCorpus)":   slopes_from_cb(args.cb_dir, "lstm"),
    }

    feng_slopes = {}
    for tsv, label in zip(args.feng_tsvs, args.feng_labels):
        feng_slopes[label] = slopes_from_sigmoid_tsv(tsv)

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(9, 4.5))
    x = np.linspace(0, 22, 800)

    def kde(arr, bw=0.6):
        if len(arr) < 5:
            return np.zeros_like(x)
        arr = np.asarray(arr)
        return np.mean(
            np.exp(-0.5 * ((x[:, None] - arr[None, :]) / bw) ** 2),
            axis=1,
        ) / (bw * math.sqrt(2 * math.pi))

    # Kids
    ax.fill_between(x, 0, kde(kid_slopes, bw=1.0), color="#c41e37",
                    alpha=0.35, label=f"Kids (English M_best, n={len(kid_slopes)})")
    # C&B (light gray bundle of 4)
    for name, arr in cb_slopes.items():
        ax.fill_between(x, 0, kde(arr), color="#33a02c", alpha=0.18)
    ax.plot([], [], color="#33a02c", alpha=0.6, lw=4,
            label=f"C&B BookCorpus LMs (n=4 architectures)")

    # Feng partial — color each seed
    seed_colors = ["#1f78b4", "#6a3d9a", "#ff7f00"]
    for i, (label, arr) in enumerate(feng_slopes.items()):
        ax.fill_between(x, 0, kde(arr, bw=0.5), color=seed_colors[i % 3],
                        alpha=0.35,
                        label=f"{label} (n={len(arr)})")

    ax.axvline(1.0, color="grey", linestyle="--", lw=0.6)
    ax.text(1.05, ax.get_ylim()[1] * 0.9, "κ = 1\n(unit accumulator)",
            fontsize=8, color="grey")

    ax.set_xlabel("Slope on log(experience): logit per natural-log unit")
    ax.set_ylabel("Density")
    ax.set_xlim(0, 22)
    ax.set_title(
        f"Per-instance scaling slopes: children vs. CHILDES-trained GPT-2 vs. BookCorpus LMs {args.title_suffix}",
        fontsize=10, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    plt.tight_layout()
    fig.savefig(args.out_png, dpi=180)
    print(f"Wrote {args.out_png}")


if __name__ == "__main__":
    main()
