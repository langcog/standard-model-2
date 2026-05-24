# Retired R scripts

Archived 2026-05-23. These have been superseded by newer scripts that
generate the slide-deck figures.

| File | Replaced by |
|---|---|
| `quantile_demo.R` | `quantile_demo_4panel_I200.R` (architecture build at I=200) + `quantile_demo_mbest_EN_NO_wordbank.R` (EN+NO side-by-side at I=500). The original was a 6-panel demo including the `s_i` variants we've now excised. |
| `theta_spaghetti.R` | `theta_spaghetti_4panel_I200.R`. Same reason — original was 5-panel with the `s_i` panel. |

Both replacements use the post-excision M_best regime (`α + ζ + δ`,
no `s`, no `s_i`) and the GAMLSS BEINF empirical quantile smoother
implemented in `model/R/empirical_xsec_helper.R`.
