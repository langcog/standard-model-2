# Retired Stan models

Archived 2026-05-23. Not used by the current slide deck or the M_best
headline pipeline. Kept for git history / reference only.

| File | Why retired |
|---|---|
| `log_irt_long_si_signed.stan` | Per-child trajectory phase `s_i` excised from the headline model (see `experiments.md` §23). The signed-normal `s_i` variant developed multi-mode posteriors at I=500. |
| `log_irt_long_si_signed_annotated.stan` | Reading copy of the above with extra comments. Same retire reason. |
| `log_irt_long_si_corr.stan` | LKJ-correlated `(ξ, ζ, s_i)` variant. Same retire reason. |
| `log_irt_long_io_accel.stan` | IO + acceleration extension; not in the current deck. |
| `log_irt_long_io_accel_si_signed.stan` | + `s_i` extension. Both retire reasons apply. |
| `log_irt_long_io_comp_si_signed.stan` | IO + comp + signed `s_i`. The non-`s_i` variants (`log_irt_long_io_comp.stan`) are still active. |
| `log_irt_long_proc_si_signed.stan` | Peekbank/proc + signed `s_i`. The non-`s_i` proc variant is still active. |
| `log_irt_long_lmm.stan` | Linear-mixed-model comparison; cited briefly in earlier deck versions but removed in the current one. |

Current production Stan files (in `../`):

- `log_irt.stan` — cross-sectional baseline.
- `log_irt_io.stan` — cross-sectional + observed input rate.
- `log_irt_long.stan` — **the M_best headline**: `η = ξ + log H + (1+δ+ζ)·log(age/a_0) − δ_j`.
- `log_irt_long_io.stan` — longitudinal + observed input.
- `log_irt_long_io_comp.stan` — + comprehension channel (used by SEEDLingS comp/std fits).
- `log_irt_long_proc.stan` — + Peekbank LWL RT channel.
