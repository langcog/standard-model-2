# CHILDES flat-vs-decay LR experiment — results

**Full writeup: [`writeup.pdf`](writeup.pdf)** (5 pp, figures). This file = quick in-repo reference.

**Question.** Does the sandbox finding — a flat, low LR gives kid-like *increasing* returns (J>0) —
hold on real child-directed data (GPT-2 small, full 24M-word CHILDES, multi-epoch), and is it a
genuine learning-to-learn effect or an optimization artifact? Readouts on a shared step axis:
val_bpb, CDI-word surprisal→J (lexical), **Zorro** (grammar, primary; CHILDES-vocab), BLiMP (grammar,
secondary; vocab-mismatched). Arms: decay 1e-4 (standard) vs flat {1e-4, 3e-5, 1e-5} × 2 seeds ×
{1× = 20 ep, 4× = 60 ep}. 12 runs, 0 failures.

## Endpoints (mean over 2 seeds)
| budget | schedule (LR) | val_bpb | Zorro | BLiMP | J |
|---|---|---|---|---|---|
| 1× | decay 1e-4 (standard) | 2.064 | 0.694 | 0.576 | +0.42 |
| 1× | flat 1e-4 | 2.225 | 0.674 | 0.576 | +0.32 |
| 1× | flat 3e-5 | 2.030 | 0.701 | 0.577 | +1.02 |
| 1× | **flat 1e-5** | **2.017** | 0.690 | 0.589 | **+1.60** |
| 4× | decay 1e-4 (standard) | **2.673** ⚠ | 0.689 | 0.576 | **−0.60** |
| 4× | flat 1e-5 | 2.106 | 0.689 | 0.580 | +1.29 |

## Four findings
1. **J signature replicates.** Lower flat LR → higher J (1.60 / 1.02 / 0.3–0.4 at flat 1e-5 / 3e-5 / 1e-4≈decay). `fig_J_dose.png`
2. **flat-low WINS val loss on CHILDES** (2.017 < decay 2.064) — flips the web-text sandbox's +0.04–0.09 bpb *cost* into a benefit, because small-data multi-epoch → gentle low LR overfits less. `fig_valbpb_1x.png`
3. **4× = decisive.** Standard decay overfits catastrophically with more epochs (val_bpb 2.06→2.67, J flips +0.42→−0.60); flat-low is robust (2.02→2.11, J stays +1.3). `fig_overfit_4x.png`
4. **Grammar is DECOUPLED.** Zorro plateaus ~0.68–0.70 across *every* arm/schedule/budget; BLiMP ~0.58 (ceiling). The schedule reshapes loss + lexical returns but **not** grammatical competence. `fig_zorro.png`

## Novel-word acquisition (the decisive learning-to-learn test — DONE)
Warm-started all 12 saved checkpoints and continued under an identical protocol (held-out CHILDES carrier
+ 24 **novel words** via nonce-substitution, constant LR 1e-5) — see [`writeup.pdf`](writeup.pdf) §6,
`analyze_nonce.R`, `nonce_by_arm.csv`. Metric = nonce surprisal **floor** after ~100 exposures (lower =
learned better). **New-word mastery tracks the REGIME, not competence:** cor(floor, J) = **−0.80** vs
cor(floor, val_bpb) = +0.63; lm(floor ~ J + val_bpb): J coef −2.51 **p=0.026**, val_bpb p=0.51. Flat-low
(high-J) checkpoints learn new words to the lowest floors; the overfit decay-4× checkpoints (high val_bpb,
J<0) are the *worst* despite the most training. Within-run κ (later cohorts faster?) is null once forgetting
is controlled (expected: converged model re-seeing known data isn't developing). Figs:
`fig_nonce_regime_J.png` (the result), `fig_nonce_regime_valbpb.png` (control), `fig_nonce_kappa.png`.

## Verdict
**Lexical learning-to-learn, syntactically inert.** The flat-low-LR schedule enhances the **lexical-learning
faculty specifically**: it gives the kid-like increasing-returns J on known words, lowers val loss (anti-
overfit in this multi-epoch regime), survives extended training, *and* makes the model a genuinely faster
learner of **new** words — tracking the regime (J), not starting competence. But it leaves **grammar
untouched** (Zorro flat across all arms). So an LM *can* be pushed toward the kid-like word-learning
signature via a training knob — but it's a lexical-optimization faculty, decoupled from structural
competence; mechanistically narrower than (not identical to) the child's acceleration.

Data: ccn2 `/data2/mcfrank/ladder/flatlr_grammar/runs/`. Scripts in this dir; `endpoints.csv` has all 12 arms.
