# glmer ladder — the simple-measurement-model benchmark (2026-06-15)

**Purpose (MCF):** fit a 3-model glmer ladder with *simple* per-child measurement
summaries as an "unbiased" reference. The Bayesian io-proc model should *improve
on* this ladder, not disagree with it. Built from `fits/joint_io_proc_mm_subset_data.rds`
(N=768,599 CDI obs, I=377 kids, J=681 items, S=6 studies). Scripts:
`/tmp/glmer_ladder.R` (common) + `/tmp/glmer_ladder_standalone.R` (full samples).

## Channel summaries (per child, within-study z-scored)
- **input_z**: mean observed log-input per kid (from `log_input_obs`), z within study.
- **proc_z**: RT measurement model `lmer(lwl_log_rt ~ log_age + (1|dataset))`, take each
  kid's **mean residual**, flip sign, z within study → higher = **faster** (better processing).
- cor(input_z, proc_z) = **0.242** on the 142 both-channel kids (weak → dissociation is real).

## Common-sample ladder (same 142 both-channel kids; AM2018+FMW2013+SEEDLingS, 326,754 obs)
Model: `produces ~ <chan>*log_age + (log_age|child) + (1|item) + (1|dataset)`, a0=21.

| term | M1 input | M2 proc | **M3 both** |
|---|--:|--:|--:|
| input_z (intercept) | 0.285 (.12) | — | 0.13 (.49) |
| input_z:log_age (slope) | **0.918 (.009)** | — | **0.845 (.018)** |
| proc_z (intercept) | — | **0.662 (<.001)** | **0.625 (<.001)** |
| proc_z:log_age (slope) | — | 0.474 (.18) | 0.25 (.50) |

child random-slope SD: M1 2.14 → M3 2.05 (channels don't compete).

### ⭐ DOUBLE DISSOCIATION (robust to joint fit + low collinearity)
- **Input → ACCELERATION** (slope 0.85–0.92, stable; intercept ~0, n.s.).
- **Processing → LEVEL** (intercept 0.63, p<.001; slope n.s.).
Each survives controlling for the other.

## Standalone anchors (full own samples)
- **Input, 193 kids:** input_z=0.351 (.021), input_z:log_age=**0.903 (.003)**.
  → input→slope reproduces on max data; the sig intercept (0.35) is **processing
    leakage** (collapses to 0.13 n.s. once proc_z enters, M3) — NOT a real input-level effect.
- **Processing, 326 kids:** _pending (689k-obs fit)_

## Implication for the SM2 io-proc spec
The SM2 gets **processing→level right** (`xi += beta_xi*rt0`, ladder β_xi≈−1.6, big).
For **input** it does the opposite of what the data want:
1. **Hard-wires input into the LEVEL** via the io identity `xi = mu_r + 1*d_i`
   (forces σ_r≈0.36 of level variance onto input) — but the controlled input→level is **~0.13, n.s.**
2. **Under-detects input→slope**: the single N(0,1)-regularized, weakly-identified `gamma_in`
   collapses to γ_in·σ_r=0.26, vs the unbiased benchmark **0.85**. See [[mm_run_state]]
   (D'0–D'3 ladder: γ_in flat at ~0.73 across rungs → NOT processing competition).

**Takeaway:** the io premise (ability = input × efficiency ⇒ input lives in the level)
is the assumption the data push back on. Input looks like an **acceleration** term, not a
level term. The Bayesian model should be re-specified to let input load on the slope
(and not be forced onto the level at coeff 1), so it reproduces this benchmark.
