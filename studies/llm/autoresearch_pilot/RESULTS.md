# Acceleration campaign — can a small LM show *increasing* returns to experience?

**Question.** Children show *increasing* returns to experience over development — harder
words are acquired faster later (returns-to-experience exponent κ ≈ 10 over development).
Language models, scaled over data, show the opposite: *diminishing* returns (the
scaling-law default; our CHILDES ladder gave effective κ falling 1.25→0.47). This campaign
asks the hypothesis-generation question: **is there any single-GPU training intervention
that flips a small from-scratch LM (nanochat-style GPT-2, ~50M params, ClimbMix, ~5-min run
on one A40) from diminishing toward kid-like *increasing* returns?**

**Metric `J`.** Per-CDI-word mean surprisal is probed at ~36 log-spaced steps during one run
(`cdi_probe.py`). `J = r_late − r_early` = (slope of −mean_surprisal vs log10(experience) in
the late half) − (same in the early half). `J < 0` = diminishing returns (gains slow on the
log-experience axis); `J ≈ 0` = unit accumulator; `J > 0` = **increasing returns** (the
kid-like acceleration signature). Scored by `returns_score.R`.

**Artifact controls (the crux).** A model that barely learns early and then "catches up"
produces `J > 0` *mechanically* without any genuine acceleration. So `J > 0` counts as a
**real hit only if** (i) `r_early ≥ ~0.4` (early learning stayed healthy) **and** (ii)
`final_surp ≤ ~6.1` (competence near the baseline 5.88, i.e. it actually learned the words).
This control is what caught grokfast: at high amplification J shot up to +2.2, but `r_early`
collapsed to −0.24 — the model *got worse* early, so the "acceleration" was a slow-start
artifact, not learning-to-learn.

Baseline (stock decayed-LR single pass, `baseline_decay_p1`): **J = −0.37, r_early = 0.98,
final_surp = 5.88** — diminishing, as expected.

79 variants across 12 families. Plots: [J vs LR](fig_J_vs_LR.png) ·
[J vs competence (artifact diagnostic)](fig_J_vs_finalsurp.png) ·
[best J per family](fig_bestJ_by_family.png).

---

## VERDICT

**One robust lever flips the sign: a flat, low learning rate.** It is the only family that
produces `J > 0` while passing *both* artifact controls, and it does so as a clean,
monotonic dose-response.

- **Genuine frontier = flat LR ×0.07–0.1**: J = +1.25 / +1.11 at near-baseline competence
  (final_surp 5.99 / 5.92) with healthy early learning (r_early 0.73 / 0.77).
- **Champion = `lr_0p1_init0p3`** (flat LR ×0.1 + small-init ×0.3): **J = +1.28, r_early =
  0.77, final_surp = 5.76** — the *best competence of any run* (below baseline) at the *highest
  genuine J*. Low-LR and small-init **stack**.
- **Competence floor**: below ×0.07 (×0.05, ×0.03) J keeps climbing (+1.29, +1.37) but
  competence breaks (final_surp 6.29, 6.35 > 6.1; r_early → 0.49) — the model undertrains in
  one pass, so that extra J is the slow-start artifact, correctly excluded by the control.

**Mechanism (from the 2×2 + dissociations): the effect needs the *conjunction* of a FLAT
(non-decaying) schedule AND a LOW GLOBAL LR magnitude. Neither alone suffices.**

| | **low magnitude** | **full magnitude** |
|---|---|---|
| **flat** | `lr_0p1` **+1.11**, `lr_0p25` **+0.46** ✅ increasing | `constlr` −0.47, `warmdown_0` −0.31 ✗ |
| **decay** | `decaylow_0p1` −0.25, `decaylow_0p25` −0.40 ✗ | `baseline` −0.37 ✗ |

- Decay-on kills it at *any* magnitude: baseline (−0.37), `matlr_0p01` (low matrix-LR but decay
  on, −0.51), and the dedicated `decaylow` cell (all 4 LRs ×0.25/×0.1 with decay on; −0.40/−0.25
  at baseline competence). The decay schedule drives LR→0 late, which kills late-experience gains.
- Flat at *full* magnitude is also insufficient (`constlr` −0.47, `warmdown_0` −0.31).
- Only flat **and** low together gives increasing returns.

**Everything else is a null or actively antagonistic** (see tables): schedule *shape*
(warmup / warmdown / WSD / final-LR), matrix-LR alone, larger batch, depth, Adam-β, attention
window, small-init alone, and weight-decay all stay diminishing. **Repetition is falsified**
(strongly diminishing — overfitting). **Grokfast is antagonistic**: its only J > 0 comes with
collapsed early learning and degraded competence, and it *destroys* the low-flat-LR effect when
combined (`lr_0p25_gf5` −0.34, `gf5_lr0p5` −0.44).

### CAVEAT (read this before believing it)

A flat low LR **mechanically sustains late-training gains** — it stops the optimizer from
"finishing" early, so the loss curve is de-front-loaded on the log-experience axis (the same
axis children are measured on). By our `J` metric and artifact controls this **is** the
genuine increasing-returns signature: real early learning, real baseline competence, and the
gain-rate genuinely rises late. But **whether this is the *same mechanism* as children's
increasing returns is open.** This is hypothesis generation, not a claim of equivalence — the
LM result may be an optimization-pacing phenomenon that happens to share the J-signature. The
decisive disambiguating experiment (a late-introduced, frequency-matched held-out lexicon) is
not yet run (see Ideas).

---

## Per-family results

`J` = returns metric (>0 increasing); `r_e` = early-half gain-rate (artifact control, want ≥0.4);
`surp` = final mean surprisal (competence, baseline 5.88, want ≤6.1). ✅ = passes both controls.

### Learning rate — flat, global magnitude (THE lever, dose-response)
| variant | LR× | J | r_e | surp | |
|---|---|---|---|---|---|
| lr_0p03 | 0.03 | +1.37 | 0.49 | 6.35 | undertrain (artifact — surp>6.1) |
| lr_0p05 | 0.05 | +1.29 | 0.63 | 6.29 | undertrain (artifact) |
| lr_0p07 | 0.07 | +1.25 | 0.73 | 5.99 | ✅ genuine (frontier) |
| **lr_0p1** | 0.10 | **+1.11** | 0.77 | 5.92 | ✅ genuine (frontier) |
| lr_0p15 | 0.15 | +0.94 | 0.81 | 5.82 | ✅ genuine |
| lr_0p25 | 0.25 | +0.46 | 0.85 | 5.90 | ✅ genuine |
| lr_0p35 | 0.35 | +0.28 | 0.86 | 5.86 | ✅ genuine |
| lr_0p5 | 0.50 | −0.01 | 0.93 | 5.92 | ✅ flat (≈unit accumulator) |
| constlr | 1.0 | −0.47 | 0.95 | 6.01 | diminishing |
| lr_2 | 2.0 | −0.62 | 1.08 | 5.81 | diminishing |
| lr_4 | 4.0 | −0.58 | 1.50 | 5.91 | diminishing |

### LR × small-init (stacking)
| variant | J | r_e | surp | |
|---|---|---|---|---|
| **lr_0p1_init0p3** | **+1.28** | 0.77 | 5.76 | ✅ CHAMPION (best J at best competence) |
| lr_0p25_init0p3 | +0.46 | 0.89 | 5.69 | ✅ best competence overall |
| lr_0p25_gf5 | −0.34 | 0.84 | 6.81 | grokfast destroys it |

### Decay × low global magnitude (the 2×2 conjunction test)
| variant | J | r_e | surp | |
|---|---|---|---|---|
| decaylow_0p1 | −0.25 | 0.95 | 5.92 | diminishing at baseline competence → **flat is necessary** |
| decaylow_0p25 | −0.40 | 0.94 | 5.94 | diminishing |

### Matrix-LR alone (dissociation — decay left on)
| variant | J | r_e | surp | |
|---|---|---|---|---|
| matlr_0p01 | −0.51 | 0.94 | 5.93 | low matrix-LR + decay = diminishing → not the matrix-LR per se |
| matlr_0p08 | −0.34 | 0.95 | 5.86 | diminishing |

### Schedule shape (NULL)
| variant | J | r_e | surp |
|---|---|---|---|
| wsd | −0.29 | 0.95 | 5.84 |
| warmdown_0 | −0.31 | 0.95 | 5.82 |
| warmup_0p2 | −0.33 | 0.94 | 5.85 |
| warmdown_1 | −0.35 | 0.95 | 5.87 |
| finallr_0p3 | −0.37 | 0.95 | 5.95 |

### Grokfast (slow-gradient amplification) — antagonistic / artifact-bound
| variant | J | r_e | surp | |
|---|---|---|---|---|
| gf_flat | −0.11 | 0.95 | 5.98 | best grokfast competence, still diminishing |
| gf_lamb10 | +0.01 | 0.65 | 6.26 | borderline, competence already degrading |
| gf_lamb20 | +0.42 | 0.44 | 6.70 | artifact (early learning + competence degrading) |
| gf_lamb40 | +1.69 | 0.06 | 6.40 | artifact (r_early ≈ 0) |
| gf_lamb80 | +2.17 | −0.24 | 6.73 | artifact (model gets *worse* early) |
| gf_rep4_l10 | −0.19 | 0.82 | 6.14 | rescues rep4 (−1.09) toward baseline but still diminishing |
| (lamb1/2/5, alpha, nowd, combos) | −0.2..−0.5 | — | — | diminishing |

*Grokfast flattens diminishing returns at moderate λ but never reaches genuine increasing
returns at baseline competence; the J>0 region is the slow-start artifact.*

### Repetition (FALSIFIED — overfitting, strongly diminishing)
| variant | J | surp | note |
|---|---|---|---|
| rep2_decay | −0.80 | 6.02 | 2 passes |
| rep4_decay | −1.20 | 6.56 | gets worse with passes |
| rep8_constlr | −1.39 | 6.52 | strongly diminishing |
| s1_rep8_constlr | −1.39 | 5.99 | (disjoint-seed control, same pattern) |

### Depth / size, batch, optimizer, init, weight-decay (all NULL)
| variant | J | r_e | surp | note |
|---|---|---|---|---|
| arch_d4 | −0.35 | 0.90 | 6.10 | smaller → more diminishing |
| arch_d2 | −0.52 | 0.89 | 6.49 | smaller → worse |
| batch_2e18 | −0.38 | 0.94 | 5.88 | larger batch diminishing |
| batch_2e20 | −0.32 | 0.95 | 5.92 | larger batch diminishing |
| betas_high | −0.44 | 0.94 | 5.83 | Adam β=(.95,.99) |
| window_full | −0.38 | 0.94 | 5.89 | all-long attention window |
| init_0p1 | −0.18 | 0.99 | 5.78 | small-init **alone** does NOT give increasing returns |
| init_0p3 | −0.19 | 1.04 | 5.75 | (rich-regime hypothesis did not pan out alone) |
| init_3 | −0.03 | 0.91 | 5.99 | large-init/lazy ≈ flat, no acceleration |
| wd_mult0 | −0.30 | 0.94 | 5.75 | no weight decay |
| wd_mult3 | −0.60 | 1.06 | 6.01 | more WD → worse |
| wd_mult10 | −0.71 | 1.11 | 5.93 | more WD → worse |

**Crashes (not retried, rc=1):** `deep_narrow`, `composite` (depth-16) and `gf5_d12`
(depth-12 + grokfast EMA buffers) OOM'd on the A40 at batch 64 — depth-up at fixed width is
untestable on this hardware without a smaller batch. `batch_2e16` was a transient GPU-race
casualty (its OOM trace named `lr_0p07`'s pid; `lr_0p07` won and completed).

---

## Ideas for Mike

1. **The decisive learning-to-learn experiment** (not yet run): introduce a *frequency-matched
   held-out lexicon late* in training and ask whether the low-flat-LR model acquires it *faster*
   than an equal-experience baseline. This separates genuine learning-to-learn from
   optimization-pacing — it is the experiment that would let us *claim* the child mechanism
   rather than just match the J-signature. Needs dataloader work.
2. **Tighter flat-vs-decay × magnitude grid** to map the full mechanism surface (we have the
   four corners; fill the interior, and test flat × {0.05, 0.07} × {with, without small-init}).
3. **Curriculum-of-capabilities** (data ordered foundational→dependent) — the "Domino"
   compounding-skill hypothesis; the most theory-aligned untried genuine lever.
4. **μP** width-scaled init+LR (pairs with Muon) — maximal-feature-learning regime + LR transfer.
5. **Sophia** optimizer swap; **ramped critical-batch** (small→large).
6. **Depth-over-width at smaller batch** — depth-16/12 OOM'd here; rerun at batch 16–32 to test
   the depth-emergence hypothesis that this hardware blocked.

Patches available (composable via comma in a manifest, env via `CDI_*`): `patch_lr` (flat LR
×mult), `patch_initscale` (init scale), `patch_const` (any module constant), `patch_arch`
(depth), `patch_grokfast`, repetition patches. Scorer: `returns_score.R`. All run data on
ccn2 `/data2/mcfrank/autoresearch/runs_bank/`; scores in `bank_results.tsv`.

---

## overnight log
- 2026-06-18 03:34 UTC: grokfast bank scored (33 variants total). **Grokfast crosses J=0 at high lamb** — lamb sweep 1/2/5/10/20 -> J −0.52 / −0.21 / −0.27 / **+0.01** / **+0.42**; gf_lamb20 shows INCREASING returns. Caveat: final_surp rises with lamb (6.70 @ lamb20 vs 5.88 baseline) — competence cost, possibly slow-start-then-catch-up not pure acceleration. Grokfast also rescues repetition (gf_rep4 J=−0.05 vs rep4 −1.09; gf_rep2 −0.21 vs rep2 −0.79). Combos: gf_flat −0.11. -> Launched bank gf2: lamb{15,40,80} + gf_flat/nowd/rep at strong lamb to test if J>0 holds WITH good competence.
- 2026-06-18 04:23 UTC: gf2 scored (41 total). **J>0 at high lamb is a SLOW-START ARTIFACT**: gf_lamb40 r_early=0.06, gf_lamb80 r_early=−0.24 (model barely learns/gets worse early -> mechanical 'acceleration'). No setting hits J≥0 at baseline competence (5.88); all grokfast final_surp worse. Real signal: grokfast FLATTENS diminishing returns in lamb 2–10 with early learning intact (gf_lamb10 J=+0.01, r_early=0.65, surp 6.26), gf_flat best competence (5.98, J=−0.11). -> gf3 confirmation: gf_flat_l5/l8 (sweet-spot?), lamb30 (boundary), gf_flat_l40.
- 2026-06-18 05:40 UTC: **REAL HIT — low flat LR gives genuine increasing returns.** lr_0p25 (flat LR x0.25): J=+0.46, r_early=0.85, final_surp=5.90 (= baseline 5.88) — passes ALL artifact controls (healthy early learning AND baseline competence), unlike grokfast (which needed early-stall). LR sweep monotonic: J = +0.46/-0.005/-0.47/-0.62/-0.58 at LR x0.25/0.5/1/2/4 -> lower LR = increasing returns, higher = diminishing. Mechanism: flat low LR ~constant per-step gain = increasing returns on log-experience axis, still reaches baseline competence; decay/high LR slows late. gf3 confirmed grokfast sweet-spot is artifact-bound (gf_flat_l40 r_early=0.06). arch: smaller (d2/d4) more diminishing. -> Launched rich bank (small-init/WSD/deep-narrow/composite) + lrlow follow-up (LR 0.1/0.15/0.35 + lr0.25 x init/grokfast).
- 2026-06-18 07:25 UTC: **LR dose-response confirmed — clean & monotonic.** Low flat LR -> increasing returns at baseline competence, all artifact controls passed: J = +1.11/+0.94/+0.46/+0.28/-0.01 at LR x0.1/0.15/0.25/0.35/0.5 (r_early 0.77-0.93 healthy, final_surp 5.8-5.9 = baseline). **Small-init STACKS**: lr_0p25_init0p3 J=+0.46 with best competence of all (5.69 < baseline 5.88). **Grokfast antagonistic**: lr_0p25_gf5 J=-0.34, surp 6.81. Small-init ALONE / WSD / depth do NOT give increasing returns (init J~-0.18, wsd -0.29) -> it's LR, not rich-regime/schedule. -> Launched lrfloor (LR 0.03/0.05/0.07 + lr0.1_init0.3: find competence floor) [2nd/final follow-up] + started batch campaign. Brief GPU race (batch_2e16 vs lr_0p07) self-resolved when deep_narrow/composite finished; both dispatchers now clean.
- 2026-06-18 08:16 UTC: **LR dose-response COMPLETE — floor found.** Full curve J vs LR x{0.03..4}: +1.37/+1.29/+1.25/+1.11/+0.94/+0.46/+0.28/-0.01/-0.47/-0.62/-0.58. **Genuine frontier = LR x0.07-0.1** (J +1.25/+1.11 at baseline competence 5.99/5.92, r_early 0.73/0.77). Below x0.07 J keeps rising but COMPETENCE BREAKS (lr_0p05/0p03 final_surp 6.29/6.35>6.1, r_early ->0.49) = undertrain artifact, correctly excluded by the control. **CHAMPION = lr_0p1_init0p3: J=+1.28, r_early=0.77, final_surp=5.76 (BEST competence, < baseline 5.88)** -- low-LR+small-init stacks. Larger batch (2e18/2e20) diminishing (J -0.37/-0.32). deep_narrow/composite/batch_2e16 OOM (depth16 too big for A40; 2e16 was race casualty vs lr_0p07) -- rc=1, not retried. -> launched sched campaign (warmup/warmdown/WSD variants).
- 2026-06-18 09:05 UTC: **Schedule family = clean NULL.** warmup_0p2/warmdown_0/warmdown_1/finallr_0p3/wsd all J=-0.29..-0.37 at baseline competence (5.82-5.95), r_early ~0.95 -> schedule SHAPE (warmup/decay/WSD/final-LR) does NOT give increasing returns. Sharpens headline: it is LR MAGNITUDE, not schedule timing. Launched optlr (matrix-LR 0.01/0.08 + betas + window) to cross-check LR via the matrix-LR knob. Fixed controller bug: 'pgrep -f bank_dispatch' self-matched the monitoring shell (always >=1) -> had silently blocked optlr launch; now use python-proc check.
- 2026-06-18 09:55 UTC: **optlr family -> DISSOCIATION that refines the mechanism.** matlr_0p01 (matrix-LR x0.25 but decay STILL ON) J=-0.51 (diminishing, worse than baseline) -- low matrix-LR alone does NOT reproduce the effect. So increasing returns requires the CONJUNCTION the patch_lr runs had: FLAT (non-decaying) schedule AND low GLOBAL magnitude. Decay dominates (LR->0 late kills late gains -> J<0 at any magnitude); flat-at-full-magnitude (constlr/warmdown_0 ~ -0.47/-0.31) also insufficient. Both necessary. betas_high J=-0.44, window_full J=-0.38 (nulls). -> Launched gfsize (grokfast x depth: gf5_d4/d12/...), the LAST manifest. Then finalize.
- 2026-06-18 10:42 UTC: gfsize (LAST manifest) scored. gf5_d12 OOM (rc=1, depth12+grokfast EMA buffers too big for A40; not retried). All campaign families now complete; box idle. Awaiting 14:00 UTC finalize.
- 2026-06-18 10:43 UTC: campaign manifests all complete (77 variants). Launched ONE focused mechanism test on the idle box (fills the 2x2 missing cell): decaylow_0p25/0p1 = DECAY schedule x low GLOBAL magnitude (all 4 LRs scaled x0.25/x0.1, warmdown stays on). Predict diminishing if flat-schedule is necessary. Then finalize (campaign exhausted; no need to idle to 14:00).
- 2026-06-18 11:32 UTC: **FINALIZED.** Conjunction confirmed: decaylow_0p1/0p25 J=-0.25/-0.40 at baseline competence (5.92/5.94, r_early 0.95) -> decay+low-global is DIMINISHING, so the 2x2 is complete and FLAT schedule is necessary (flat+low is the only increasing-returns cell). 79 variants scored. Report rewritten; 3 plots (fig_J_vs_LR, fig_J_vs_finalsurp, fig_bestJ_by_family) generated. Verdict: low FLAT LR (frontier x0.07-0.1) is the one robust genuine increasing-returns lever; champion lr_0p1_init0p3 J=+1.28 @ surp 5.76; mechanism = conjunction of flat schedule + low global magnitude; all other families null/antagonistic. Caveat recorded (genuine by our J metric+controls, but child-mechanism equivalence is open). Task #22 complete.
