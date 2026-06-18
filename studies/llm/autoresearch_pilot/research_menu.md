# Interventions to increase per-word acquisition acceleration (κ) — research menu

Generated 2026-06-17 by a background research agent (web search + primary-source
verification). Seeds the experiment bank. Two flavors kept strictly separate:

- **(i) SHARP** — a single word's curve is abrupt (high per-word κ); grokking / phase-transition-like.
- **(ii) LEARN-TO-LEARN** — the *rate* of acquisition improves over training, so words
  first encountered later are acquired faster/steeper than early ones (rate-of-rate).

**Master discriminator (free, run on every condition):** the **onset→κ correlation**
at matched frequency. Flavor (i) = high κ scattered across training with no onset
dependence. Flavor (ii) = **positive onset→κ slope** (later words steeper). This is
the model analogue of children's "later/harder words learned faster."

**Adversarial framing:** at the *aggregate*, loss **decelerates** over training
(Src 9, zero-sum learning) — the population rate-of-rate is *negative*. So flavor (ii)
is a local/minority effect fighting a headwind; a credible claim must beat a
same-frequency early-word control, not just exist. This matches the κ<1 / =1 / >1
returns trichotomy: LMs default to **diminishing returns (κ≤1)**; learning-to-learn
would be a subset flipping to **increasing returns (κ>1)**.

**Validity (good news):** mean *surprisal* is the continuous metric Schaeffer et al.
(Src 2) say should be smooth, so a high κ measured on surprisal **cannot** be a
thresholding artifact; and per-direction loss is sigmoidal with staggered inflections
hidden in smooth aggregate loss (Src 8) — so high *per-word* κ is expected and
recoverable even when total loss is smooth. Our method measures the right thing.

## Top 6–10 to run first (ranked; flavor tagged; cost)

1. **Per-word sigmoid-κ on surprisal [analysis, enabling]** — lock the readout; defines success.
2. **Weight-decay sweep {0,0.01,0.1,0.3,1.0} [trivial, i]** — most-supported sharpness lever (Src 1,6).
3. **Grokfast slow-gradient amplification [trivial, i]** — one line, ×40–50 grokking speedup (Src 5,7).
4. **Corpus-shrink + repetition vs single-pass-fresh [moderate, i + data-regime]** — likely *why our κ sits ≈1*: single-pass over huge fresh web data keeps every word in the smooth "data-surplus" regime that suppresses sharpness (Src 6,14). Strongest candidate to move κ off 1.
5. **Late-introduced held-out target lexicon [moderate, ii — decisive causal test]** — withhold a frequency-matched word set until post-induction-bump; compare its κ / exposures-to-threshold vs early-matched controls (Src 10,11).
6. **Induced burstiness + preserved Zipfian tail in packing [moderate, ii]** — pack target words to recur within-context; don't dedup-flatten the tail (Src 11).
7. **WSD vs cosine schedule [trivial, i — but a CONFOUND]** — decay-phase sharpening can *manufacture* both (i) and (ii); run only with a constant-LR control (Src 12).
8. **AoA→κ regression at matched frequency [analysis, ii — free]** — does later onset predict steeper κ controlling for freq/length/MLU? C&B never asked; our data can (Src 15, E2).
9. *(stretch)* **budget × L2 for ICL persistence [trivial, ii failure-mode]** — ICL is transient; over-training can *reverse* the later-vs-early effect (Src 13).

## Distinguishing (i) from (ii)
- **onset→κ correlation (#8)** — master discriminator, free on every condition.
- **late held-out lexicon (#5)** — the causal version (manipulate onset at fixed frequency).
- **WSD (#7)** — the confound to neutralize with a constant-LR control.
- **aggregate deceleration (Src 9)** is the null for (ii): report κ relative to a same-frequency early-word control.

## Key mechanism notes
- **Grokking knobs (flavor i):** weight decay (Src 1,6), AdamW decoupling (Src 1,3), Grokfast EMA of slow gradients (Src 5,7), StableMax/⊥Grad to avoid Softmax Collapse (Src 4), Critical-Data-Size / repetition regime (Src 6,14).
- **Learning-to-learn (flavor ii):** induction-head phase change = "majority of ICL acquired" at a bump, needs ≥2 layers, causally tied to architecture (Src 10); emerges only with **bursty + many-rare-class, Zipf≈1** data (Src 11); **transient** — decays to in-weights learning, L2 prolongs it (Src 13); moving the bump earlier maximizes the late-word population that benefits (Src 10,12).
- **What predicts AoA in LMs:** log-frequency dominates (R²≈0.91–0.94), MLU+, n-chars− (opposite of kids), concreteness null (Src 15) — frequency is the confound to control everywhere.

## Sources
1. Grokking review (Power et al. 2022 + follow-ups) — https://en.wikipedia.org/wiki/Grokking_(machine_learning)
2. Schaeffer et al. 2023, Emergent Abilities a Mirage? — https://arxiv.org/abs/2304.15004
3. Curriculum Learning pretraining (EACL 2026) — https://arxiv.org/abs/2506.11300
4. Grokking at the Edge of Numerical Stability (ICLR 2025) — https://arxiv.org/html/2501.04697v1
5. Grokfast (Lee et al.) — https://arxiv.org/abs/2405.20233
6. Critical Data Size of LMs from a Grokking Perspective — https://arxiv.org/html/2401.10463v1
7. Grokfast repo — https://github.com/ironjr/grokfast
8. Hidden Breakthroughs in LM Training (POLCA) — https://arxiv.org/html/2506.15872v1
9. Loss Deceleration & Zero-Sum Learning — https://arxiv.org/html/2506.05447v2
10. Olsson et al. 2022, In-context Learning & Induction Heads — https://arxiv.org/abs/2209.11895
11. Chan et al. 2022, Data Distributional Properties Drive ICL — https://arxiv.org/abs/2205.05055
12. Why Warmup the Learning Rate? (NeurIPS 2024) + WSD — https://arxiv.org/html/2410.05192v1
13. Singh et al. 2024, Transient Nature of Emergent ICL — https://arxiv.org/abs/2311.08360
14. Datasets, Documents, Repetitions — https://arxiv.org/html/2503.07879v2 ; Scaling Data-Constrained LMs — https://arxiv.org/abs/2305.16264
15. Chang & Bergen 2022, Word Acquisition in Neural LMs — https://arxiv.org/abs/2110.02406
