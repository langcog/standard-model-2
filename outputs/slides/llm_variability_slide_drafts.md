# Slide drafts: LLM variability (slides 45 & 46)

Drafted 2026-05-19 for insertion into `standard_model_pptx.pdf`.

Slide 45 replaces the currently-empty "Minimal variance for LLMs"
placeholder. Slide 46 is new (future-work). Both use the deck's
existing vocabulary ($\kappa$, $\sigma_\zeta$, "between-instance").

Numbers below are sourced from:
- `data/feng_2026/gpt2_childes_seed{0,42,123}_sigmoids.txt` (per-seed C&B fits, ecstatic-gould branch)
- `outputs/figs/longitudinal/chang_bergen_slope_summary.csv` (C&B BookCorpus reference)
- Long-format kid fits: $\sigma_\zeta \approx 3.5$ from `long_no_freq_slopes_si_signed`

---

## Slide 45 — Minimal variance for LLMs

**Title:** Minimal variance for LLMs

**Subtitle / first-line setup (small text):**
The kid-side $\sigma_\zeta \approx 3.5$ is between-*child*. The right
LM-side analog is between-*LM*. Two preliminary checks:

**Bullets:**

- **Between architectures (Chang & Bergen 2022, BookCorpus).** Four
  LMs — BERT, GPT-2, biLSTM, LSTM — give median per-word slopes of
  0.76, 0.81, 0.87, 0.96. Between-architecture SD ≈ 0.09.
- **Between seeds (this work, CHILDES-trained GPT-2 small).** Three
  seeds ({0, 42, 123}) on the same 24.5M-token CHILDES split give
  median per-word slopes of **0.72, 0.74, 0.74**. Between-seed SD
  ≈ 0.01.
- **The gap stays an order of magnitude even at the most permissive
  reading.** Pooling both sources: between-LM SD ≤ ~0.1, vs.
  $\sigma_\zeta \approx 3.5$ for children. That's a **~35×
  difference** in between-instance variability, on top of the
  ~10× difference in *mean* slope.

**Footer / caveat line (small text):**
Caveats: $n = 3$ seeds is small for a true bound; "between
architectures" mixes Transformer and recurrent families. Future
work (next slide) widens the sample.

**Suggested figure (optional, right-hand side):**
A small inset table or strip plot:

| Source | n LMs | Median $\kappa$ | between-LM SD |
|---|---|---|---|
| BookCorpus (C&B) | 4 | 0.76–0.96 | 0.09 |
| CHILDES seeds (this work) | 3 | 0.72–0.74 | 0.01 |
| **Children** | **5000** | **10.3** | **3.5** |

---

## Slide 46 — Future work: a proper test of LM variability

**Title:** Next: a proper test of LM variability

**Subtitle / first-line setup (small text):**
The preliminary checks are suggestive but underpowered. We have
applied for Marlowe (Stanford GPU cluster) compute to run a
formal between-LM variance decomposition.

**Bullets:**

- **Pilot (now, 2 runs).** Train GPT-2 small on the two disjoint
  10M-token chunks extractable from CHILDES, fix everything else.
  Asks: *does data identity move per-LM $\kappa$ at the
  CDI-completion scale?* Decides the data substrate for the main
  study.
- **Main study (90 runs, 30 per axis).** Three sources of between-LM
  variability, measured independently:
    - **Seed variability** — same data, same architecture, n=30
      seeds. Pins the irreducible randomness floor.
    - **Data variability** — varied subsamples (CHILDES or
      TinyDialogues per pilot outcome), fixed seed and architecture.
    - **Architecture variability** — within-GPT-2 sweep over
      learning rate, width, depth. Within-family, *not* between
      Transformer/RNN as in Chang & Bergen.
- **What the experiment can show.** Combined upper bound
  $\sigma_\kappa^{\text{LM}} = \sqrt{\sigma_\text{seed}^2 +
  \sigma_\text{data}^2 + \sigma_\text{arch}^2}$. Two clean outcomes:
    - $\sigma_\kappa^{\text{LM}} \ll \sigma_\zeta^{\text{kids}}$ →
      strengthens the structural-difference claim with a decomposed
      source-of-variance story.
    - $\sigma_\kappa^{\text{LM}}$ reaches kid scale on some axis →
      "LMs match human variability *if* you vary the right
      thing" — a richer story, identifies which axis matters.

**Footer / caveat line (small text):**
Pipeline already runs end-to-end on Sherlock (3 CHILDES seeds
above). The Marlowe ask is for portability + the throughput
needed for 90 parallel runs.

---

## Notes for Mike (not on the slides)

- The "between-architecture SD ≈ 0.09" on slide 45 mixes Transformer
  and recurrent families — that's a *generous* upper bound; the
  within-Transformer arch SD measured on slide 46's main study is
  almost certainly smaller. If you want a tighter pre-main-study
  bound, drop the biLSTM/LSTM from the C&B comparison and the SD
  becomes 0.04 (just BERT vs. GPT-2). I left both in for slide 45 to
  steelman the LM side.
- The pilot's "data substrate decision" framing on slide 46 mirrors
  what's in [`outputs/marlowe_application_draft.md`](../marlowe_application_draft.md)
  — keep these aligned if either evolves.
- $\sigma_\zeta$ vs $\sigma_\kappa$ notation: the deck uses
  $\sigma_\zeta$ for the kid-side per-child slope SD elsewhere. I
  used $\sigma_\kappa^{\text{LM}}$ on slide 46 to make the analogy
  explicit without overloading $\sigma_\zeta$. Adjust to whichever
  is cleaner for your audience.
