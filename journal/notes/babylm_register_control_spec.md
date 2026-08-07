# Register / composition control — spec (the "affirmative leg" for L5)

**Status:** specced 2026-07-17, decisions locked with MCF. Prereqs open (data staging +
CDI-frequency go/no-go). Launch after the CHILDES archive re-training sweep frees the ccn2 GPUs.

## Why

L4/L5 showed LM "individuals" (seed = init + input stream) converge with input — between-instance
σ_κ ≈ 0.10 (CV 8%) vs children's ≈ 3.5 (CV 35%). But L5 split **CHILDES** halves, which are
distributionally near-identical → **low power to *find* divergence**. This is the affirmative leg:
does *compositionally different* disjoint data induce child-like between-instance spread? MCF is
committed to genuine interpretation — if it diverges, we report it (LM individual-variation is
composition-driven), not just confirm the null.

## Design (locked)

**Fixed — identical measurement to L6:** `GPT2_CHILDES` tokenizer; CHILDES-val CDI probe (611
words); GPT-2-small (125.8M); early-stop patience 4 / ceiling 40 epochs; nested-cumulative budgets;
best-val readout; budget in **whitespace words**.

**Varied — two maximally-different corpora × disjoint subsets × seed:**
- **6 datasets** = **3 BabyLM (CHILDES excluded)** + **3 ClimbMix**, **random-disjoint within each
  corpus**. The headline lever is the **between-corpus** contrast (child-oriented vs web); the 3
  within-corpus subsets give the sample-variance baseline. Variance ladder:
  **σ_seed (init) < σ_subset (sample) < σ_corpus (composition).**
- **8 seeds** per dataset (MCF: fine with the cost). **12 rungs to 24M**:
  `[0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24]M` (supports per-word 4-PL κ for panel A + the
  aggregate trajectory).
- **Total: 6 × 8 × 12 = 576 runs.**

**Why exclude CHILDES from BabyLM:** BabyLM-2024 100M is ~29% CHILDES (CHILDES ~29M, Gutenberg ~26M,
OpenSubtitles ~20M, Simple Wiki ~15M, BNC ~8M, Switchboard ~1M). A naive subsample would re-use our
main-experiment data (~7M CHILDES/subset). Dropping CHILDES leaves ~70M non-CHILDES → 3 disjoint
~23M subsets (reaches the 24M top). ClimbMix (400B tok, `karpathy/climbmix-400b-shuffle`) trivially
gives 3 disjoint 25M slices.

## Go/no-go — CDI-word frequency check (no GPU; MCF's idea)

LM per-CDI-word surprisal is frequency-driven (≈ log-freq, r≈−0.6). So **before spending GPU time,
prep the subsets and tabulate each of the 611 CDI words' frequency in each subset.** Confirms the
probe has a learnable signal in ClimbMix (web text far out-of-domain for a CHILDES-val CDI probe —
floor-bound words would make κ unfittable), and simultaneously **quantifies the composition contrast**
(how different are the CDI-word frequency profiles across corpora/subsets). Proceed to training only
if CDI-word coverage is adequate in both corpora.

## Analysis

- Per (corpus, subset, seed) ladder: aggregate dev-axis κ (0.434 × slope of −mean CDI surprisal vs
  log10-budget) + per-word 4-PL κ across 12 rungs.
- `kappa ~ 1 + corpus + (subset|corpus) + (seed|subset)` → σ_corpus, σ_subset, σ_seed.
- Overlay γ-trajectory spaghetti vs CHILDES + children; compare all σ's to CHILDES σ_seed≈0.10 and
  children σ_κ≈3.5.
- **Pre-committed interpretation:** converge → null bulletproofed ("even child-oriented vs web data,
  LM individuals develop at the same rate"); fan out → composition-driven LM variation, reported.

## Compute

576 runs, 12 rungs to 24M. ~486 GPU-hr × (48/32 ladders) ≈ **~730 GPU-hr**; +20–40% for higher
fertility under the CHILDES tokenizer (ClimbMix web vocab especially) → **~900–1000 GPU-hr ≈ ~5 days
at 8-wide** on ccn2. MCF: nothing else important on the box. Runs after the CHILDES archive sweep.
Pipeline: `worker.sh` reparameterizes — new chunk dir, HF path `composition-control/{corpus}/S{k}/seed{s}/rung{M}`.

## Caveats

- ClimbMix / non-CHILDES models scored on the CHILDES-val CDI probe sit at higher absolute surprisal
  (domain mismatch) — expected; we read **slope + spread**, not level. (Gated by the freq check above.)
- 6 datasets ⇒ σ_corpus is a 2-level contrast, σ_subset ~3 levels/corpus — wide CIs; the visual +
  variance components carry it with honest caveats.

## Open / prereqs

1. **Stage BabyLM-2024** + slice ClimbMix; confirm exact non-CHILDES token counts.
2. **CDI-word frequency go/no-go** across subsets (above).
3. **Oak** re-transfer of babyview (tar in progress; retry after the partial is cleared — git-annex
   read-only objects needed `chmod -R u+w` before rm).
4. Archive the 576 control models' weights (HF, +~140 GB) or just CSVs+recipe — decide later.
