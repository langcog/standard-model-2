# Marlowe project application — draft

**Drafted 2026-05-19**, revised after design discussion. This is a
draft, not a submitted application. Edit before submitting via
<https://datascience.stanford.edu/marlowe/project-application-guide>.

The full project ask is **~920 H100-hours** for one cycle: a 2-run
portability pilot followed by a 90-run main study (described
below). The pilot is the first phase, not a separate application.

---

## Section 1 — Project & PI information

- **Title.** Sources of between-run variability in language-model
  vocabulary learning: a comparison to variability in children
- **PI.** Michael C. Frank, Department of Psychology, Stanford
- **Group.** Language and Cognition Lab (`langcog`)
- **Project type.** Full project (one cycle), with a 2-run
  portability pilot as the first phase
- **Urgent timeline.** None binding. Project outcomes feed a
  follow-up manuscript to Frank, Tan, & Feng (in preparation) on
  individual variability in early word learning — target submission
  Q4 2026.

## Section 1b — Abstract

Children differ markedly from one another in how quickly they
acquire words. Do language models? For children we measure
between-child variability in their vocabulary, but for LMs there is
very little effort to measure this kind of variability. We attempt
to do this. A prior run on Sherlock already trained 3 GPT-2-small
seeds on child-directed speech (25M words). In this current
experiment, we will expand this to measure seed-related,
architecture-related, and data-related variability, training GPT-2
instances across a range of different conditions and measuring
between-run variation in vocabulary knowledge. These runs will be
compared to variability in human learners, providing insight into
similarities and differences between human and machine language
acquisition.

## Section 2 — Computational Suitability Statement (CSS)

### Background and significance

Children's early vocabulary grows very rapidly, with large
variation between children. Some children know a few words by their
first birthday; others know hundreds. This variability is one of
the most robust features of early language development, and
understanding where it comes from has theoretical, practical, and
clinical importance: it bears on the mechanisms of word learning,
it predicts later academic achievement, and it shapes how we
identify and support children with language delays.

A common view in both developmental psychology and machine learning
treats word learning as a process of accumulation — more exposure
to language yields more vocabulary in a predictable, additive way.
This view is the conceptual bridge to large language models, whose
learning behavior is famously well-described by smooth scaling laws
in the amount of training data they see. If children are pure
accumulators, then modern language models may be reasonable
computational accounts of early language acquisition; if they are
not, then the apparent similarity is superficial.

Our prior work suggests the second answer. Across large samples of
children's vocabulary growth, we find that vocabulary acquisition
*accelerates* over development in a way that pure accumulation
cannot explain, and that variation in how much language children
hear accounts for only about 10% of the variability between
children. Most between-child variability in vocabulary growth comes
from sources other than input quantity.

Language models, in contrast, *do* appear to behave like pure
accumulators: their learning is well-fit by smooth scaling
relationships. But there is a missing piece in this comparison. For
children, we measure variability *between individuals*. For
language models, almost nobody has measured the corresponding
quantity. The typical evaluation reports performance for a single
trained model (or averages over a few seeds), without asking how
much two language models trained under nominally identical
conditions actually differ in what they have learned. Without that
measurement, the apparent contrast between "variable children" and
"deterministic LMs" cannot be sharpened into a real claim.

This project measures between-LM variability in vocabulary
learning, decomposes it into contributions from different sources
of variation, and compares it directly to the between-child
variability we measure in CDI data (over 5,000 children).

### Building on existing infrastructure

A prior run on Sherlock has already trained 3 GPT-2-small models on
the same set of child-directed speech (a filtered 25M-token English
subset of the CHILDES corpus). All three models gave nearly
identical vocabulary-learning curves on the 611 words in the
MacArthur-Bates Communicative Development Inventory. This tells us
that seed-only variability at full data scale is small.

The next question is whether *other* sources of variation —
particularly variation in the data each model sees, and variation
in architecture — produce larger differences between models, and
how those differences compare to the between-child variability we
see in CDI data. Answering that requires training a population of
models, not just a few seeds.

The training and evaluation pipeline is already built, runs
end-to-end on Sherlock, and produced the 3 completed seeds. The
proposed Marlowe project ports it to Marlowe and scales it out.

### Phase 1: portability pilot (2 runs)

The first phase trains **2 GPT-2-small models** on two non-overlapping
10M-token chunks of CHILDES (the maximum number of fully disjoint
~10M chunks the corpus permits). Same architecture, same random
seed, same training recipe — the only intentional difference is the
data each model sees.

This 2-run pilot has two purposes:

1. **Validate the Sherlock-to-Marlowe port.** Confirm that the
   training pipeline reproduces our existing 3-seed CHILDES results
   on Marlowe before scaling to 90 runs.
2. **Decide the data substrate for Phase 2.** If the two 10M-CHILDES
   models give noticeably different vocabulary-learning curves, the
   90-run study's data-variation axis uses CHILDES at a slightly
   smaller scale (to support more disjoint subsamples). If they
   look identical, the data axis uses TinyDialogues (Feng et al.,
   2024) — a 200M-word synthetic child-directed-speech corpus that
   permits many non-overlapping training subsamples.

Compute: 2 runs × ~10 H100-hours = **~20 H100-hours**.

### Phase 2: main study (90 runs)

The main study trains 90 GPT-2 models in a factorial design that
isolates three sources of between-run variability:

- **Seed variability (30 runs).** Same training data, same
  architecture, different random seeds. Establishes the irreducible
  randomness in what a language model learns under fixed inputs.
- **Data variability (30 runs).** Same architecture, same random
  seed, different disjoint training subsamples (substrate chosen
  by the Phase 1 pilot). Measures how much of the between-run
  variability arises from differences in the training input — the
  language-model analog of how much between-child variability is
  driven by input differences.
- **Architecture variability (30 runs).** Same training data,
  varied architectural choices within the GPT-2 family (learning
  rate, model width, model depth). Measures how much between-model
  variability arises from architectural and hyperparameter choices
  that nominally do not change the learning problem.

For each of the 90 models, we save ~200 intermediate checkpoints
across training. At each checkpoint, we measure the probability the
model assigns to each of the 611 CDI words in context — yielding
per-model word-learning trajectories directly comparable to the
per-child vocabulary trajectories in our CDI data.

Compute: 90 runs × ~10 H100-hours per run = **~900 H100-hours**.

### Total project compute

| Phase | Runs | H100-hours |
|---|---|---|
| Pilot | 2 | ~20 |
| Main study | 90 | ~900 |
| **Total** | **92** | **~920** |

Per-run cost is anchored on our group's measured training profile
for GPT-2-small on a 25M-token corpus on Sherlock A100s. 10 H100-hours
per run is a conservative ceiling; at 10M tokens we expect
4–6 H100-hours per run in practice.

### Storage and I/O

- **Per-run checkpoint footprint.** ~200 checkpoints × ~500 MB
  (GPT-2-small + optimizer state) ≈ 100 GB per run.
- **Project total.** 92 runs × 100 GB ≈ **9 TB** project scratch.
- **Per-job I/O.** Marlowe's per-node NVMe (30 TB local) covers
  any individual training run with significant headroom.

### Software stack

PyTorch + HuggingFace transformers (matching the published Feng et
al. 2026 pipeline), with HuggingFace `datasets` for corpus
preprocessing. Post-training analysis (sigmoid fitting,
between-model variability estimates, comparison to CDI data) runs
in R on local hardware and is not part of the Marlowe ask. No exotic
dependencies; we will request a standard Marlowe Python/PyTorch
environment.

### Readiness

- **Training pipeline.** Already runs end-to-end on Sherlock; 3
  CHILDES seeds completed successfully under this pipeline (see
  `model/scripts/feng_eval/` in `standard-model-2`).
- **Data.** CHILDES corpus already preprocessed in our group; the
  10M-token chunks for the pilot will be drawn from the same 25M
  filtered English subset used in Feng et al. (2026). TinyDialogues
  (the candidate Phase 2 substrate) is publicly available from
  Feng et al. (2024).
- **Evaluation.** Per-word vocabulary-learning curves and the
  comparison to between-child variability are already implemented
  and validated in the parent project (`model/scripts/`).
- **Reproducibility.** Random seeds and training configurations
  will be logged for every run; per-run checkpoints, configs, and
  evaluation outputs will be retained for the full project lifetime.

### Why Marlowe

- Our group's primary GPU resource (Sherlock fairshare) is currently
  saturated by an unrelated Bayesian-modeling campaign for the
  parent project. Running the 90 training runs there would
  significantly delay both lines of work.
- The 90-run main study is embarrassingly parallel across runs.
  Marlowe's H100 DGX nodes give the throughput needed to complete
  the campaign in a manageable wall-clock window (target: ~4 weeks
  end-to-end).
- The 9 TB project storage and high-bandwidth NVMe simplify
  capture of the ~18,400 intermediate checkpoints the analysis
  requires.

## Section 6 — Expected outcomes & acknowledgments

- **Outputs.** (1) A decomposition of between-language-model
  variability in CDI-word vocabulary learning into contributions
  from random seed, training data, and architecture. (2) A direct
  comparison between this between-LM variability and the
  between-child variability we measure in the same 611 CDI words
  across 5,000+ children. (3) Conclusions about whether language
  models and children show qualitatively similar or different
  patterns of variability in early vocabulary acquisition.
- **Publication.** Results feed a follow-up manuscript to the
  in-preparation *Standard Model 2* paper (Frank, Tan, & Feng);
  expected target venues include *Open Mind* and *Cognitive
  Science*, with all Marlowe-supported runs acknowledged.
- **Code release.** Training and analysis code released under MIT
  license on submission.
- **Acknowledgment text.** "This work used the Marlowe GPU
  computational instrument at Stanford University, supported by the
  Office of the Vice Provost and Dean of Research."

---

## Notes for Mike (not part of the submission)

- The pipeline (`model/scripts/feng_eval/`) currently lives on the
  `claude/ecstatic-gould-922579` branch of `standard-model-2`. Once
  it's merged to master the CSS reference path becomes the master
  one — adjust before submitting.
- The 10 H100-hr/run estimate is from the Sherlock A100 runs at
  the full 25M-token corpus; at 10M tokens (Phase 1) it should be
  ~4–6 H100-hrs per run, so the 920 H100-hr ceiling has slack.
- "Section 6" numbering is from a partial reading of the project
  application guide; the live form may have intermediate sections
  (data management plan, security/compliance, fairshare
  integration) that aren't drafted here. Pull the live form before
  submitting.
- The full-scope ask is sized as one project allocation cycle. If
  reviewers prefer a phased ask, the natural split is Phase 1
  (~20 H100-hrs) approved immediately and Phase 2 (~900 H100-hrs)
  released on pilot completion — both within the same project, no
  re-application.
