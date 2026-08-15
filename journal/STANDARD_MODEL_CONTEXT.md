# The "standard model" projects — orientation for a new session

**What this is.** A single page of context for someone (or some agent) starting cold on
either of the two repos that came out of this work. It explains what the research program
is, what was tried and what came of it, which repo now holds what, and where the detailed
records live. It is deliberately short; the primary sources are listed at the end.

**This file is maintained in BOTH repos and should stay identical in both.** If you change
it in one, copy it to the other in the same session.

---

## The research program

The umbrella question: **can children's word learning and language-model training be put
on a single quantitative axis, and if so, do they behave the same way?**

The shared object is an **accumulator model** of vocabulary growth. A child's latent
ability rises with the log of accumulated input, and a word is produced once ability
exceeds that word's difficulty:

$$\theta_i(t) = \xi_i + \kappa_i \log(t/a_0) + \log H, \qquad P(\text{produce } j) = \mathrm{logit}^{-1}(\theta_i(t) - \delta_j)$$

Three parameters carry the whole program. $\xi_i$ is a child's **efficiency** (level),
$\kappa_i$ their **acceleration** (how fast ability grows per log unit of input), and
$\delta_j$ a word's **difficulty**. The null of interest is $\kappa = 1$: the "pure
accumulator", in which learning is proportional to input and nothing about the learner
changes. Children come out around $\kappa \approx 11$–13. Language models come out around
$\kappa \approx 1.2$.

That gap is the headline of the acceleration paper. A second, related program asks what
*drives* $\xi_i$ and $\kappa_i$ — input quantity, processing speed, demographics — and
became the input paper.

## How the work actually went, in five movements

The detailed log is `journal/experiments.md` (41 numbered experiments) and
`journal/experiments_llm.md` (L1–L7, the GPU work). The arc:

**1. Cross-sectional beginnings, and why they failed (exp. 1–13).** Early fits were
cross-sectional. They recover *a* $\kappa$, but one snapshot per child cannot distinguish
"children differ in level" from "children differ in slope" — the two generate identical
clouds. This is why the program moved to longitudinal data, and the point survives as a
figure in the SI.

**2. The nested model family and its identifiability problems (exp. 14–27).** A ladder
M0–M5 was built and compared by LOO. Most of the work here was fighting parameters that
would not identify: a developmental-onset term `s_i` consumed four experiments (21–23)
before being excised, and $\sigma_r$ sensitivity took three (4, 18, 36). The recurring
lesson — a parameter that fits better but does not identify is worse than no parameter —
shaped everything after.

**3. Input and processing channels (exp. 8–9, 30–39).** Whether input quantity (BabyView,
SEEDLingS, LENA) and processing speed (looking-while-listening RT) explain individual
differences. Produced the joint `io_proc` models and the variance partition. **This became
the input paper and stays in `standard_model_2`.**

**4. Language models as a comparison system (L1–L7).** GPT-2 trained on CHILDES at
increasing data budgets, with per-word learning curves fitted the way Chang & Bergen
(2022) do. The key methodological move was the **developmental ladder** — separate models
each trained to convergence on a larger subsample — rather than reading checkpoints from
one training run, because the latter confounds data quantity with optimization time.
A composition control (BabyLM, ClimbMix) followed.

**5. The Bayesian by-dataset ladder, and the acceleration paper (exp. 40–41).** The
current `bayes_long` family: M0 ($\kappa=1$) → M1 (population $\kappa$) → M2 (per-child
$\xi_i$) → M3 (per-child $\kappa_i$), fitted separately to five longitudinal samples in
three languages. **This became the acceleration paper and now lives in `acceleration`.**

## Things that were settled, and are worth not relitigating

- **$\kappa$ is not an artifact of the CDI.** Refitting on deliberately narrowed item sets
  moves $\kappa$ by 9–22%, upward, where a pure instrument-rescaling account requires an
  ~80% fall. Size-matched random controls sit within a few percent of baseline.
- **$\kappa$ is not an artifact of the estimator.** Running the language models' own
  four-parameter-logistic estimator on child data, with no IRT model involved, gives 6.3–10.2
  against 1.16 for LMs.
- **Acceleration forecasts; per-child acceleration does not.** Prospective cross-validation
  on held-out administrations: pinning $\kappa=1$ costs +92 to +294 nats per child, while
  per-child $\kappa_i$ never beats the shared value. Both are true and the paper says so.
- **$\sigma_b$ is inflated by sparse data.** Two administrations give a per-child slope with
  zero residual df, so the main analysis uses 3+.
- **Constant input rate is defensible over 8–36 months.** Within-child slopes in four
  longitudinal input datasets are flat.

## Which repo holds what

| | `acceleration` | `standard_model_2` |
|---|---|---|
| paper | children vs. LMs; $\kappa \gg 1$ vs $\kappa \approx 1$ | what drives $\xi_i$ and $\kappa_i$ |
| models | `bayes_long` M0–M3, 2PL, forward CV | `io_proc`, `proc_dp`, glmer ladder |
| studies | `bayes_long`, `llm`, `gamlss` | `io_proc_glmer`, `cross_sectional_demographics`, `input_estimation` |
| fits | `fits/bayes_long`, `fits/llm` | `fits/io_anchored`, `fits/glmer_ladder`, `fits/summaries` |

Both repos carry a copy of the experiment logs, because the history is genuinely shared:
the acceleration paper's model family grew out of experiments that also produced the input
paper's. **The copies are marked as copies at the top and are not maintained in parallel —**
see the header of each log for which repo is authoritative for new entries.

## Conventions that outlive any one session

- **Prose is human-written.** Agents do not edit prose, captions, section headers, or
  citation markers in `.qmd`/`.tex`. Code chunks, YAML, `.bib`, and helper scripts are fair
  game. Manuscript changes go through a PR on a git worktree, never a branch checkout in
  the author's tree. See `~/.claude/skills/rmd-qmd-manuscript-collab.md`.
- **Bibliography entries are verified against Crossref/arXiv, never written from memory.**
- **Cluster work**: Sherlock (`-p mcfrank` owned node, `-p owners` preemptible with
  `--requeue`) for Stan; ccn2 A40s and Marlowe for GPU. See
  `~/.claude/skills/sherlock-stan-fitting.md`.
- **Fits are gitignored**; small text exports (`summaries/*.csv`) and `paper/cache/*.rds`
  are committed so the manuscript builds without cluster access.

## Primary sources

| file | what it holds |
|---|---|
| `journal/experiments.md` | 41 numbered experiments: hypothesis, what ran, what came of it |
| `journal/experiments_llm.md` | L1–L7, the GPU/LM work |
| `journal/PROVENANCE.md` | data sources and how they were pulled |
| `journal/paper_models_provenance.md` | fit → cache → manuscript chain (audited 2026-06-13; **describes the input paper**) |
| `journal/HANDOFF.md` | a mid-2026 state snapshot; historical |
| `studies/*/README.md` | per-study specifics |

**A caution on the older docs.** `paper_models_provenance.md` and `HANDOFF.md` were
written for the long/input paper before the split and were last audited in June 2026.
Treat them as historical for anything concerning the acceleration paper.
