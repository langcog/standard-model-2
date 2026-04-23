# Formal model summaries and a proposed unified Bayesian model

## Goal

A simple, analytically tractable Bayesian model of vocabulary growth that
directly connects three numerical quantities in real units:

1. **Input quantity** — tokens of child-directed (and possibly overheard)
   speech per unit time ($r_i$, tokens/hour for child $i$).
2. **Per-word exposure** — expected tokens of word $j$ heard by child $i$ by
   age $a$: $E_{ij}(a) = r_i \cdot p_j \cdot H \cdot a$, where $p_j$ is the
   corpus probability of word $j$ and $H$ is waking hours per month.
3. **Acquisition** — whether child $i$ has produced word $j$ by age $a_i$,
   as measured by CDI.

The model must let us push real input estimates in and get real vocabulary
predictions out (not standardized IRT units).

---

## 1. McMurray (2007), *Science* — the "standard" accumulator

**Setup.** Each of $N$ words has a *time-to-acquisition* threshold
$\tau_j$ drawn from a population distribution $g(\tau)$. At each discrete
time step every word accumulates one point. Word $j$ is *learned* at the
first $t$ with $t \ge \tau_j$.

**Growth curve.**
$$L(t) = \int_0^t g(\tau)\,d\tau = G(t).$$

Acceleration ( $L''(t) > 0$ ) iff $g'(t) > 0$: as long as there are fewer
easy than moderately-hard words, growth accelerates. Gaussian $g$ with
mean $\mu$ produces the classic "vocabulary explosion" centered on $\mu$.

**Assumptions.** Deterministic, parallel, shared accumulation rate across
words, difficulty $\tau_j$ abstract (Gaussian by CLT over many independent
word properties). Frequency enters only implicitly through $\tau_j$.

---

## 2. Mitchell & McMurray (2009), *Cognitive Science* — leveraged learning

Continuous-time generalization of (1) with a feedback term $c$:
$$a(t) = t + c\,L(t), \quad \frac{dL}{dt} = \frac{g(t+cL)}{1 - c\,g(t+cL)}, \quad L(t) = G(t + cL).$$

**Key result.** Leverage shifts curves left/right but *cannot create*
acceleration that isn't already in $g$. Acceleration requires
$g'(\tau) > 0$ over the relevant range.

**Treatment of frequency.** Tried two mappings:
- Additive: difficulty $= f_1 - f_j$. Yields almost-step-function growth
  (too extreme).
- Multiplicative: difficulty $\propto 1/f_j$. With Zipf $s=1$, difficulties
  are uniformly distributed → *linear* growth, no acceleration.

This is the first serious attempt to plug real frequencies in and reveals
that frequency alone + deterministic accumulation doesn't reproduce the
spurt.

---

## 3. Mayor & Plunkett (2010), *CogSci Proc.* — "Are infants full of Zipf?"

Sharpens the critique. Stochastic version: a word $i$ with probability
$f(i)$ is heard every $\sim 1/f(i)$ time units; if a fixed number of
presentations $\kappa$ is needed, $T(i) \propto \kappa / f(i)$. Under Zipf
$f(i) \propto 1/i$, times-to-acquire $T(i) \propto i$ are uniform on rank,
so $L(t)$ is *linear*. Simulations with the CHILDES Parental Corpus (~24k
types, 2.6M tokens) confirm no spurt.

**Conclusion.** You *need* either (a) a change in learning capacity over
development or (b) non-Zipfian variability in word difficulty (phonology,
syntax, context) to produce acceleration. Their proposed fix: let the
per-presentation learning probability grow with age, e.g. $p(t) = (t/t_0)^3$.

---

## 4. Hidaka (2013), *PLOS ONE* — AoA distributions → learning process

Population-level approach: shape of the across-child AoA distribution
diagnoses the underlying learning process.

Three model families, each yielding a closed-form AoA distribution:

| Process | Dynamics | AoA distribution |
|---|---|---|
| Cumulative | $N$-accumulator, constant rate $\delta$ | Gamma($N, \delta$) |
| Rate-change | 1-accumulator, rate $\delta \cdot t^D$ | Weibull |
| Cumulative-and-rate-change | $N$-accumulator, rate grows | Weibull-Gamma |

Fit to MCDI monthly acquisition curves (654 words, ~1000 children each).
**Findings.** Closed-class words best fit by rate-change (Weibull);
nouns/verbs/adjectives best by cumulative (Gamma). Accumulation parameter
$N$ correlates with word frequency; change-of-rate parameter correlates
with imageability. BIC favored cumulative-and-rate-change overall.

---

## 5. Mollica & Piantadosi (2017), *Open Mind* — the waiting-time Bayesian model

The most directly useful precedent. Formalized Hidaka's cumulative model
as a Poisson process of *effective learning instances* (ELIs).

**Generative model for word $w$.**
- Accumulation begins at start age $s_w$.
- ELIs arrive as a Poisson process with rate $\lambda_w$ (ELIs/month).
- The word is acquired after $k_w$ ELIs.
- Acquisition age $\sim s_w + \text{Gamma}(k_w, \lambda_w)$.
- Proportion of children who know $w$ at age $a$:
$$F(a - s_w; k_w, \lambda_w) = \int_0^{a-s_w} \frac{t^{k_w - 1} e^{-t\lambda_w}\lambda_w^{k_w}}{\Gamma(k_w)}\,dt.$$
- Data model: $x_{aw} \sim \text{Binomial}(F(a-s_w;k_w,\lambda_w), N_a)$.

Priors: $k, \lambda \sim \text{Uniform}(0, 10^4)$, $s \sim \text{Uniform}(0, 10^3)$.
Fit per-word with JAGS across 13 languages.

**Key findings.** $k \approx 10$ ELIs, $\lambda \approx 0.5$ ELIs/month,
$s \approx 2$ months for comprehension. ELIs correlate only weakly with
log-frequency ($r \approx -0.14$ for $k$), meaning "an ELI is not just a
token."

**Limitations for our purposes.** No child-level parameters; $\lambda$ and
$k$ are jointly non-identifiable (only their product matters for the mean,
only their ratio $k/\lambda$ for variance/shape); no explicit link to
input rate $r_i$ or to corpus frequency $p_j$.

---

## 6. Kachergis, Marchman, & Frank (2021), *CDPS* — accumulator ↔ IRT

Reframes the family as a 1-PL Rasch model:
$$P(y_{ij}=1 \mid \theta_i, d_j) = \frac{1}{1 + e^{\theta_i + d_j}}.$$

Fit to Wordbank W&S (5,492 English children, 16–30 mo, 680 items) via
`lme4::glmer`. Item-level covariate: expected tokens/month $= p_j \cdot 1200 \cdot 12 \cdot 30.44$
(using the Sperry/Hart-Risley/Weisleder-Fernald pooled CDS rate of ~1200 tok/hr).

Progressive models:
- `m_1pl`: random intercepts for item and person.
- `m_tokmo`: + `tok_per_mo` item-level covariate.
- `m_lc`: + lexical class.
- `m2`: lexical class × lifetime_tokens.
- `m3`: lexical class × tok_per_mo + age (best balance).
- `m4`: three-way interactions.

**What worked.** IRT framing, fits to big data, shows age and lexical
class swamp raw frequency.

**What didn't.** Parameters live on logit scale, not in units of
tokens. You could not invert the mapping back to "how many tokens does
this child need for word $j$." Stan attempts foundered on the
numerical scale of token counts (hundreds of thousands per month) and
on posterior geometry around joint $(k,\lambda)$-like parameters.

---

## What the literature leaves unsolved

Across these six papers, no model simultaneously has:

1. A child-level exposure rate $r_i$ in tokens/hour.
2. A word-level frequency $p_j$ tied to corpus counts.
3. A threshold/ELI parameter that is interpretable in *tokens*.
4. Individual differences in both children and words.
5. Bayesian inference that actually runs.

Mollica & Piantadosi has (3) in ELI units but not (1) or (4). Kachergis
et al. has (1)–(4) but only in logit units. Mitchell & McMurray has (1)
and (2) but is deterministic and has no (4). This is the gap.

---

## Sketch of a unified model

Let $T_{ij}(a) = r_i \cdot p_j \cdot H \cdot a$ be the expected number of
tokens of word $j$ heard by child $i$ by age $a$, where

- $r_i$ = child $i$'s input rate (tokens of CDS per hour);
- $p_j$ = corpus probability of word $j$ (from CHILDES);
- $H$ = waking hours per month (≈ 365);
- $a$ = age in months.

**Option A — log-linear (simplest, likely tractable).**

Threshold per word in log-token units: $\psi_j = \log T^*_j$ where
$T^*_j$ is the expected number of tokens needed to learn word $j$. Then
$$P(\text{produces}_{ij} \mid r_i, a_i) = \Phi\!\left( \frac{\log r_i + \log p_j + \log H + \log a_i - \psi_j}{\sigma} \right).$$

This is a probit-IRT with ability $= \log(r_i \cdot a_i) + \log H$ and
difficulty $= \psi_j - \log p_j$. But because everything is on the
log-token scale, we can back-transform: $\psi_j$ is literally "log
expected tokens to acquisition" for word $j$. Priors on $\psi_j$ go on
the log scale (e.g. $\psi_j \sim N(\log 1000, 2)$), which fixes the
"big numbers" Stan problem.

**Option B — Poisson/Gamma waiting time (Mollica & Piantadosi extended).**

Each child has an ELI *filtering* rate $\alpha_i$: a fraction of tokens
become ELIs. Then for child $i$, word $j$:
- ELI rate $\lambda_{ij} = \alpha_i \cdot r_i \cdot p_j \cdot H$ per month.
- Threshold $k_j$ ELIs per word.
- P(learned by $a$) $= F_\text{Gamma}(a - s_j; k_j, \lambda_{ij})$.

This keeps the Mollica & Piantadosi interpretability but adds a child
factor ($\alpha_i, r_i$). Trade-off: harder to fit, more parameters,
$(k_j, \lambda_{ij})$ weakly identifiable.

**Option A is probably the right place to start.** It is equivalent to
the Mollica & Piantadosi model in the limit where the Gamma is
well-approximated by a lognormal (large $k$), but with an explicit
per-child input rate.

---

## Observables-to-parameters mapping (what we can ground)

| Quantity | Source | Estimate |
|---|---|---|
| $r_i$ distribution | Sperry/Hart-Risley/Weisleder-Fernald | $N(\sim 1200, \sim 840)$ tok/hr |
| $p_j$ | CHILDES American English | 26k cleaned types |
| $H$ | Assume 12 waking hrs × 30.44 days | 365 hrs/mo |
| $\psi_j$ prior | M&P: ~10 ELIs × unknown tokens/ELI | very wide |
| Age $a_i$ | Wordbank | 16–30 months |
| $y_{ij}$ | Wordbank CDI:WS | 5,492 children × 680 words |

---

## Open questions for the user

Recorded in the companion file `open_questions.md` — I need answers
before writing the Stan model and simulation code.
