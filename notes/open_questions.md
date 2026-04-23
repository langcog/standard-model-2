# Questions to resolve before writing the model

Grouped by how consequential the answer is. I'll hold off on coding until
you've weighed in on the starred ones.

## A. Scope and modeling choices

**★ A1. Model family.** I see three candidates (details in
`model_summaries.md`):

- **(A) Log-linear probit/logit IRT** with ability $= \log(r_i \cdot a_i)$
  and difficulty $= \psi_j - \log p_j$. Analytically closed, fits in Stan
  with routine priors on log-token scale, back-transforms to "expected
  tokens to acquisition" per word. Closest to your 2021 paper but with
  units.
- **(B) Mollica-&-Piantadosi-style cumulative Gamma** with a per-child
  ELI-filter rate $\alpha_i$. More mechanistically interpretable (ELIs,
  thresholds) but harder to fit, with known $(k,\lambda)$ ridge.
- **(C) Stochastic accumulator** simulated forward (non-closed-form),
  used for prediction/exploration only.

My default recommendation is (A) as the primary analytic/fitting model
and (C) as a forward simulator to cross-check. Do you want (B) instead,
or (A)+(B) both?

**★ A2. What failed before in Stan.** You mentioned "working with big
numbers" and "couldn't back-transform." Can you say which specifically:

- Divergences / non-identifiability between threshold and rate parameters?
- Numerical overflow when exponentiating token counts?
- Poor priors that forced you into standardized units?
- Slow sampling on Wordbank scale (5k children × 680 items)?

The fix depends on which of these bit you.

**A3. Individual differences.** What do children vary in?

- Input rate $r_i$ only (exogenous)?
- Input rate *and* a learning-efficiency multiplier $\alpha_i$?
- A single latent ability that absorbs both?

M&P collapse the two; Kachergis et al. have a latent $\theta_i$ but no
explicit $r_i$. Teasing them apart requires at least some children to
have measured input. Do we have per-child input measurements we can use
for anchor, or must $r_i$ be a latent variable with a population prior?

**A4. Word-level structure.** Should $\psi_j$ (or $k_j$) depend on:

- Corpus frequency $p_j$ only?
- Lexical class (nouns/verbs/function words/adjectives)?
- Phonological/semantic covariates?

Your 2021 paper found lexical class explains more than frequency. Do we
want to keep it in the baseline model?

## B. Data and scope

**B1. Wordbank subset.** English CDI:WS production, 16–30 mo,
5,492 children, 680 items — i.e. the exact dataset in your 2021 paper?
Or add comprehension / other languages?

**B2. CHILDES frequencies.** Use the cleaned
`childes_english_word_freq_cleaned_noHapaxes.csv` you already built? Any
reason to redo or to restrict to a specific age range of the input?

**B3. Input rate prior.** Use the Sperry/Hart-Risley/Weisleder-Fernald
pooled $N(1198, 840)$ for $r_i$, or something else?

**B4. Overheard speech.** Model CDS only (as in the 2021 paper), or try
the v2.0 CDS + ODS extension sketched in your appendix? The latter
requires an extra parameter (relative value of an ODS token) and I don't
see a dataset that identifies it; probably defer.

**B5. Start time $s$.** Include a Mollica-style start time
(child/word-specific age at which accumulation begins), or set it to 0?
Including it buys better fits at early ages but adds identifiability
headaches.

## C. Deliverables

**★ C1. Primary use case.** Which drives the design more:

- **Fit to Wordbank** and report posterior distributions on $\psi_j$,
  $r_i$, etc., in real units?
- **Simulate** — plug in an input distribution $r_i \sim p(r)$ and a
  frequency distribution $p_j$ and predict vocabulary-size distributions
  at each age?
- **Counterfactuals** — "if this child's input doubled, what happens?"

I think all three are end points, but which is the *acceptance test*
for the paper?

**C2. Observable app.** What should the reader be able to manipulate?
My guess: sliders for (mean input rate, variance of input rate, Zipf
exponent of $p_j$, threshold scale $\psi$) with live plots of
(vocabulary growth curves, age distribution for a chosen word, item-level
learning curves). Confirm or redirect.

**C3. Scale of fits.** Full Wordbank fit in Stan on 5k × 680 = 3.7M
Bernoulli observations is doable but heavy. Acceptable to start with a
smaller subset (say 1000 children × 200 items) to iterate quickly, then
scale?

## D. Specification details I'd otherwise have to guess

**D1. Threshold parametrization.** In Option A, is $\psi_j$ better
parametrized as (a) log-tokens-to-learn, or (b) log-ELI-rate per token?
They differ by how priors interact with $p_j$.

**D2. Likelihood.** Bernoulli with probit link (Option A), or the true
cumulative Gamma (Option B)? The probit is a good approximation of the
Gamma for $k \gtrsim 5$ (CLT) and is much faster.

**D3. Does the rate change with age?** M&P found constant $\lambda$ fits
fine. Your 2021 paper found older children accumulate faster than
younger. Bake in a $\lambda(a) = \lambda_0 \cdot a^\delta$ growth, or
keep it constant for the baseline?

---

Once you've answered A1, A2, C1 at minimum, I can write the Stan code
and a simulation notebook. The rest can be pinned down as we go.
