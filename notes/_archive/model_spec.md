# Model specification: Log-linear IRT with absolute units

## Resolved design choices

| # | Choice |
|---|---|
| A1 | Log-linear probit IRT (everything on log-token scale). |
| A2 | Fix for Stan convergence: all latent quantities parametrized in log-scale with O(1) priors; raw counts never enter the likelihood. |
| A3 | **Variance decomposition** is the headline scientific question. Handled by fixing the prior on input rate $r_i$ from Sperry/HR/WF and freely estimating a residual learning-efficiency term $\alpha_i$. Proportion of outcome variance attributable to input = $\text{Var}(\log r)/(\text{Var}(\log r) + \text{Var}(\log \alpha))$. |
| A4 | Word intercepts grouped by lexical class, with a class-specific frequency slope (so noun bias ≈ class mean × class frequency slope). |
| B1–B4 | English CDI:WS Wordbank + CHILDES frequencies as they stand; no ODS. |
| B5 | Start time $s_j$ per word with prior centered at ~4.5 mo. |
| C1 | Fit Wordbank is the acceptance test. |
| C3 | Subsample (500 children × 200 items) first; scale after it fits. |
| D3 | Age-rate-change $\delta$ included as a testable parameter with $\delta = 0$ as the null. |

---

## The model

### Quantities (all absolute, all on log scale when they enter $\eta$)

- $r_i$ — child $i$'s CDS input rate in tokens/hour.
- $\alpha_i$ — child $i$'s learning-efficiency multiplier (dimensionless; $\log \alpha_i$ centered at 0).
- $p_j$ — corpus probability of word $j$ (CHILDES, fixed data).
- $H$ — waking hours/month (fixed constant, 365).
- $a_i$ — age at CDI assessment in months.
- $s_j$ — start time (months of age) at which accumulation for word $j$ begins.
- $\psi_j$ — log of the expected cumulative tokens of word $j$ needed for acquisition (i.e. the threshold on the log-token scale).
- $\delta$ — age rate-change exponent ($\delta = 0$: constant rate; $\delta > 0$: older children accumulate faster).

### Linear predictor

For the constant-rate case ($\delta = 0$):
$$\eta_{ij} = \log r_i + \log \alpha_i + \log p_j + \log H + \log(a_i - s_j) - \psi_j$$

With age rate-change:
$$\eta_{ij} = \log r_i + \log \alpha_i + \log p_j + \log H + (1+\delta)\log\!\left(\frac{a_i - s_j}{a_0}\right) - \psi_j$$

where $a_0 = 20$ months is a fixed reference age to keep $\delta$'s interpretation clean.

Acquisition probability (probit link):
$$P(y_{ij} = 1 \mid \ldots) = \Phi(\eta_{ij}/\sigma)$$

with $\sigma$ a noise parameter (I will start with $\sigma$ fixed to 1 so the log-token scale is the scale; if the model needs more flex we can free it).

**Edge case.** When $a_i \le s_j$ the log is undefined. We handle this by setting $P(y_{ij}=1) = \epsilon$ (a small floor, e.g. $10^{-4}$) and not using those child-word cells in the log-age term — i.e., conditioning on $a_i > s_j$.

### Priors (all in units we can reason about)

**Input rate.** Informed entirely by the Sperry/HR/WF pooled distribution, *not* estimated from Wordbank:
$$\log r_i \sim N(\mu_r, \sigma_r), \quad \mu_r = \log 1198, \quad \sigma_r = \text{sd}(\log r)_\text{Sperry}.$$

This is the key move for A3: $r_i$'s variance is *pinned by external data*, so any residual child-level variance in outcomes must flow through $\log \alpha_i$.

**Learning efficiency.**
$$\log \alpha_i \sim N(0, \sigma_\alpha), \quad \sigma_\alpha \sim \text{HalfNormal}(0, 1).$$

$\sigma_\alpha$ is the quantity we most want to recover.

**Word thresholds by lexical class.**
$$\psi_j \mid \text{class}(j) = c \sim N(\mu_c, \tau_c),$$

with $\mu_c \sim N(\log 1000, 2)$, $\tau_c \sim \text{HalfNormal}(0, 1)$.

So the mean threshold is ~1000 tokens ± an order of magnitude per class. Noun bias falls out of $\mu_\text{noun} < \mu_\text{verb}$.

Optional extension (probably important given your 2021 findings): a class-specific frequency slope $\beta_c$ so word-$j$ difficulty is $\psi_j - \beta_{c(j)} \log p_j$. If $\beta_c < 1$, frequency matters less than one would naïvely expect from token accumulation. See Question 1 below.

**Start time.**
$$s_j \sim \text{Normal}(4.5, 2), \text{ truncated to } [0, 15].$$

**Age rate-change.**
$$\delta \sim N(0, 0.5).$$

Null hypothesis $\delta = 0$ is inside the prior; prior allows but doesn't favor rate growth.

### Likelihood

$$y_{ij} \sim \text{Bernoulli}(\Phi(\eta_{ij})).$$

---

## Identifiability sketch (A3)

Only $\log r_i + \log \alpha_i$ enters the likelihood at the child level. Without external information, these are not separately identifiable. We identify them by:

1. Placing a **strong informative prior** on $\log r_i$ from external data (Sperry/HR/WF). The variance $\sigma_r$ is known, not estimated.
2. Putting a vague prior on $\log \alpha_i$ so its variance $\sigma_\alpha$ is free.

Posterior samples of $\sigma_\alpha$ then directly inform:
- Expected proportion of vocabulary variance due to learning efficiency: $\sigma_\alpha^2 / (\sigma_r^2 + \sigma_\alpha^2)$.
- 95% CrI on that proportion.

**Sensitivity analysis.** Refit with $\sigma_r$ set to 0.5×, 1×, 2× the empirical value — this tests how much our answer depends on the input-rate variance estimate. If the answer is robust, great; if not, we know the analysis is bottlenecked on better input measurement (which is itself publishable).

---

## Implementation plan

### Phase 1 — write the Stan model and simulate
1. Code `models/log_irt.stan`.
2. Write a simulator in R that generates $(r_i, \alpha_i, p_j, \psi_j, s_j, y_{ij})$ from the model with known parameters.
3. Recover parameters. Confirm:
   - $\psi_j$ recovered in log-token units.
   - $\sigma_\alpha$ recovered.
   - $\delta$ recovered.
   - Variance-decomposition proportion recovered.

### Phase 2 — fit subset of Wordbank
1. 500 children × 200 items (stratified by age bin and lexical class).
2. Run 4 chains, 1000 warm-up + 1000 sample.
3. Diagnostics: Rhat, divergences, ESS. Posterior predictive check against the held-out 2-year-old curves.

### Phase 3 — scale to full Wordbank
1. All 5,492 children × 680 items.
2. If too slow, use `cmdstanr` + reduce_sum threading or switch to Laplace / VI as a warm start.

### Phase 4 — planned comparisons
1. Constant-rate ($\delta = 0$) vs. age-rate-change ($\delta$ free).
2. With vs. without class-specific frequency slope $\beta_c$.
3. With vs. without word-specific start times $s_j$.

### Phase 5 — simulations and the Observable app
1. Counterfactual simulator: plug in $p(r)$ and $p(p)$, get predicted vocab-size distributions by age.
2. Observable app: sliders for ($\mu_r$, $\sigma_r$, $\sigma_\alpha$, $\psi$-scale, $\delta$), live plots of growth curves and word-learning curves.

---

## Remaining questions (just two)

**Q1. Class-specific frequency slopes?** Your 2021 paper showed frequency explains less than it "should" — nouns are more frequency-sensitive than function words. Concretely, should the model be:

- **(simpler)** $\eta$ contains $\log p_j$ with unit coefficient (the "one token counts as one" assumption), and class only shifts $\psi_j$?
- **(richer)** $\eta$ contains $\beta_{c(j)} \log p_j$ with a class-specific slope? Then $\beta_c < 1$ means frequency matters sub-linearly for that class.

Richer is more informative and directly tests the "pure accumulator" assumption. I lean toward richer, but it adds 4 parameters and potentially weakens frequency-threshold identifiability.

**Q2. $s_j$ — per word, per class, or global?** Three levels:

- Global $s$ (one number): simplest.
- Per-class $s_c$: captures that e.g. function words start late.
- Per-word $s_j$ with class-mean hierarchy: richest, but with 680 items this could soak up a lot of degrees of freedom and trade off with $\psi_j$.

M&P used per-word. I lean toward **per-class** start times in the baseline (5 numbers) and make per-word a reachable extension.

---

Once you confirm Q1/Q2, I'll write the Stan file and the simulator.
