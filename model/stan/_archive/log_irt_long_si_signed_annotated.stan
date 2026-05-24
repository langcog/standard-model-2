// ============================================================================
// log_irt_long_si_signed_annotated.stan
//
// READING COPY — do NOT use this file for fits. The working file is
// `log_irt_long_si_signed.stan`. This annotated copy is byte-identical in
// model semantics but with extensive inline comments explaining every
// design choice and every line of nontrivial math. Use it to read the
// model carefully and ask questions. Keep the two in sync if you change
// the working model.
// ============================================================================
//
// PURPOSE OF THIS MODEL
// ---------------------
// This is the **headline** longitudinal Stan model in standard_model_2 —
// the one labeled `long_no_freq_slopes_si_signed` in the 5-variant family
// build, used for the main results on English and Norwegian Wordbank.
// It is the most-fully-specified model: every parameter that earlier
// variants in the family pin (σ_α, κ_pop, σ_ζ, σ_s) is freed here.
//
// Quick statement of what's modeled:
//   - Each observation y_n ∈ {0,1} is whether child i produced word j on
//     admin a (the CDI form they filled out at age admin_age[a]).
//   - We model logit P(y_n = 1) = child intercept + log H + (population
//     growth rate + per-child slope deviation) * log(age - s - s_i[i])
//     - word difficulty.
//   - This is a 1PL Rasch IRT model (per-item difficulty delta_j; all
//     items share unit discrimination) with a power-law time predictor,
//     fit hierarchically.
//   - (Earlier versions of this file included 2PL discrimination
//     lambda_j and a per-class frequency slope beta_c[cc[j]] * log p_j.
//     Both were pinned off in every active variant of the 5-panel
//     family build, and were removed on 2026-05-21 for readability.
//     If you ever want them back, see the git history of this file.)
//
// WHAT'S NEW IN si_signed vs the base model (log_irt_long.stan)
// -------------------------------------------------------------
// Two changes, both about how per-child onset s_i is parameterized.
// The base model used HALF-NORMAL s_i ≥ 0 (delays only) and sampled
// σ_ζ and σ_s as direct free parameters. This caused two distinct
// sampling pathologies that we fixed here:
//
//   1) SIGNED-NORMAL s_i with sum-to-zero centering.
//      Old: s_i ~ Normal_+(0, σ_s), so s_i ≥ 0 always. Interpretation
//           was "child i is delayed by s_i months relative to population
//           start s." A child with s_i = 0 was exactly at the
//           population reference; a child could only be late, never early.
//      New: s_i ~ Normal(0, σ_s), unconstrained sign, but with the
//           constraint that the empirical mean of s_i across kids is
//           exactly zero (enforced in transformed parameters by
//           subtracting the empirical mean).
//      Why: Under the half-normal version, the population mean of s_i
//           was σ_s * sqrt(2/π) ≈ 0.8 σ_s, which moves with σ_s. The
//           population's "true" onset is s + E[s_i], so changing σ_s
//           changes the population onset. Sampler exploited this by
//           trading off between σ_s (large) and s (small) along a
//           diagonal compensation ridge in the joint (σ_s, s) posterior.
//           Result: poor mixing on s and σ_s individually, even though
//           the data identified σ_xi tightly.
//      Fix: with signed s_i + sum-to-zero, E[s_i] = 0 exactly regardless
//           of σ_s. So σ_s only controls the spread of s_i around the
//           population, and s carries the population onset alone. The
//           (σ_s, s) ridge is broken.
//      Interpretation change: s_i is now a developmental DEVIATION from
//           the population reference, signed. s_i < 0 means "ahead of
//           the population reference" (their developmental clock at
//           calendar age t reads more than t-s; equivalently log_age is
//           larger). s_i > 0 means "behind." This is analogous to how
//           ζ_i is a deviation from the population growth rate, not a
//           one-sided modifier.
//
//   2) (σ_total, p_zeta) REPARAMETERIZATION of (σ_ζ, σ_s).
//      Old: σ_ζ and σ_s sampled directly with independent half-normal
//           priors.
//      New: Sample
//             σ_total ~ HalfNormal(0, 5)             (broad)
//             p_zeta  ~ Beta(2, 2)                   (weakly uniform)
//           and back-transform in transformed parameters:
//             σ_ζ = σ_total * sqrt(p_zeta)
//             σ_s = σ_total * sqrt(1 - p_zeta)
//           This is a one-to-one map; σ_ζ and σ_s retain their original
//           interpretation and are reported in posterior summaries.
//      Why: σ_ζ and σ_s both inflate the variance of η across (i, t)
//           observations (each per-child random effect "dilates" the
//           per-child predicted trajectory in slightly different ways).
//           The data tightly identify σ_ξ² + per-child-variance, but
//           the SPLIT between σ_ζ² and σ_s² is weak — they trade off
//           against each other along a diagonal in the joint posterior.
//           Direct sampling of (σ_ζ, σ_s) thus has a diagonal ridge.
//           The reparam rotates coordinates 45°: σ_total moves along
//           the ridge (data-identified), p_zeta moves across the ridge
//           (poorly identified, but prior-dominated by Beta(2, 2)).
//           Sampler now sees axis-aligned posterior geometry: easy to
//           explore.
//      Why this works: HMC efficiency depends on the posterior having
//           "round" level sets, not elongated ones. The reparam doesn't
//           reduce posterior uncertainty (you can't; the data doesn't
//           identify p_zeta well) — it just changes coordinates so the
//           uncertainty isn't expressed as a diagonal correlation that
//           the sampler has to chase.
//      Dimensional caveat: σ_ζ scales log-time (multiplicative on
//           log_age), σ_s scales age-months (additive inside log()).
//           Treating σ_total² = σ_ζ² + σ_s² is a SAMPLING choice — it
//           gives Stan an axis-aligned coordinate system — not a claim
//           that σ_ζ and σ_s are dimensionally commensurable. What the
//           data identifies is the marginal posterior on each.
//
// LINEAR PREDICTOR (single equation form, identical to base model)
// ---------------------------------------------------------------
//   η_{i,a,j} = λ_j * [ ξ_i
//                       + β_c[cc[j]] * log p_j
//                       + log H
//                       + (1 + δ + ζ_i) * log((age_a - s - s_i[i])/a0)
//                       - δ_j ]
//
// Each term:
//   ξ_i        = child intercept on the latent ability scale.
//                Includes both log r_i (per-child input rate) and log α_i
//                (per-child efficiency). Drawn from N(μ_r, σ_ξ²) with
//                σ_ξ² = σ_r² + σ_α² (so σ_α absorbs the variance not
//                attributable to input rate). When we report π_α =
//                σ_α² / (σ_α² + σ_r²) in §19, this is the share.
//   β_c[cc[j]] = per-class log-p slope. Pinned at 1 in the "unit
//                accumulator" baseline (DEFAULT_PRIORS, beta_c_prior_sd
//                = 0.001). The `no_freq_*` family of variants pins it at
//                0 instead, dropping the explicit log-p frequency term
//                (item difficulty then absorbs the frequency-related
//                difficulty into δ_j alone).
//   log H      = log(365 hr/mo). A unit conversion: when we say "input
//                rate r_i is in tokens/hr," the ability accumulator wants
//                tokens, so we multiply r_i by hours-per-month-of-life.
//                This is a constant offset on the ability scale.
//   δ          = population shift away from the unit-accumulator
//                exponent. The full exponent on log_age is 1 + δ. When
//                δ = 0 we have pure linear-in-log-time accumulation
//                (κ_pop = 1). When δ > 0 the population is super-linear
//                (κ_pop = 1 + δ > 1). The headline result of the paper
//                is δ ≈ 9, i.e. κ_pop ≈ 10.
//   ζ_i        = per-child deviation from population growth rate.
//                Child i has growth rate (1 + δ + ζ_i) = κ_pop + ζ_i =
//                κ_i. Drawn from a bivariate N with ξ_i via Cholesky and
//                LKJ correlation. Centered to mean zero across kids.
//   s          = global population onset offset. Time before s has
//                no learning (or undefined log). Pinned at 6 mo by
//                DEFAULT_PRIORS in this work; previously was 0.5 mo
//                — the change was made so that signed s_i deviations
//                around s = 6 yield positive effective onsets ~6 ± 2 mo,
//                which is interpretively cleaner.
//   s_i[i]     = per-child trajectory phase, signed and sum-to-zero
//                centered (E[s_i] = 0 exactly). σ_s controls the
//                spread.
//   a0         = age normalization constant (the "reference age" at
//                which we say "0 effective experience"). Set to 19 mo
//                by DEFAULT_PRIORS — roughly the median CDI:WS age. The
//                log(age/a0) form means at a0 the time term is 0 and the
//                per-child effects only see ξ_i and δ_j.
//   δ_j        = word-level difficulty, hierarchical by lexical class.
//                δ_j ~ N(μ_c[class], τ_c[class]).
//   λ_j        = word-level discrimination (2PL). Hierarchical log-normal
//                via λ_j = exp(σ_λ * z_j) where z_j ~ N(0,1). When
//                σ_λ_prior_sd is pinned tiny, λ_j ≈ 1 and we have 1PL Rasch.
//
// SAMPLING NOTES
// --------------
// reduce_sum threading: the per-observation likelihood is parallelized
// across threads via `reduce_sum(partial_sum_lpmf, ...)`. Each thread
// computes η locally for its slice of y_slice, then evaluates
// bernoulli_logit_lpmf on that slice. Linear speedup with
// STAN_NUM_THREADS up to ~16 threads per chain on physical cores.
//
// Per-observation log-likelihood for PSIS-LOO is computed in `generated
// quantities` after sampling, in a single thread, using the same η
// formula (without reduce_sum). The unused `eta` local in `model` would
// be redundant; we structure both the partial_sum_lpmf body and the
// generated_quantities loop to compute η identically, so log_lik from
// LOO matches the in-sample likelihood exactly.
//
// ============================================================================


// ----------------------------------------------------------------------------
// FUNCTIONS BLOCK
// ----------------------------------------------------------------------------
//
// Stan's reduce_sum lets us parallelize a sum-of-likelihood-terms across
// threads. The convention: write a function that takes a SLICE of the
// outcome vector y_slice (a chunk of consecutive observations) and
// returns the SUM of log-likelihoods over that slice. Stan then chops
// the full y vector into chunks, evaluates each chunk on a separate
// thread, and adds them up.
//
// Why this is fast: each thread reads its slice of y, computes η for
// that slice, and adds bernoulli_logit_lpmf(slice | η) to target. With
// N = 1.1M observations and 16 threads, each thread processes ~70K obs.
// Threads only share read-only references to parameters; no contention.
//
// The function signature takes:
//   y_slice: chunk of outcomes (passed by reduce_sum)
//   start, end: 1-indexed bounds of the slice in the original y
//   aa, jj, admin_to_child: integer index arrays (same as in data)
//   admin_age, s, a0, time_baseline, delta, log_H: scalars + age vector
//   delta_j: word-level difficulty
//   xi, zeta, s_i: child-level random effects
functions {
  real partial_sum_lpmf(array[] int y_slice,
                        int start, int end,
                        array[] int aa, array[] int jj,
                        array[] int admin_to_child,
                        vector admin_age, real s, real a0,
                        real time_baseline, real delta,
                        real log_H,
                        vector delta_j,
                        vector xi, vector zeta, vector s_i) {

    int n_slice = end - start + 1;
    vector[n_slice] eta_slice;

    // Loop over the slice. For each observation:
    //   - Look up which admin (a), word (j), child (ch) it's on
    //   - Compute effective age ae = max(t - s - s_i[ch], 0.01)
    //     The fmax with 0.01 prevents log(<=0); when (t - s - s_i) is
    //     negative or tiny (e.g., a child with s_i > t - s) we floor it
    //     so log_age stays well-defined. The Stan sampler can still try
    //     such samples (e.g., during early warmup) without crashing.
    //   - Compute log_age relative to the reference age a0
    //   - Compute slope = time_baseline + delta + zeta[ch]
    //     time_baseline is 1 by default; this is the "(1 + delta + ζ_i)"
    //     coefficient on log_age. (When time_baseline = 0 in m0/m1
    //     variants, we are testing models without time at all.)
    //   - Assemble eta = xi[ch] + log_H + slope*log_age - delta_j[j]
    //     This is the "ability minus difficulty" term, directly used
    //     as the logit (no lambda multiplier; this is a 1PL Rasch model).
    for (i in 1:n_slice) {
      int n = start + i - 1;
      int ch = admin_to_child[aa[n]];                          // child index for this obs
      real ae = fmax(admin_age[aa[n]] - s - s_i[ch], 0.01);    // effective age (>= 0.01)
      real log_age_n = log(ae / a0);                            // log-time on reference scale
      real slope_n = time_baseline + delta + zeta[ch];          // (1 + delta + ζ_i) = κ_i
      eta_slice[i] = xi[ch] + log_H + slope_n * log_age_n - delta_j[jj[n]];
    }

    // Sum log Bernoulli pmf over the slice. Returns to Stan's target += .
    return bernoulli_logit_lpmf(y_slice | eta_slice);
  }
}


// ----------------------------------------------------------------------------
// DATA BLOCK
// ----------------------------------------------------------------------------
//
// All quantities here are SUPPLIED BY R when stan_data is built. None
// of these are inferred — they're either observed (y, admin_age,
// log_p) or pinned (a0, log_H, prior hyperparameters).
data {
  // ---- Counts (used as array dimensions) ----
  int<lower=1> N;                          // observations
  int<lower=1> A;                          // CDI administrations (one per child×age)
  int<lower=1> I;                          // children
  int<lower=1> J;                          // distinct words (CDI items)
  int<lower=1> C;                          // lexical classes (e.g., nouns/verbs/adj/...)

  // ---- Index arrays ----
  // aa[n] tells us which admin observation n is part of. jj[n] tells
  // us which word. admin_to_child[a] tells us which child filled out
  // admin a. cc[j] tells us which class word j belongs to.
  array[N] int<lower=1, upper=A> aa;
  array[N] int<lower=1, upper=J> jj;
  array[A] int<lower=1, upper=I> admin_to_child;
  array[J] int<lower=1, upper=C> cc;

  // ---- Outcomes ----
  // y[n] = 1 if child produced (or comprehended) word j on admin a, else 0.
  array[N] int<lower=0, upper=1> y;

  // ---- Continuous data ----
  vector[A] admin_age;                      // age in months for each admin
  vector[J] log_p;                          // log token-frequency for each word

  // ---- Pinned constants ----
  real log_H;                               // log of hours/month = log(365) ≈ 5.9
  real<lower=0> a0;                         // reference age in months (typically 19)

  // ---- Population input-rate prior (externally pinned) ----
  // We don't infer the population input rate from the data; it's pinned
  // externally from observational studies (Sperry et al., Hart & Risley, etc.)
  // because the CDI data alone can't distinguish input variance from
  // efficiency variance (the variance partition is not data-identified
  // without an external pin). σ_r = 0.534 is the central Western-CDS
  // empirical estimate. See input_estimation/README.md for derivation.
  real mu_r;
  real<lower=0> sigma_r;

  // ---- Hierarchical word-difficulty prior (per-class) ----
  // delta_j ~ N(mu_c[class], tau_c[class]). mu_c is itself N(mu_mu_c, sigma_mu_c).
  real mu_mu_c;
  real<lower=0> sigma_mu_c;

  // ---- Variant toggle priors ----
  // Each of these widths can be tiny (e.g. 0.001) to "pin" the
  // corresponding random effect at its prior mean (effectively
  // disabling it). variant_hyperpriors() in helpers.R sets these
  // depending on which model variant we're fitting (demo_pure,
  // demo_alpha, ..., no_freq_slopes_si_signed).
  real s_prior_mean;                        // population onset s (typically 6 mo)
  real<lower=0> s_prior_sd;
  real delta_prior_mean;                    // population κ_pop shift (typically 0)
  real<lower=0> delta_prior_sd;
  real<lower=0> sigma_lambda_prior_sd;      // UNUSED after 2026-05-21 cleanup; kept for bundle compat
  // sigma_zeta_prior_sd and sigma_s_prior_sd are NOT used in this file —
  // the (σ_total, p_zeta) reparam replaces them with sigma_total and
  // p_zeta. Kept in the data block for compatibility with bundles built
  // for the base log_irt_long.stan file.
  real<lower=0> sigma_zeta_prior_sd;
  real<lower=0> sigma_alpha_prior_sd;
  real<lower=0> sigma_s_prior_sd;
  real<lower=0> beta_c_prior_sd;            // UNUSED after cleanup; kept for bundle compat
  real beta_c_prior_mean;                   // UNUSED after cleanup; kept for bundle compat
  real time_baseline;                       // typically 1; m0/m1 variants set 0
}


// ----------------------------------------------------------------------------
// PARAMETERS BLOCK
// ----------------------------------------------------------------------------
//
// These are the unknowns Stan samples. Every parameter here gets a prior
// in the model block (or in transformed parameters, derived from these
// via fixed transformations).
//
// Non-centered parameterization: where possible we sample raw N(0,1)
// quantities here and apply the scale-shift in transformed parameters.
// This avoids funnel geometry where a tightly-shrunken random effect
// has a hard-to-sample joint posterior with its scale parameter.
parameters {
  // ---- Bivariate (xi, zeta) per child via Cholesky non-centered ----
  // z_child is a 2xI matrix of std-normal raw values. The Cholesky
  // factor L_child encodes the correlation between xi and zeta. The
  // actual (xi, zeta) values are built in transformed parameters by
  // applying L_child * z_child plus the appropriate scales (sigma_xi
  // for the xi axis, sigma_zeta for the zeta axis).
  //
  // Note: only L_child[2,1] is meaningfully free (the lower triangle
  // of a 2x2 Cholesky of a correlation matrix). The other entries are
  // pinned by the LKJ prior + the constraint that diagonals come from
  // sqrt(1 - subdiagonal²).
  matrix[2, I] z_child;
  cholesky_factor_corr[2] L_child;

  // sigma_alpha: SD of the per-child efficiency component of xi.
  // sigma_xi is derived in transformed parameters from sqrt(sigma_r² + sigma_alpha²).
  // When sigma_alpha_prior_sd is tiny (demo_pure variant), this is pinned
  // near 0 and σ_ξ ≈ σ_r (all between-child variance comes from input rate).
  real<lower=0> sigma_alpha;

  // Word-level: FLAT (non-class-stratified) hyperprior on δ_j.
  // Previously this was class-stratified (vector[C] mu_c, vector<lower=0>[C] tau_c
  // with δ_j ~ N(μ_c[class], τ_c[class])); collapsed to a single global
  // hyperprior on 2026-05-21 because (a) the headline params don't
  // depend on which class shrinkage δ_j gets, (b) it's fewer parameters
  // and simpler to explain in the model build. cc[] and C remain in
  // data because β_c uses them (when un-pinned).
  vector[J] delta_j_raw;                    // non-centered, scaled by tau_delta in TP
  real mu_delta;                            // global word-difficulty mean
  real<lower=0> tau_delta;                  // global word-difficulty SD

  // Population structural parameters.
  real<lower=0, upper=15> s;                // onset offset (months). [0, 15] is a hard bound.
  real delta;                               // shift from unit-accumulator exponent
                                            // (κ_pop = 1 + delta)

  // (2PL discrimination lambda_j and per-class frequency slope beta_c
  // were removed 2026-05-21; both were pinned off in every active
  // variant. The active model is 1PL Rasch with no frequency term.)

  // ---- s_i: per-child trajectory phase, SIGNED ----
  // s_i_raw is a vector of unconstrained reals (one per child). It's
  // centered to mean zero in transformed parameters so the empirical
  // mean of s_i across kids is exactly 0. This breaks the (s_i mean,
  // population s) compensation ridge that plagued the half-normal
  // version. Sampling is on s_i_raw directly; the sum-to-zero is a
  // deterministic transform.
  vector[I] s_i_raw;

  // ---- REPARAMETERIZED (σ_ζ, σ_s) channels ----
  // sigma_total = sqrt(sigma_zeta² + sigma_s²): the "total" between-kid
  // scale on the joint (zeta, s) axis. Data identifies this tightly.
  real<lower=0> sigma_total;

  // p_zeta = sigma_zeta² / sigma_total² ∈ (0, 1): the fraction of total
  // between-kid scale carried by zeta (vs. s). Weakly identified by the
  // data; Beta(2,2) prior centers it near 0.5 with broad spread.
  real<lower=0, upper=1> p_zeta;
}


// ----------------------------------------------------------------------------
// TRANSFORMED PARAMETERS BLOCK
// ----------------------------------------------------------------------------
//
// Deterministic transforms of `parameters`. Stan computes these every
// leapfrog step (along with the gradient) so they participate in HMC.
// Anything you can compute as a function of `parameters` and `data`
// goes here, especially the non-centered scale-shift steps and any
// quantities you want reported in posterior summaries.
transformed parameters {
  // ---- Derive sigma_xi from input variance + efficiency variance ----
  // sigma_xi² = sigma_r² + sigma_alpha². This is the total scale on the
  // per-child intercept ξ_i. The xi axis of the (xi, zeta) bivariate
  // gets this scale; the zeta axis gets sigma_zeta below.
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));

  // ---- BACK-TRANSFORM the reparameterized variance channels ----
  // We sampled sigma_total and p_zeta (axis-aligned coordinates). Now
  // recover the human-interpretable sigma_zeta and sigma_s for use in
  // the model block (in the (xi, zeta) bivariate scaling) and for
  // posterior summaries (so we can report sigma_zeta and sigma_s
  // posteriors, not sigma_total and p_zeta).
  //
  // Reversal of the reparam:
  //   sigma_zeta² = sigma_total² * p_zeta       => sigma_zeta = sigma_total * sqrt(p_zeta)
  //   sigma_s²    = sigma_total² * (1 - p_zeta) => sigma_s    = sigma_total * sqrt(1 - p_zeta)
  // These are the same sigma_zeta and sigma_s you'd find in the base
  // model (log_irt_long.stan); only the SAMPLING coordinates differ.
  real<lower=0> sigma_zeta = sigma_total * sqrt(p_zeta);
  real<lower=0> sigma_s    = sigma_total * sqrt(1 - p_zeta);

  // ---- Build the (xi, zeta) random effects via Cholesky ----
  // The 2x2 covariance matrix Σ for (ξ_i, ζ_i) has:
  //   Σ[1,1] = sigma_xi²
  //   Σ[2,2] = sigma_zeta²
  //   Σ[1,2] = Σ[2,1] = sigma_xi * sigma_zeta * rho
  // where rho is the correlation between xi and zeta (drawn from LKJ).
  //
  // The Cholesky decomposition is Σ = L_scaled * L_scaled' where
  //   L_scaled[1,1] = sigma_xi
  //   L_scaled[1,2] = 0                                        (upper tri)
  //   L_scaled[2,1] = sigma_zeta * L_child[2,1]                (= sigma_zeta * rho)
  //   L_scaled[2,2] = sigma_zeta * L_child[2,2]                (= sigma_zeta * sqrt(1-rho²))
  // where L_child is the Cholesky factor of the *correlation* matrix
  // (sampled directly from LKJ).
  //
  // Then we get (xi_raw, zeta_raw) by multiplying L_scaled times the
  // std-normal raw z_child. Result is a 2xI matrix; we transpose so
  // child_effs is Ix2 with column 1 = xi-axis, column 2 = zeta-axis.
  matrix[I, 2] child_effs;
  {
    matrix[2, 2] L_scaled;
    L_scaled[1, 1] = sigma_xi * L_child[1, 1];
    L_scaled[1, 2] = 0;
    L_scaled[2, 1] = sigma_zeta * L_child[2, 1];
    L_scaled[2, 2] = sigma_zeta * L_child[2, 2];
    child_effs = (L_scaled * z_child)';
  }

  // ---- Sum-to-zero centering ----
  // Three random-effect vectors get centered so their empirical means
  // across kids are exactly zero. Without centering, each random-effect
  // mean partly absorbs the corresponding population-level fixed effect:
  //   mean(xi) and mu_r both contribute to E[xi]
  //   mean(zeta) and (1 + delta) both contribute to E[slope]
  //   mean(s_i) and s both contribute to E[onset]
  // Centering pins each fixed effect to carry its full intended role.
  vector[I] xi   = mu_r + child_effs[, 1] - mean(child_effs[, 1]);
  vector[I] zeta = child_effs[, 2] - mean(child_effs[, 2]);
  vector[I] s_i  = s_i_raw - mean(s_i_raw);

  // ---- Word-level: build delta_j from raw + hyperparams ----
  // Flat hyperprior: delta_j[j] = mu_delta + tau_delta * delta_j_raw[j].
  // All words share one mean and one SD; no class-stratified pooling.
  //
  // We briefly tried sum-to-zero centering here (delta_j = mu_delta +
  // delta_j_unc - mean(delta_j_unc)) to match the kid-level random
  // effects, hoping to break the small mu_delta-vs-globals soft
  // correlations seen in pilot v1. Result: chain agreement on rho_xi_zeta
  // and sigma_alpha got worse, because mean(delta_j_unc) introduces a
  // J-wide coupling between every delta_j_raw[k] that NUTS finds harder
  // to traverse. Reverted. The mu_delta soft correlations (~0.3-0.4)
  // remain but are mild enough to live with.
  vector[J] delta_j = mu_delta + tau_delta * delta_j_raw;
}


// ----------------------------------------------------------------------------
// MODEL BLOCK
// ----------------------------------------------------------------------------
//
// Where we add log-priors and log-likelihood to `target`. Stan then
// samples from the resulting joint posterior via HMC/NUTS.
model {
  // ---- Priors on the reparameterized variance channels ----
  // sigma_total has a broad half-normal. At N = 1.1M obs the data
  // dominates this prior easily; we just need to keep the parameter
  // bounded above zero and finite. SD=5 is much wider than any
  // plausible posterior.
  sigma_total ~ normal(0, 5);

  // p_zeta has a weakly-uniform Beta(2,2) prior. Beta(2,2) is the
  // simplest conjugate prior that's smooth on (0,1), centered at 0.5,
  // and not bunched at either endpoint. Beta(1,1) would be exactly
  // uniform but allows the sampler to drift to endpoints; Beta(2,2)
  // gently pulls away from 0 and 1.
  p_zeta ~ beta(2, 2);

  // ---- Bivariate (xi, zeta) prior ----
  // Non-centered: z_child is std-normal, L_child is LKJ(2).
  // LKJ(2) prior on the 2x2 correlation matrix means corr(xi, zeta) is
  // pulled gently toward 0 (uniform on [-1, 1] would be LKJ(1); larger
  // shape parameter shrinks toward identity).
  to_vector(z_child) ~ std_normal();
  L_child           ~ lkj_corr_cholesky(2);

  // sigma_alpha: half-normal with width from data block. When this is
  // tight (demo_pure variant), sigma_alpha is pinned near 0.
  sigma_alpha ~ normal(0, sigma_alpha_prior_sd);

  // ---- Word-level priors ----
  delta_j_raw ~ std_normal();               // non-centered
  mu_delta  ~ normal(mu_mu_c, sigma_mu_c);  // hyperprior on global mean difficulty
                                            // (reuses mu_mu_c, sigma_mu_c numerics from
                                            // the old class-stratified prior)
  tau_delta ~ normal(0, 1);                 // hyperprior on global difficulty SD

  // ---- Population structural priors ----
  s     ~ normal(s_prior_mean, s_prior_sd);                  // onset: N(6, 0.05) post-2026-05-21
  delta ~ normal(delta_prior_mean, delta_prior_sd);          // κ_pop shift

  // ---- Per-child onset s_i ----
  // We SAMPLE s_i_raw (the uncentered version) with mean zero and SD
  // sigma_s. Then in transformed parameters, s_i is the sum-to-zero-
  // centered version. The prior is conceptually on s_i_raw, not on
  // the centered s_i: the centering is a deterministic post-hoc step
  // that removes one degree of freedom from the sample (the empirical
  // mean), but since the empirical mean was already prior-mean-zero,
  // this doesn't bias anything.
  s_i_raw ~ normal(0, sigma_s);

  // ---- Likelihood, parallelized via reduce_sum ----
  // The full N-observation Bernoulli likelihood is chunked across
  // threads. Each thread evaluates partial_sum_lpmf on its slice.
  // Stan adds up the partial sums into target.
  target += reduce_sum(partial_sum_lpmf, y, 1,
                       aa, jj, admin_to_child,
                       admin_age, s, a0,
                       time_baseline, delta, log_H,
                       delta_j,
                       xi, zeta, s_i);
}


// ----------------------------------------------------------------------------
// GENERATED QUANTITIES BLOCK
// ----------------------------------------------------------------------------
//
// Computed AFTER sampling, per draw. These don't affect the posterior;
// they're derived quantities we want in the posterior summary or for
// downstream use (e.g., PSIS-LOO).
generated quantities {
  // π_α: the share of intercept-variance attributable to per-child
  // efficiency (sigma_alpha²) vs. input rate (sigma_r²). Headline result
  // of the paper: π_α ≈ 0.91 in English, 0.94 in Norwegian.
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));

  // Correlation between ξ_i and ζ_i (the off-diagonal of the LKJ-derived
  // correlation matrix). Typically positive (~0.35 in our fits): kids
  // with higher baselines tend to also have steeper growth rates.
  real rho_xi_zeta = L_child[2, 1];

  // ---- Per-observation log-likelihood for PSIS-LOO ----
  // Single-threaded; mirrors the η computation in partial_sum_lpmf above.
  // For large N (e.g., Norwegian N = 2.18M), this matrix is huge —
  // sherlock/extract_loo_thinned.R thins to 2000 draws before computing
  // PSIS-LOO in chunks.
  vector[N] log_lik;
  {
    // Pre-compute log_age once (small saving but cleaner)
    vector[N] log_age;
    for (n in 1:N) {
      int ch_n = admin_to_child[aa[n]];
      log_age[n] = log(fmax(admin_age[aa[n]] - s - s_i[ch_n], 0.01) / a0);
    }
    for (n in 1:N) {
      int ch = admin_to_child[aa[n]];
      real slope_n = time_baseline + delta + zeta[ch];
      real eta_n = xi[ch] + log_H + slope_n * log_age[n] - delta_j[jj[n]];
      log_lik[n] = bernoulli_logit_lpmf(y[n] | eta_n);
    }
  }

  // ---- Per-child "efficiency" projection ----
  // xi_i = log r_i + log α_i, and we've pinned σ_r externally. So we
  // can decompose each child's xi into the input-explained part
  // (proportional to σ_r²/σ_ξ²) and the efficiency-explained part
  // (proportional to σ_α²/σ_ξ²). This projection gives the
  // efficiency-attributed component:
  //   log_alpha_mean[i] = π_α * (xi[i] - mu_r)
  // (since π_α = σ_α²/σ_ξ²).
  //
  // We report this so per-child efficiency posteriors can be examined.
  // Note: this is the POSTERIOR MEAN projection; individual log_α_i
  // is not separately identified from log_r_i in this model (we only
  // observe xi).
  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
}
