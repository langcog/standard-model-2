// ============================================================================
// log_irt_long_annotated.stan
//
// READING COPY — do NOT use this file for fits. The working file is
// `log_irt_long.stan`. This annotated copy is byte-identical in model
// semantics but with extensive inline comments explaining every design
// choice and every line of nontrivial math. Use it to read the model
// carefully and ask questions. Keep the two in sync if you change the
// working model.
// ============================================================================
//
// PURPOSE OF THIS MODEL
// ---------------------
// This is the **headline** longitudinal Stan model in standard_model_2
// (post-2026-05-22 cleanup). It is the file used to fit `M_best`
// = `no_freq_slopes` = the model labeled "α + ζ + δ" in the slide
// deck, on English Wordbank (I=500, J=671, N=1.1M) and Norwegian
// (I=500, J=674, N=2.06M).
//
// Quick statement of what's modeled:
//   - Each observation y_n ∈ {0,1} is whether child i produced word j
//     on admin a (the CDI form they filled out at age admin_age[a]).
//   - We model logit P(y_n = 1) =
//       child intercept ξ_i
//       + log H (unit conversion: tokens/month from tokens/hour)
//       + (1 + δ + ζ_i) · log((age_a − s − s_i[i]) / a_0)
//       − word difficulty δ_j .
//   - This is a 1PL Rasch IRT model (per-item difficulty δ_j; all items
//     share unit discrimination) with a power-law time predictor, fit
//     hierarchically.
//
// WHY THIS FILE STILL HAS s, s_i, AND sigma_s
// -------------------------------------------
// The current headline M_best regime pins `s` tight at 0 (via
// DEFAULT_PRIORS: s_prior_mean = 0, s_prior_sd = 0.001) and pins
// `sigma_s` tight near 0 (via sigma_s_prior_sd = 0.001), which in
// turn pins `s_i ≈ 0` for every child. So in production runs of this
// file, the `s` and `s_i` terms are essentially zero and the linear
// predictor reduces to:
//   η = ξ_i + log H + (1 + δ + ζ_i) · log(age_a / a_0) − δ_j .
//
// Why keep them in the Stan file then? Two reasons:
//   1. Backward compatibility with stan_data bundles built for older
//      variants. The data block must still accept the s_prior_*,
//      sigma_s_prior_sd, beta_c_prior_*, etc. fields even though
//      they're now pinned tight.
//   2. Variants that DO want a free s or s_i (e.g., the
//      `free_s_no_freq_slopes` robustness check in §20 of
//      experiments.md) can run this file with looser priors — no
//      separate Stan file needed.
// See model/R/helpers.R variant_hyperpriors() for the variant grammar.
//
// HISTORY (what was removed in the 2026-05-21 cleanup)
// ----------------------------------------------------
// Earlier versions of this file included two now-removed terms in the
// linear predictor:
//
//   - 2PL per-word discrimination λ_j (hierarchical log-normal):
//     η = λ_j × [...] (multiplicative front of the bracket).
//     Pinned at 1 in every active variant (sigma_lambda_prior_sd ≈
//     0.001). Removed to simplify the predictor and the Stan code.
//   - Per-class log-p frequency slope β_c[cc[j]] × log p_j inside the
//     bracket. Pinned at 0 in `no_freq_*` variants and at 1 in the
//     unit-accumulator default. Removed for the same reason; per-word
//     frequency now lives entirely inside δ_j.
//
// Per-class hierarchical δ_j prior (μ_c[class], τ_c[class]) was also
// flattened to a global (μ_δ, τ_δ) on the same date — see the comment
// in the parameters block.
//
// The data block retains the deprecated fields (sigma_lambda_prior_sd,
// beta_c_prior_*, log_p, cc, C) so bundles built for older variants
// remain accepted. They are unused in the current model.
//
// LINEAR PREDICTOR (single equation form, post-cleanup)
// -----------------------------------------------------
//   η_{i,a,j} = ξ_i + log H + (1 + δ + ζ_i) · log((age_a − s − s_i[i]) / a_0) − δ_j
//
// Each term:
//   ξ_i        = child intercept on the latent ability scale.
//                Conceptually ξ_i = log r_i + log α_i. Drawn from
//                N(μ_r, σ_ξ²) with σ_ξ² = σ_r² + σ_α². When we report
//                π_α = σ_α² / (σ_α² + σ_r²) we mean: of the kid-level
//                intercept variance, this share is "efficiency" vs
//                "input rate." σ_r is pinned externally; σ_α is freely
//                estimated.
//   log H      = log(365 hr/mo). Unit conversion: r_i is in tokens/hr;
//                multiplying by hours-per-month converts to a per-month
//                effective experience rate. Constant offset.
//   δ          = population shift away from the unit-accumulator
//                exponent. The full exponent on log_age is 1 + δ.
//                δ = 0 is pure linear-in-log-time (κ_pop = 1, McMurray
//                pure-accumulator). δ > 0 is super-linear (Hidaka
//                rate-change). EN posterior median ≈ 10.3 (κ_pop ≈ 11.3).
//   ζ_i        = per-child deviation from population growth rate. The
//                kid's effective exponent is κ_i = 1 + δ + ζ_i.
//                Drawn jointly with ξ_i from a bivariate Normal via
//                Cholesky + LKJ correlation. Centered to mean zero
//                across kids.
//   s          = global population onset offset. Time before s has
//                no learning (the fmax(t − s − s_i, 0.01) floor
//                kicks in if anyone tries it). Pinned at 0 in M_best
//                (was 0.5 historically, was briefly pinned at 6 during
//                the si_signed experiment, finally settled at 0 after
//                that experiment was backed out — see §23 of
//                experiments.md).
//   s_i[i]     = per-child trajectory phase, half-normal in this file
//                (s_i ≥ 0). Pinned at ~0 in M_best via tight σ_s.
//                Note: the experimental signed-s_i variant lived in
//                log_irt_long_si_signed.stan (now in _archive/) where
//                s_i could be negative; that file is no longer used.
//   a_0        = age normalization constant ("reference age"). Set
//                from the bundle as the dataset's median admin age:
//                ~18 mo for the EN longitudinal bundle, ~19 mo for
//                Norwegian. The log(age/a_0) form means at a_0 the
//                time term is 0 and per-child effects only see ξ_i
//                and δ_j.
//   δ_j        = word-level difficulty. Flat hierarchical prior
//                δ_j ~ N(μ_δ, τ_δ), shared across lexical classes
//                (post-2026-05-21 flattening; see comment in
//                parameters block).
//
// SAMPLING NOTES
// --------------
// reduce_sum threading: the per-observation likelihood is parallelized
// across threads via `reduce_sum(partial_sum_lpmf, ...)`. Each thread
// computes η locally for its slice of y_slice, then evaluates
// bernoulli_logit_lpmf on that slice. Linear speedup with
// STAN_THREADS_PER_CHAIN up to ~16 threads per chain on physical
// cores. (Beyond ~16 the TBB scheduler overhead starts to dominate
// the gain.)
//
// Per-observation log-likelihood for PSIS-LOO is computed in
// `generated quantities` after sampling, in a single thread, using
// the same η formula (without reduce_sum). The η computation in
// partial_sum_lpmf and the GQ loop are structured identically, so
// log_lik from LOO matches the in-sample likelihood exactly.
//
// IDENTIFIABILITY NOTES
// ---------------------
// 1. (μ_r, mean(ξ)). The data only sees ξ_i, not its mean separately
//    from μ_r. We pin μ_r externally and sum-to-zero center ξ_i so
//    mean(ξ) = μ_r exactly.
// 2. (δ, mean(ζ)). Same story for the population slope: the
//    likelihood only sees (1 + δ + ζ_i), so the mean of ζ could
//    absorb part of δ. Sum-to-zero centering of ζ pins δ to carry
//    the full population slope.
// 3. (σ_r, σ_α). NOT identified separately from data — both add to
//    σ_ξ². σ_r is pinned externally from observational studies
//    (Sperry/HR/W&F), letting σ_α be identified as σ_ξ² − σ_r²
//    residual. See §27 in experiments.md and `outputs/input_rate_table.md`.
// 4. (s, δ). When s is near zero and fmax doesn't activate, the
//    (s, δ) gradient correlation is moderate; with s pinned at 0 and
//    σ_s ≈ 0 the ridge is effectively switched off.
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
// N = 1.1M observations and 16 threads, each thread processes ~70K
// obs. Threads only share read-only references to parameters; no
// contention.
//
// The function signature takes:
//   y_slice         : chunk of outcomes (passed by reduce_sum)
//   start, end      : 1-indexed bounds of the slice in the original y
//   aa, jj          : observation -> admin, observation -> word indices
//   admin_to_child  : admin -> child index
//   admin_age, s, a0, time_baseline, delta, log_H : scalars + age vec
//   delta_j         : word-level difficulty
//   xi, zeta, s_i   : child-level random effects
functions {
  real partial_sum_lpmf(array[] int y_slice,
                        int start, int end,
                        // observation indices
                        array[] int aa, array[] int jj,
                        array[] int admin_to_child,
                        // global / scalar
                        vector admin_age, real s, real a0,
                        real time_baseline, real delta,
                        real log_H,
                        // per-item / per-child
                        vector delta_j,
                        vector xi, vector zeta, vector s_i) {

    int n_slice = end - start + 1;
    vector[n_slice] eta_slice;

    // Loop over the slice. For each observation:
    //   - Look up which admin (a), word (j), child (ch) it's on.
    //   - Compute effective age ae = max(t − s − s_i[ch], 0.01). The
    //     fmax with 0.01 prevents log(<= 0); when (t − s − s_i) is
    //     negative or tiny we floor it so log_age stays well-defined.
    //     The sampler can still propose such samples during warmup
    //     without crashing. In M_best with s ≈ 0 and s_i ≈ 0 this is
    //     basically a no-op.
    //   - Compute log_age relative to the reference age a_0.
    //   - Compute slope = time_baseline + delta + zeta[ch], the
    //     "(1 + δ + ζ_i)" coefficient on log_age. time_baseline is 1
    //     by default; m0/m1 variants set it to 0 to test "no time
    //     term" reductions.
    //   - Assemble eta = xi[ch] + log_H + slope * log_age − delta_j[j].
    //     This is the "ability minus difficulty" logit, with no λ_j
    //     multiplier (1PL Rasch).
    for (i in 1:n_slice) {
      int n = start + i - 1;
      int ch = admin_to_child[aa[n]];                       // child index
      real ae = fmax(admin_age[aa[n]] - s - s_i[ch], 0.01); // effective age (≥ 0.01)
      real log_age_n = log(ae / a0);                         // log-time on reference scale
      real slope_n = time_baseline + delta + zeta[ch];       // (1 + δ + ζ_i) = κ_i
      eta_slice[i] = xi[ch] + log_H + slope_n * log_age_n - delta_j[jj[n]];
    }

    // Sum log Bernoulli pmf over the slice. reduce_sum will add this
    // to Stan's `target` via the parent target += statement.
    return bernoulli_logit_lpmf(y_slice | eta_slice);
  }
}


// ----------------------------------------------------------------------------
// DATA BLOCK
// ----------------------------------------------------------------------------
//
// All quantities here are SUPPLIED BY R when stan_data is built. None
// of these are inferred — they're either observed (y, admin_age),
// pinned (a0, log_H, prior hyperparameters), or deprecated-but-kept
// (log_p, cc, C, beta_c_*, sigma_lambda_prior_sd) for backward
// compatibility with older bundles.
data {
  // ---- Counts (used as array dimensions) ----
  int<lower=1> N;                       // observations
  int<lower=1> A;                       // CDI administrations (one per child × age)
  int<lower=1> I;                       // children
  int<lower=1> J;                       // distinct words (CDI items)
  int<lower=1> C;                       // lexical classes — UNUSED after 2026-05-21
                                        // cleanup; kept so bundles built for
                                        // earlier versions remain accepted.

  // ---- Index arrays ----
  // aa[n]            tells us which admin observation n is part of
  // jj[n]            tells us which word
  // admin_to_child[a] tells us which child filled out admin a
  // cc[j]            tells us which class word j belongs to (UNUSED here;
  //                  was used by the removed β_c log-p frequency term)
  array[N] int<lower=1, upper=A> aa;
  array[N] int<lower=1, upper=J> jj;
  array[A] int<lower=1, upper=I> admin_to_child;
  array[J] int<lower=1, upper=C> cc;    // UNUSED (see above)

  // ---- Outcomes ----
  // y[n] = 1 if child produced word j on admin a, else 0.
  array[N] int<lower=0, upper=1> y;

  // ---- Continuous data ----
  vector[A] admin_age;                  // age in months for each admin
  vector[J] log_p;                      // UNUSED after cleanup; kept for bundle compat

  // ---- Pinned constants ----
  real log_H;                           // log of hours/month ≈ log(365) ≈ 5.9
  real<lower=0> a0;                     // reference age in months
                                        // (~18-19 for the longitudinal bundles)

  // ---- Population input-rate prior (externally pinned) ----
  // We don't infer the population input rate from the CDI data; it's
  // pinned externally from observational studies (Sperry/HR/W&F).
  // Reason: σ_r and σ_α are not separately data-identified — both
  // contribute to σ_ξ. Pinning σ_r lets σ_α be recovered as the
  // residual. σ_r = 0.534 is the central Western-CDS empirical pin
  // (Sperry within-pool). See input_estimation/README.md and
  // outputs/input_rate_table.md for the sensitivity sweep + sources.
  real mu_r;
  real<lower=0> sigma_r;

  // ---- Word-difficulty hyperprior (flat, global) ----
  // mu_mu_c and sigma_mu_c are the hyperprior on μ_δ (the global
  // mean of δ_j). The names are historical: they were originally
  // for the class-stratified mu_c[class] vector. After the 2026-05-21
  // flattening, only the global μ_δ uses these numerics; the names
  // stick for backward bundle compatibility.
  real mu_mu_c;
  real<lower=0> sigma_mu_c;

  // ---- Variant toggle priors ----
  // Each of these widths can be made tiny (e.g., 0.001) to "pin" the
  // corresponding random effect at its prior mean (effectively
  // disabling it), or set to a regular value (~1) to let the data
  // inform the parameter freely. variant_hyperpriors() in helpers.R
  // sets these depending on which model variant is being fit
  // (demo_pure, demo_alpha, demo_kappa, no_freq_slopes, etc.).
  real s_prior_mean;                    // population onset s
                                        // (M_best: 0; was 6 during si_signed era)
  real<lower=0> s_prior_sd;             // tight (0.001) pins s
  real delta_prior_mean;                // population κ_pop shift (typically 0)
  real<lower=0> delta_prior_sd;         // wide (0.5) lets data find δ
  real<lower=0> sigma_lambda_prior_sd;  // UNUSED after cleanup; kept for bundle compat
  real<lower=0> sigma_zeta_prior_sd;    // tight pins ζ off; ~1 frees per-kid slopes
  real<lower=0> sigma_alpha_prior_sd;   // tight pins σ_α near 0 (pure-accumulator demos);
                                        //   ~1 in M_best frees efficiency variation
  real<lower=0> sigma_s_prior_sd;       // per-child onset s_i. Tight (0.001) pins
                                        //   s_i ≈ 0 (M_best default).
  real<lower=0> beta_c_prior_sd;        // UNUSED; kept for bundle compat
  real beta_c_prior_mean;               // UNUSED; kept for bundle compat
  real time_baseline;                   // 1 = unit-rate accumulator; 0 = no time
                                        //   term (m0 variant)
}


// ----------------------------------------------------------------------------
// PARAMETERS BLOCK
// ----------------------------------------------------------------------------
//
// These are the unknowns Stan samples. Every parameter here gets a
// prior in the model block (or in transformed parameters, derived
// from these via fixed transformations).
//
// Non-centered parameterization: where possible we sample raw N(0,1)
// quantities here and apply the scale-shift in transformed parameters.
// This avoids funnel geometry where a tightly-shrunken random effect
// has a hard-to-sample joint posterior with its scale parameter.
parameters {
  // ---- Bivariate (ξ_i, ζ_i) per child via Cholesky non-centered ----
  // z_child is a 2×I matrix of std-normal raw values. The Cholesky
  // factor L_child encodes the correlation between ξ and ζ. The
  // actual (ξ_i, ζ_i) values are built in transformed parameters by
  // applying a scaled L (with σ_ξ on the ξ axis and σ_ζ on the ζ axis)
  // times z_child.
  //
  // Note: only L_child[2,1] is meaningfully free (the off-diagonal of
  // a 2×2 Cholesky of a correlation matrix). The other entries are
  // pinned by the LKJ prior + the constraint that diagonals come
  // from sqrt(1 − subdiagonal²).
  matrix[2, I] z_child;

  // sigma_child is a length-2 vector with two roles. (And this is
  // confusing — see the comment block in transformed parameters.)
  //   sigma_child[1] : technically the "ξ-axis SD" but we DON'T use
  //                    it as such; sigma_xi is derived deterministically
  //                    from σ_r² + σ_α² instead. So sigma_child[1] just
  //                    gets a harmless N(0, 1) prior and is otherwise
  //                    ignored. (We keep it in the vector because Stan
  //                    likes the type to be vector[2] for the bivariate
  //                    construction.)
  //   sigma_child[2] : the actual sigma_zeta = SD of ζ_i. The variant
  //                    grammar routes the sigma_zeta_prior_sd here so
  //                    that `baseline` (tight prior) pins ζ off and
  //                    `slopes` (~1) frees it.
  vector<lower=0>[2] sigma_child;

  // L_child is the Cholesky factor of the (ξ, ζ) correlation matrix.
  // LKJ(2) prior in the model block pulls the correlation gently
  // toward 0.
  cholesky_factor_corr[2] L_child;

  // σ_α: SD of the per-child efficiency component of ξ. When
  // sigma_alpha_prior_sd is tight (demo_pure variant), this is pinned
  // near 0 and σ_ξ ≈ σ_r (all between-child variance comes from
  // input rate). When loose (M_best), σ_α is data-identified as the
  // residual variance σ_ξ² − σ_r².
  real<lower=0> sigma_alpha;

  // ---- Word-level: FLAT (global) hyperprior on δ_j ----
  // Previously this was class-stratified (vector[C] mu_c, vector<lower=0>[C]
  // tau_c with δ_j ~ N(μ_c[class], τ_c[class])); collapsed to a single
  // global hyperprior on 2026-05-21 because:
  //   (a) the headline params (σ_α, κ_pop, σ_ζ, σ_s) don't depend on
  //       which class shrinkage δ_j gets,
  //   (b) it's fewer parameters and simpler to explain in the model build.
  // cc[] and C remain in the data block since β_c (now removed) used
  // them, and bundles built for earlier versions still pass them.
  vector[J] delta_j_raw;                // non-centered std-normal; scaled
                                        // by tau_delta in TP
  real mu_delta;                        // global word-difficulty mean
  real<lower=0> tau_delta;              // global word-difficulty SD

  // ---- Population structural parameters ----
  real<lower=0, upper=15> s;            // onset offset (months). [0, 15] hard bound.
                                        // M_best pins this at 0 via s_prior_sd = 0.001.
  real delta;                           // shift from unit-accumulator exponent.
                                        // The full kappa_pop is 1 + delta.

  // (2PL discrimination λ_j and per-class frequency slope β_c were
  // removed 2026-05-21. Both were pinned off in every active variant.
  // Removing them simplified the linear predictor to:
  //   eta = xi + log_H + (1 + delta + zeta) * log_age − delta_j .)

  // ---- s_i: per-child onset, HALF-NORMAL (s_i ≥ 0) ----
  // CENTERED parameterization: we sample s_i directly (not via raw +
  // scale). The lower=0 bound makes this a half-normal.
  //
  // Non-centered was tried (s_i = sigma_s × s_i_raw) but produced a
  // bimodal posterior: one mode had σ_s ~ 36 with σ_α ~ 0.9. The
  // multiplicative factorization lets the model trade kid-onset
  // variation for kid-efficiency variation across distant
  // configurations. Centered ties each s_i directly to the data,
  // which empirically gives single-mode posteriors. When
  // sigma_s_prior_sd is tight (~0.001), σ_s collapses to ~0 and
  // s_i ~ 0 — the M_best default.
  //
  // Note: this file uses HALF-NORMAL s_i ≥ 0. The experimental
  // SIGNED-NORMAL variant lived in
  // _archive/log_irt_long_si_signed.stan (now retired — see §23 of
  // experiments.md for why).
  vector<lower=0>[I] s_i;
  real<lower=0> sigma_s;
}


// ----------------------------------------------------------------------------
// TRANSFORMED PARAMETERS BLOCK
// ----------------------------------------------------------------------------
//
// Deterministic transforms of `parameters`. Stan computes these every
// leapfrog step (along with the gradient) so they participate in HMC.
// Anything you can compute as a function of `parameters` and `data`
// goes here, especially:
//   - Non-centered scale-shift steps for random effects
//   - Quantities you want reported in posterior summaries
//   - Constraints expressed via deterministic relationships
transformed parameters {
  // ---- Derive sigma_xi from input variance + efficiency variance ----
  // σ_ξ² = σ_r² + σ_α². This is the total scale on the per-child
  // intercept ξ_i. The ξ axis of the (ξ, ζ) bivariate gets this
  // scale; the ζ axis gets sigma_child[2] (= σ_ζ).
  //
  // We DON'T declare σ_ξ as a parameter and try to enforce the
  // constraint via a prior — instead we compute it deterministically
  // here so the constraint is exact, and the upstream sampler only
  // sees σ_α (and σ_r as data). This is the "σ_r externally pinned"
  // convention.
  real<lower=0> sigma_xi = sqrt(square(sigma_r) + square(sigma_alpha));

  // ---- Build the (ξ, ζ) random effects via Cholesky ----
  // The 2×2 covariance matrix Σ for (ξ_i, ζ_i) has:
  //   Σ[1,1] = σ_ξ²
  //   Σ[2,2] = σ_ζ²
  //   Σ[1,2] = Σ[2,1] = σ_ξ · σ_ζ · ρ
  // where ρ is the correlation between ξ and ζ (drawn from LKJ).
  //
  // The Cholesky decomposition is Σ = L_scaled · L_scaled' where
  //   L_scaled[1,1] = σ_ξ
  //   L_scaled[1,2] = 0                                  (upper triangle)
  //   L_scaled[2,1] = σ_ζ · L_child[2,1]                 (= σ_ζ · ρ)
  //   L_scaled[2,2] = σ_ζ · L_child[2,2]                 (= σ_ζ · sqrt(1 − ρ²))
  // where L_child is the Cholesky factor of the *correlation* matrix
  // (sampled directly from LKJ).
  //
  // Then we get (ξ_raw, ζ_raw) by multiplying L_scaled times the
  // std-normal raw z_child. Result is a 2×I matrix; we transpose so
  // child_effs is I×2 with column 1 = ξ-axis, column 2 = ζ-axis.
  //
  // Subtlety: sigma_child[1] is NOT used here — we use σ_ξ instead
  // (the deterministic version that depends on σ_r and σ_α). So
  // sigma_child[1] is effectively a "phantom" parameter that exists
  // only because the vector[2] type requires both entries. Its
  // normal(0, 1) prior in the model block is harmless (it doesn't
  // couple to anything else).
  matrix[I, 2] child_effs;
  {
    matrix[2, 2] L_scaled;
    L_scaled[1, 1] = sigma_xi * L_child[1, 1];
    L_scaled[1, 2] = 0;
    L_scaled[2, 1] = sigma_child[2] * L_child[2, 1];
    L_scaled[2, 2] = sigma_child[2] * L_child[2, 2];
    child_effs = (L_scaled * z_child)';
  }

  // ---- Sum-to-zero centering on ξ and ζ ----
  // Without centering, each random-effect mean partly absorbs the
  // corresponding population-level fixed effect:
  //   mean(ξ) and μ_r both contribute to E[ξ] from the model's POV.
  //   mean(ζ) and (1 + δ) both contribute to E[slope].
  // Centering pins each fixed effect to carry its full intended
  // role: μ_r is the population intercept, (1 + δ) is the population
  // slope, and the random effects deviate symmetrically around them.
  //
  // (Note: s_i is NOT sum-to-zero centered in this file because it's
  // half-normal — its mean is positive by construction (≈ σ_s ·
  // sqrt(2/π) ≈ 0.8 σ_s). That positivity is one of the structural
  // reasons we backed away from the signed-s_i / sum-to-zero
  // approach — see §23 of experiments.md.)
  vector[I] xi   = mu_r + child_effs[, 1] - mean(child_effs[, 1]);
  vector[I] zeta = child_effs[, 2] - mean(child_effs[, 2]);
  real<lower=0> sigma_zeta = sigma_child[2];

  // ---- Word-level: build δ_j from raw + hyperparams ----
  // Flat hyperprior: δ_j[j] = μ_δ + τ_δ · δ_j_raw[j].
  // All words share one mean and one SD; no class-stratified pooling.
  //
  // We briefly tried sum-to-zero centering here on 2026-05-21 to
  // match the kid-level random effects, hoping to break the small
  // (μ_δ vs μ_r, μ_δ vs δ) soft correlations seen in pilot v1.
  // Result: chain agreement on rho_xi_zeta and σ_α got WORSE, because
  // mean(delta_j_unc) introduces a J-wide coupling between every
  // delta_j_raw[k] that NUTS finds harder to traverse than the small
  // (~0.3-0.4) soft correlations it was meant to fix. Reverted.
  vector[J] delta_j = mu_delta + tau_delta * delta_j_raw;
}


// ----------------------------------------------------------------------------
// MODEL BLOCK
// ----------------------------------------------------------------------------
//
// Where we add log-priors and log-likelihood to `target`. Stan then
// samples from the resulting joint posterior via HMC/NUTS.
model {
  // ---- (ξ, ζ) bivariate prior, non-centered ----
  // z_child is std-normal raw values; the scale-shift is in TP.
  to_vector(z_child) ~ std_normal();

  // sigma_child[1]: harmless N(0, 1) prior. This parameter is
  // declared because the bivariate construction needs vector[2], but
  // it's not actually used as the ξ-axis scale — σ_ξ is derived
  // deterministically from σ_r² + σ_α². Keeping a proper prior here
  // ensures the sampler doesn't drift to extreme values that would
  // slow down adaptation.
  sigma_child[1] ~ normal(0, 1);

  // sigma_child[2] = σ_ζ : the actual slopes-toggle parameter.
  // Variant grammar routes the sigma_zeta_prior_sd here so:
  //   - `baseline` (tight prior ~0.001) pins ζ off
  //   - `slopes`   (~1)                  frees per-kid slopes
  sigma_child[2] ~ normal(0, sigma_zeta_prior_sd);

  // L_child : LKJ(2) prior. Concentration parameter 2 gently pulls
  // ρ(ξ, ζ) toward 0 (uniform would be LKJ(1); larger η shrinks
  // toward identity).
  L_child ~ lkj_corr_cholesky(2);

  // σ_α : half-normal with width from data block. When
  // sigma_alpha_prior_sd is tight (demo_pure variant), σ_α is pinned
  // near 0.
  sigma_alpha ~ normal(0, sigma_alpha_prior_sd);

  // ---- Word-level priors ----
  delta_j_raw ~ std_normal();                // non-centered raw
  mu_delta  ~ normal(mu_mu_c, sigma_mu_c);   // hyperprior on global mean difficulty
                                             // (reuses mu_mu_c / sigma_mu_c numerics
                                             // from the old class-stratified prior)
  tau_delta ~ normal(0, 1);                  // hyperprior on global difficulty SD

  // ---- Population structural priors ----
  s     ~ normal(s_prior_mean, s_prior_sd);            // onset: N(0, 0.001) in M_best
  delta ~ normal(delta_prior_mean, delta_prior_sd);    // κ_pop shift

  // ---- Per-child onset s_i (half-normal via lower=0 on s_i) ----
  // sigma_s_prior_sd ~ 0.001 (M_best default) effectively pins s_i at 0;
  // sd = 2 frees per-child onset variation. Centered parameterization:
  // data directly informs each s_i. (Non-centered was tried and
  // produced multimodal posteriors; centered is single-mode.)
  s_i     ~ normal(0, sigma_s);
  sigma_s ~ normal(0, sigma_s_prior_sd);

  // ---- Likelihood, parallelized via reduce_sum ----
  // The full N-observation Bernoulli likelihood is chunked across
  // threads. Each thread evaluates partial_sum_lpmf on its slice.
  // grainsize = 1 lets Stan's TBB scheduler auto-tune the slice size
  // (it'll pick something like N / (4 × n_threads) under the hood).
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
  // π_α: the share of intercept variance attributable to per-child
  // efficiency (σ_α²) vs input rate (σ_r²). Headline result of the
  // paper: π_α ≈ 0.92 (EN) and 0.94 (NO).
  real pi_alpha = square(sigma_alpha) / (square(sigma_alpha) + square(sigma_r));

  // Correlation between ξ_i and ζ_i (the off-diagonal of the
  // LKJ-derived correlation matrix). EN M_best posterior: ρ ≈ −0.09
  // [−0.18, −0.01] (essentially zero); NO posterior: ρ ≈ −0.13.
  // (The pre-excision si_signed-era ρ was +0.35 in some fits;
  // changing s and excising s_i flipped the structure.)
  real rho_xi_zeta = L_child[2, 1];  // since L_child is 2×2 Cholesky

  // ---- Per-observation log-likelihood for PSIS-LOO ----
  // Single-threaded; mirrors the η computation in partial_sum_lpmf
  // above. For large N (e.g., Norwegian N = 2.06M), this matrix is
  // huge — sherlock/extract_loo_thinned.R thins to 2000 draws before
  // computing PSIS-LOO in column chunks. After the disk-full crash
  // experience (§24 of experiments.md) this is mostly a safety net;
  // we don't routinely compute LOO for the headline I=500 fits.
  vector[N] log_lik;
  {
    // Pre-compute log_age once per observation (small saving, cleaner)
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
  // ξ_i = log r_i + log α_i, and we've pinned σ_r externally. So we
  // can decompose each child's ξ into the input-explained part
  // (proportional to σ_r² / σ_ξ²) and the efficiency-explained part
  // (proportional to σ_α² / σ_ξ²). This projection gives the
  // efficiency-attributed component:
  //   log_alpha_mean[i] = π_α · (ξ_i − μ_r)
  // (since π_α = σ_α² / σ_ξ²).
  //
  // Reported so per-child efficiency posteriors can be examined
  // downstream. Note: this is the POSTERIOR MEAN projection;
  // individual log_α_i is not separately identified from log_r_i
  // in this model (we only observe ξ_i).
  vector[I] log_alpha_mean;
  for (i in 1:I) {
    log_alpha_mean[i] = (square(sigma_alpha) / square(sigma_xi)) * (xi[i] - mu_r);
  }
}
