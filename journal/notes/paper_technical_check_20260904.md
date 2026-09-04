# Technical check of the PNAS draft (sm2_input.pdf, compiled 2026-09-04)

Statements in the draft that do not match what the code does, with the draft's
margin line numbers. Text/flow/argument deliberately not reviewed. Verified
against: `fits/joint_io_proc_english_count_subset_data.rds`,
`model/stan/log_irt_long_proc_bilingual.stan`, `model/stan/log_irt_long.stan`,
`studies/glmer_ladder/02_fit_one.R`, `studies/cross_sectional_demographics/00_build.R`
+ README, `paper/build_fig_io_cache.R`, `model/scripts/data_check_io_panels.R`,
`journal/experiments.md` #43–45.

## Errors (wrong as written)

1. **Eq. 1 sign (l. 291).** `P = 1/(1 + e^{θ − δ})` decreases in θ; the model is
   `logit⁻¹(θ_i(t) − δ_j) = 1/(1 + e^{−(θ − δ)})`.
2. **"Four include recordings from day-long audio recordings using LENA" (l. 665–667).**
   Three LENA datasets (SEEDLingS, AM2018, FMW2013; 808 recordings, instrument 2) plus
   BabyView head-cam video (5,739 videos, instrument 1). Four datasets have input; three
   are LENA.
3. **Fig. 3 caption, "non-parametric GAMLSS smoother" (l. 533–534).** It is a loess
   smoother (`geom_smooth(method = "loess", span = 1)`), one per dataset.
4. **Fig. 4 caption, "95% confidence intervals" (l. 626–627).** Posterior 2.5/97.5%
   quantiles (credible intervals) from the draws; "confidence" is right for Fig. 2 only.
5. **Fig. 1 caption, "processing trials n" (l. 276–277) and the plate label.** The
   processing observations are per-administration mean log RTs from Peekbank
   (`1_d_sub.Rds`; 1,374 rows for 350 children, 2.1–6.6 per child), not trials. The TikZ
   plate now says "session n"; the caption should too.
6. **Fig. 1 caption, "fitted population parameters from the English (American) data"
   (l. 279–280).** Panel B uses the joint io/proc fit on the six English datasets
   (κ = 11.3, σ_ξ = 2.02, σ_κ = 4.40, a₀ = 21 mo, log H = 5.9), not the Wordbank English
   fits; the reference child is 0.5 SD below mean efficiency and the two comparison lines
   add +1 SD of ξ and +1 SD of κ.
7. **Fig. 2 caption (l. 428–430).** "Blue bands" are now shaded in the trait colours
   (blue efficiency, vermillion acceleration); "error bars show 95% CIs" holds for the
   cross-sectional points only (longitudinal diamonds carry no CI in A/B; they do in C).
8. **"cross-language meta-analytic estimates … (which largely agree)" (l. 555–557).**
   Three of four cells agree; sex → acceleration does not: cross-sectional 0.35
   [0.24, 0.45] vs longitudinal −0.11 [−0.47, 0.24]. The study README also records that
   cross-sectional acceleration effects ran 1.5–2.3× the longitudinal ones in validation
   (an upper bound); the paper's "sex differences manifest primarily in efficiency"
   rests on the longitudinal arm.
9. **"a small but significant effect on acceleration" (l. 685–686; also abstract l. 20–22,
   intro l. 233–235).** With the 95% intervals now used everywhere: input → acceleration
   0.72 [0.02, 1.42] in the joint model, 0.76 [−0.11, 1.55] fit alone; P(> 0) = 0.98 /
   0.96. "Positive but imprecise" (as the Discussion already says at l. 756–759).

## Descriptions that are misleading or incomplete

10. **"a moderate effect [of input] on efficiency" (l. 684–685).** Not an estimated
    coupling: ξ_i = μ_ξ + d_i + …, coefficient fixed at 1 (footnote †), so the per-SD
    "effect" plotted is σ_r itself (0.35 [0.31, 0.40]) — the SD of log input rate across
    children. Say so where it is read as a finding.
11. **Longitudinal demographics, "We fit the accelerating accumulator (without
    measurement models) to data from each language" (l. 541–547).** These are lme4
    `glmer` fits (model D_log: `produces ~ log_age + (1 + log_age | child) + (1 | item)`,
    `nAGQ = 0`, a₀ = median age of the dataset), fit per *dataset* (Thal and Marchman
    separately for English; Norwegian; Japanese), not Bayesian and not per language.
    BLUPs are then regressed on the predictor per dataset (`lm(xi ~ sex)`,
    `lm(zeta ~ scale(ed_years))`); Fig. 2A/B pools the two English datasets by inverse
    variance, and the longitudinal meta uses k = 4 (sex) / k = 2 (maternal ed) units.
12. **Cross-sectional demographics, "apply the same model to recover population
    parameters" (l. 547–552).** Per language, one administration per child, WG + WS
    (incl. short forms), item-level Rasch GLMM `produces ~ predictor × log(age/a₀) +
    (1 | item) + (1 | child)` (`glmer`, `nAGQ = 0`, a₀ = median age); main effect =
    efficiency, interaction = acceleration. Undocumented inclusion rules: monolingual-TD
    datasets only (bilingual + clinical datasets dropped), ≥ 300 children per language
    displayed (28 languages sex, 12 maternal ed), English (British) and French (Québec)
    excluded, French (French) dropped from maternal ed (skewed education distribution).
    Maternal education = caregiver education mapped to years, z-scored (caption correct).
13. **Measurement models (l. 350–352).** Input SD is per instrument (LENA vs head-cam
    transcripts), study means per study; LWL has one residual SD, study intercepts τ_s and
    one global log-age slope ψ (= −0.37). "general standard deviations" undersells this.
14. **Count administrations not described (Model section).** 95 FMW2013 ELENA
    administrations enter as vocabulary sumscores through a Poisson-binomial
    (normal-approximation) likelihood over the form's items, not item by item
    (1,520 item-level + 95 count administrations).
15. **Units in Fig. 2C / Fig. 4 and the text's "moderate vs small" comparison
    (l. 683–686).** Efficiency effects are logits per SD; acceleration effects are κ units
    (per unit log age). At 30 months a κ effect converts by × log(30/21) = 0.36, so
    input → acceleration 0.72 ≈ 0.26 logits. State the units in both captions; the
    processing effects are also sign-flipped to "faster processing → +".
16. **σ_r band provenance (l. 719–722, Fig. 5).** The plausible band [0.32, 0.58] runs
    from the direct fit's own lower 95% bound on observed recordings (a lower bound on
    population σ_r) to the literature upper value; anchors 0.35 (direct estimate), 0.44
    (population re-anchor, the dashed line), 0.58 (literature). Not purely "estimates in
    the published literature".
17. **"strong congruence" with 4–7% (l. 726–731).** At σ_r = 0.44: English 5.7%,
    Norwegian 3.1%; across the band English 3.0–9.9%, Norwegian 1.7–5.5%. Overlap,
    not congruence, for Norwegian.
18. **"focused exclusively on monolingual, typically-developing children" (l. 779–781).**
    True for the Wordbank arms (bilingual + clinical datasets excluded). Nothing in the
    io-proc bundle filters on term status or late-talker status: AM2018 is the
    full-term + preterm cohort (ref. 50's title) and FM2012 includes late talkers
    (ref. 17's title). Verify against the source datasets before keeping the sentence.

## Checked and correct

N = 413 children, six datasets, CDI administrations 6–32 mo (LWL 8–30 mo, count
administrations 17–25 mo); five datasets with multiple LWL sessions per child; log RT
declining linearly in log age (ψ); no age term in the input model; Eq. 3–4 match the Stan
structural block (independent residuals); imputed model ξ_i = μ_r + log r_i + log α_i with
σ_r pinned, share = σ_r²/σ_ξ², three anchors per language; separate-vs-joint arms on the
same bundle with near-identical estimates; repo URL; raw io/proc source files and the
Wordbank pull via wordbankr are in the repository.

## Methods / SI checklist — things the pipeline does that the draft does not describe

Facts verified against the scripts, bundles and summaries on 2026-09-04.

### Data
1. **Dataset table.** None in the draft. `paper/cache/table1_datasets.csv` has it: Wordbank
   longitudinal Thal 653 / Smith 316 / Marchman 2,136 / Norwegian 1,630 / Japanese 187
   children (admins, age ranges); the six io/proc datasets (children, CDI admins, LWL
   sessions, recordings, ages); 31 cross-sectional languages with N. `si_io_data_table.csv`
   is the io/proc detail (413 kids, 1,520 item-level admins, 1,374 LWL sessions, 219 with
   input).
2. **Wordbank extraction.** wordbankr 2.0.0 via Redivis, Wordbank v2.0; forms WG + WS
   incl. short forms (`FORM_KEEP`); child key = `dataset::study_internal_id` (links WG and
   WS administrations); uni_lemma cross-form item harmonisation; crater/jump QC on each
   child's trajectory (`QC_REL_TOL` 0.25, `QC_PEAK` 0.10, `QC_DROP` 0.05; removed 90
   Marchman and 50 Norwegian admins); longitudinal = ≥ 2 distinct ages; monolingual-TD:
   "Bilingual" dataset origins and the clinical datasets (Edgin, Byers) excluded.
3. **English longitudinal pooling.** Thal + Smith + Marchman form one English dataset for
   the imputed fits (I = 3,105, a₀ = 18 mo); Norwegian I = 1,630 (a₀ = 23); Japanese
   (187) enters only the demographic BLUP arm.
4. **io/proc CDI.** Stanford lab item files (TL2, TL3, TLO + ELENA sumscores; TL3 WS only),
   BabyView webCDI (item keys mapped), SEEDLingS CDI; 681 WS items; four lexical classes
   for the difficulty hierarchy (101 / 318 / 97 / 165 items); 95 ELENA administrations as
   sumscore counts.
5. **LWL.** Peekbank 2022.1 pull keyed by lab_subject_id; RT = latency of the first shift
   from distractor to target (D–T trials, 40 Hz, +1 sample); one row per administration =
   mean (log) RT; SEEDLingS from EyeLink reports (Zhu et al.): Fernald RT window 300–1800
   ms, accuracy window 367–3500 ms, low-look trials excluded, per-session mean over valid
   trials. Citation for the SEEDLingS LWL data is Zhu et al. (2026), not ref. 49.
6. **Input.** LENA adult word count per hour at fixed ages (AM2018 16 + 18 mo; FMW2013
   18 + 24 mo; SEEDLingS monthly 6–17 mo), 808 recordings; BabyView per-video adult
   (FEM + MAL) alphabetic tokens at transcript confidence ≥ 0.3, videos > 60 s, per hour
   of video, 5,739 videos. Each recording/video is one observation of the measurement
   model; measurement SD per instrument (LENA vs head-cam); study means per dataset.
7. **Demographic variables.** Wordbank `sex`; caregiver education → years (None 0,
   Primary 6, Some Secondary 9, Secondary 12, Some College 14, College 16, Some Graduate
   18, Graduate 20), z-scored within language.
8. **Sample composition.** No preterm / late-talker filter exists for the io/proc
   datasets (see error 18); state the actual composition.

### Models and estimation
9. **Stan model.** η = ξ_i + log H + κ_i log(t/a₀) − δ_j with H fixed (log H = 5.9),
   κ = 1 + δ, non-centred child effects, lexical-class hierarchy on δ_j, LWL and input
   measurement blocks as in Fig. 1A, Poisson-binomial (normal-approx.) count likelihood.
   Priors: δ ~ N(·, 10); coupling coefficients scaled to reference SDs (prior SD = 1 /
   reference SD: σ_r,ref = 0.44, σ_rt0,ref = 0.143, i.e. τ = 1 per reference SD);
   half-normal σ priors with bundle scales; LWL priors (τ_s, ψ, σ). Ladder rungs
   D′0 (input coupling only) → D′1 (+ β_ξ) → D′2 (+ β_κ); reported: joint D′2, io-only
   D′0 of the io model, proc-only D′2 of the proc model.
10. **Sampling.** CmdStan 2.38.0 / cmdstanr 0.9.0, 4 chains × 1,000 warmup + 1,000
    sampling, adapt_delta 0.95, max treedepth 10, reduce_sum threading, seed 1, informed
    inits. Imputed fits: 2,000 + 1,000, adapt_delta 0.97, warm-started from 1k runs.
11. **Convergence.** Joint: max r̂ 1.03, min bulk ESS 163 (structural parameters);
    io-only 1.03 / 88; proc-only 1.04 / 112. Imputed: σ_ξ / π_α r̂ 1.02–1.16, ESS 19–152
    (Norwegian mixes slowly on the (α, ζ) ridge); whole-summary max r̂ up to 2.0 is the
    prior-pinned onset nuisance σ_s only. Report honestly in the SI.
12. **Per-SD effect definitions.** Input → efficiency per SD = σ_r (coefficient fixed at
    1); input → acceleration = γ_in σ_r; processing effects = β × σ_rt0, sign-flipped so
    faster processing is positive; acceleration in κ units, × log(30/21) = 0.36 for the
    level-equivalent effect at 30 months.
13. **Three-arm design.** Same bundle, three Stan files (`log_irt_long_io_count`,
    `log_irt_long_proc_count`, `log_irt_long_proc_bilingual`); the "largely unchanged"
    claim should carry the table in experiments.md #43 (95% addendum).
14. **Imputed model.** `log_irt_long.stan`, variant `no_freq_slopes` (per-child κ, no
    word-frequency term, onset pinned); ξ_i = μ_r + log r_i + log α_i, σ_r fixed by
    override at 0.35 / 0.44 / 0.58; π_α = σ_α² / (σ_α² + σ_r²), share = σ_r² / σ_ξ²; σ_ξ
    for the analytic curve from the 0.44 fits (EN 1.85, NO 2.49; shares EN 3.6 / 5.7 /
    9.9%, NO 2.0 / 3.1 / 5.4%); anchor provenance (0.35 direct estimate, 0.44 population
    re-anchor, 0.58 literature; band lower edge 0.32 = direct estimate's lower 95%).
15. **Longitudinal demographics.** lme4 1.1.37 `glmer`, D_log: `produces ~ log_age +
    (1 + log_age | child) + (1 | item)`, binomial, nAGQ = 0, a₀ = median age per dataset;
    BLUPs = conditional modes; per-dataset `lm(BLUP ~ predictor)` (sex; z-scored
    education years); English datasets pooled by inverse variance in Fig. 2A/B; metafor
    4.8.0 REML random-effects meta with k = 4 (sex) / 2 (maternal ed).
16. **Cross-sectional demographics.** Per language `produces ~ predictor × log(age/a₀) +
    (1 | item) + (1 | child)`, nAGQ = 0, one randomly chosen administration per child
    (seed 2026), a₀ = median age; language list = monolingual-TD with ≥ 200 children with
    sex (31 fitted), ≥ 300 displayed (28 sex / 12 maternal ed); English (British) and
    French (Québec) excluded; French (French) dropped from maternal ed; sex reported as
    female − male; maternal ed needs ≥ 1,000 item rows; REML meta across languages; the
    README's validation note that cross-sectional acceleration runs 1.5–2.3× the
    longitudinal (upper bound).
17. **Figure conventions.** Fig. 2: acceleration CIs truncated at ±(max |estimate| + 0.7)
    in the forests; bands = cross-sectional meta 95% CI; panel C on its own scale.
    Fig. 1B: reference child 0.5 SD below mean efficiency; +1 SD of ξ (2.02) and κ (4.40);
    words produced = 680 × logistic((θ − 15) / 2.5), an illustrative item bank. Fig. 3:
    loess (span 1); BabyView videos averaged to monthly means; RT panel = session means on
    a log axis; input in study-specific units. Fig. 4 / Fig. 5: 95% posterior intervals.
18. **Software.** R 4.5.1, cmdstanr 0.9.0 / CmdStan 2.38.0, lme4 1.1.37, metafor 4.8.0,
    wordbankr 2.0.0 (Redivis), peekbank 2022.1; repo langcog/standard-model-2.
