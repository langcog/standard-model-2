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
