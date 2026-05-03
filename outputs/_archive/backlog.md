# Future experiments / backlog

Not urgent — things to consider once the core model story settles.

- **Within-person input variation.** Use BabyView or Seedlings corpora,
  which have dense at-home recordings over time, to estimate the
  within-person σ_r (i.e., how much a single child's input rate fluctuates
  across days/weeks). This is a separate quantity from the cross-person σ_r
  we currently pin from Sperry/HR/WF. Having both would let us decompose
  total input variance into within-child and between-child components.

- **Non-English robustness fits.** Norwegian longitudinal WS data is
  available (1,562 children with ≥2 WS admins in the Kristoffersen
  dataset) — already in our sights via the longitudinal analysis. Could
  also pull cross-sectional subsets for other languages where CDI +
  CHILDES-style frequency data exist, and compare π_α estimates.

- **Comprehension vs. production (WG joint model).** WG form has both
  comp and prod per word per child. Bivariate (ξ_comp, ξ_prod) with
  estimated correlation ρ would tell us whether "ability" is a single
  construct or has modality-specific components.

- **Correlated random effects (ξ, ζ).** The longitudinal LMM shows
  corr(intercept, slope) ≈ +0.72 on log-age scale. The current Stan
  model has them independent. Add an LKJ prior on the correlation if we
  revisit per-child slopes with proper longitudinal data.
