## Longitudinal test of per-child growth-rate variance using Wordbank
## admin-level production counts as a direct measure of ability.
##
## Strategy:
##   - Read admins.feather. Filter to longitudinal children in WS forms for
##     English (American) and Norwegian.
##   - Use logit(production_proportion) as the outcome (ability on logit
##     scale). Standardize proportion to avoid extremes.
##   - Fit two lme4 models:
##       m_int: logit_prop ~ log_age + (1 | child_id)                # intercepts only
##       m_slp: logit_prop ~ log_age + (log_age | child_id)          # + random slope
##     LRT tests whether children have significant variance in growth rate.
##   - Report SD(intercept), SD(slope), correlation, and p-value for LRT.
##
## This is the cleanest available test of the sigma_zeta hypothesis using
## real longitudinal data (no reliance on cross-sectional heteroskedasticity).

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(lme4)
})

ADMINS_PATH <- "/Users/mcfrank/Projects/wordbank/wordbank-book/data/_common/admins.feather"
# WS has 680 items in English; Norwegian WS is similar. We'll use form-
# specific item counts via `items.feather`.
ITEMS_PATH  <- "/Users/mcfrank/Projects/wordbank/wordbank-book/data/_common/items.feather"

admins <- read_feather(ADMINS_PATH)
items  <- read_feather(ITEMS_PATH)

# Item count per (language, form, item_kind=word)
word_counts <- items %>%
  filter(item_kind == "word") %>%
  count(language, form, name = "n_items")

adm <- admins %>%
  inner_join(word_counts, by = c("language", "form"))

# Identify longitudinal children (WS only for now) and filter
long_ids <- adm %>%
  filter(form == "WS",
         language %in% c("English (American)", "Norwegian")) %>%
  count(child_id, language, name = "n_admins") %>%
  filter(n_admins >= 2) %>%
  pull(child_id)

d_long <- adm %>%
  filter(form == "WS",
         language %in% c("English (American)", "Norwegian"),
         child_id %in% long_ids,
         !is.na(production), !is.na(age)) %>%
  mutate(prop = (production + 0.5) / (n_items + 1),
         logit_prop = qlogis(pmin(pmax(prop, 1e-4), 1 - 1e-4)),
         # Center log-age at log(24) so the intercept is interpretable
         # as ability at 24 months and isn't an extrapolation.
         log_age    = log(age) - log(24))

cat(sprintf("\nLongitudinal WS data:\n"))
cat(sprintf("  rows: %d\n", nrow(d_long)))
cat(sprintf("  children: %d\n",
            length(unique(paste(d_long$language, d_long$child_id)))))
cat("  by language:\n")
print(d_long %>% distinct(language, child_id) %>% count(language))

cat("\n  age range:", range(d_long$age), "\n")
cat("  production range:", range(d_long$production), "\n")

# Fit by language
results <- list()
for (lang in c("English (American)", "Norwegian")) {
  dl <- d_long %>% filter(language == lang)
  if (nrow(dl) < 30) next
  cat(sprintf("\n===== %s (n_admins=%d, n_children=%d) =====\n",
              lang, nrow(dl), length(unique(dl$child_id))))

  m_int <- lmer(logit_prop ~ log_age + (1 | child_id),
                data = dl, REML = FALSE)
  m_slp <- lmer(logit_prop ~ log_age + (log_age | child_id),
                data = dl, REML = FALSE)

  cat("\nFixed effects (intercept + slope model):\n")
  print(fixef(m_slp))

  cat("\nRandom-effect SDs (m_slp):\n")
  print(VarCorr(m_slp))

  lrt <- anova(m_int, m_slp)
  cat("\nLRT for adding random slope:\n"); print(lrt)

  vc <- as.data.frame(VarCorr(m_slp))
  results[[lang]] <- list(
    lang = lang,
    n_admins = nrow(dl),
    n_children = length(unique(dl$child_id)),
    sd_int = vc$sdcor[vc$grp == "child_id" & vc$var1 == "(Intercept)" &
                      is.na(vc$var2)],
    sd_slp = vc$sdcor[vc$grp == "child_id" & vc$var1 == "log_age" &
                      is.na(vc$var2)],
    cor    = vc$sdcor[vc$grp == "child_id" & !is.na(vc$var2)],
    sd_res = vc$sdcor[vc$grp == "Residual"],
    lrt_p  = lrt$`Pr(>Chisq)`[2],
    fixed  = fixef(m_slp)
  )
}

cat("\n\n===== Summary =====\n")
for (r in results) {
  cat(sprintf("\n%s  (n_admins=%d, n_children=%d)\n",
              r$lang, r$n_admins, r$n_children))
  cat(sprintf("  Population intercept: %.2f, slope (log age): %.2f\n",
              r$fixed[1], r$fixed[2]))
  cat(sprintf("  SD(random intercept) = %.3f\n", r$sd_int))
  cat(sprintf("  SD(random slope)     = %.3f\n", r$sd_slp))
  cat(sprintf("  cor(intercept,slope) = %.3f\n", r$cor))
  cat(sprintf("  residual SD          = %.3f\n", r$sd_res))
  cat(sprintf("  LRT p (add random slope): %.3g\n", r$lrt_p))

  # Interpretation aid: slope SD relative to |fixed slope|
  cat(sprintf("  |SD(slope)/fixed slope| = %.3f (0 = no variation, ~1 = ±100%% around fixed slope)\n",
              r$sd_slp / abs(r$fixed[2])))
}

saveRDS(results,
        "/Users/mcfrank/Projects/standard_model_2/model/fits/longitudinal_slopes.rds")
cat("\nSaved: model/fits/longitudinal_slopes.rds\n")
