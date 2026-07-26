## Prototype SI figure: non-parametric GAMLSS quantile fan overlaid on the Bayesian
## M3 *parametric* fan, across the five datasets. Point: the parametric accumulator is
## a faithful compression of the non-parametric quantile structure of the sumscores.
##
## GAMLSS = beta regression (BE) on per-administration proportion-produced, the
## CDI-norms model (Wordbank): a monotone P-spline on the mean, a P-spline on the SD
##   prop ~ pbm(age),  sigma ~ pb(age),  family = BE
## (no survey weighting -- just the gamlss, per MCF). Fit on the SAME per-admin
## sumscores the Bayesian model sees (recovered from the bayes_long bundles), so the
## two fans are a true apples-to-apples overlay on the existing Fig-1 fan.
##
## Usage:  Rscript studies/gamlss/01_gamlss_overlay.R
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(here); library(gamlss)
})
set.seed(1)
SFX   <- Sys.getenv("SLUG_SUFFIX","_a3")
N_SIM <- 500; N_SPAG <- 120; M_ITEM <- 500
QS <- c(.1,.25,.5,.75,.9); PCT <- c(10,25,50,75,90)
LANGS <- c(thal="English (Thal)", smith="English (Smith)", marchman="English (Marchman)",
           norwegian="Norwegian", japanese="Japanese")

one <- function(slug, label) {
  sf <- here("fits","bayes_long","summaries", paste0(slug,SFX,"_m3.summary.rds"))
  bf <- here("fits","bayes_long", paste0("bundle_",slug,SFX,".rds"))
  if (!file.exists(sf) || !file.exists(bf)) { cat("skip", slug, "(missing fit/bundle)\n"); return(NULL) }
  s <- as.data.frame(readRDS(sf)); g <- function(v) s$median[s$variable==v]
  b <- readRDS(bf); sd <- b$stan_data
  mu_xi<-g("mu_xi"); delta<-g("delta"); sa<-g("sigma_a"); sb<-g("sigma_b")
  rho<-g("rho_ab"); tau<-g("tau_delta"); log_H<-sd$log_H; a0<-sd$a0; J<-sd$J

  ## shared data: per-administration proportion produced (the sumscore / J)
  emp <- tibble(aa=sd$aa, y=sd$y) |> group_by(aa) |> summarise(prop=mean(y), .groups="drop") |>
    mutate(age=sd$admin_age[aa], child=sd$admin_to_child[aa])
  aw <- quantile(emp$age, c(.02,.98)); ages <- seq(floor(aw[1]), ceiling(aw[2]), by=0.5)

  ## ---- Bayesian M3 parametric fan (synthetic population, as in 03_fan.R) ----
  Sig <- matrix(c(sa^2, rho*sa*sb, rho*sa*sb, sb^2), 2)
  Z <- matrix(rnorm(N_SIM*2), N_SIM, 2) %*% chol(Sig)
  xi <- mu_xi + Z[,1]; kappa <- 1 + delta + Z[,2]
  base_j <- log_H - rnorm(M_ITEM, 0, tau)
  bayes <- lapply(ages, function(t){ A<-log(t/a0)
    v <- rowMeans(plogis(outer(xi + kappa*A, base_j, "+")))
    data.frame(age=t, pct=PCT, vocab=quantile(v, QS, names=FALSE)) }) |> bind_rows() |>
    mutate(model="Bayesian M3 (parametric)")

  ## ---- GAMLSS non-parametric fan (BE beta regression on the same props) ----
  ## squeeze exact 0/1 to the nearest representable interior value (BE needs open (0,1))
  ## gamlss's predict/centiles.pred re-evaluate the training data BY NAME in the
  ## global env, so it must live there (not as a function-local variable).
  assign("GDF", emp |> mutate(yp = pmin(pmax(prop, 0.5/J), 1 - 0.5/J)) |> as.data.frame(),
         envir = .GlobalEnv)
  gm <- tryCatch(gamlss(yp ~ pbm(age, lambda=10000), sigma.formula = ~ pb(age),
                        family = BE, data = GDF, trace = FALSE),
                 error=function(e){ cat("gamlss failed for", slug, ":", conditionMessage(e), "\n"); NULL })
  gaml <- NULL
  if (!is.null(gm)) {
    cp <- centiles.pred(gm, cent=PCT, xname="age", xvalues=ages)   # cols: x, "10","25",...
    gaml <- cp |> as_tibble() |> pivot_longer(-x, names_to="pct", values_to="vocab") |>
      transmute(age=x, pct=as.integer(pct), vocab, model="GAMLSS (non-parametric)")
  }

  fan <- bind_rows(bayes, gaml) |>
    mutate(pct=factor(pct, levels=PCT), lang=factor(label, levels=LANGS))
  kids <- sample(unique(emp$child), min(N_SPAG, n_distinct(emp$child)))
  spag <- emp |> filter(child %in% kids) |> mutate(lang=factor(label, levels=LANGS))
  cat(sprintf("%-10s kappa=%.1f  gamlss=%s\n", slug, 1+delta, ifelse(is.null(gm),"FAILED","ok")))
  list(fan=fan, spag=spag)
}

res  <- Filter(Negate(is.null), Map(one, names(LANGS), LANGS))
fan  <- bind_rows(lapply(res, `[[`, "fan"))
spag <- bind_rows(lapply(res, `[[`, "spag"))

pal <- c("10"="#2c7fb8","25"="#7fcdbb","50"="#238b45","75"="#fe9929","90"="#c41e37")
p <- ggplot() +
  geom_line(data=spag, aes(age, prop, group=child), color="grey35", alpha=0.10, linewidth=0.2) +
  geom_line(data=fan, aes(age, vocab, color=pct, linetype=model, linewidth=model)) +
  facet_wrap(~lang, nrow=1) +
  scale_color_manual(values=pal, name="Percentile") +
  scale_linetype_manual(values=c("Bayesian M3 (parametric)"="solid","GAMLSS (non-parametric)"="21"), name=NULL) +
  scale_linewidth_manual(values=c("Bayesian M3 (parametric)"=0.9,"GAMLSS (non-parametric)"=0.6), guide="none") +
  scale_y_continuous(limits=c(0,1)) +
  labs(x="Age (months)", y="proportion of items produced",
       title="Parametric (Bayesian M3) vs non-parametric (GAMLSS beta-regression) quantile fans",
       subtitle="solid = Bayesian accumulator fan; dashed = GAMLSS centiles on the same sumscores") +
  theme_minimal(base_size=10) +
  theme(strip.text=element_text(face="bold"), legend.position="bottom",
        panel.grid.minor=element_blank(), plot.title=element_text(size=11, face="bold"))
out <- here("studies","gamlss",paste0("fig_gamlss_overlay", if(SFX=="") "_2plus" else SFX, ".png"))
ggsave(out, p, width=3.1*length(res)+0.6, height=3.8, dpi=150)
cat("wrote", out, "\n")
