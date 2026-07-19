## build_cache_si_gamlss.R -- SI figures for the GAMLSS/fitted-kappa comparison.
## Reads the 2+ M3 fits (bundles + *_m3_{child,psi}.csv + summaries) and caches everything
## the SI needs, so the qmd renders without the fits:
##   paper/cache/si_gamlss.rds = list(overlay, kappa, contrast)
## Run:  Rscript paper/build_cache_si_gamlss.R
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(here); library(gamlss)})
set.seed(1)
SUMM <- here("fits","bayes_long","summaries"); BL <- here("fits","bayes_long")
LANGS <- c(thal="English (Thal)", smith="English (Smith)", marchman="English (Marchman)",
           norwegian="Norwegian", japanese="Japanese")
QS <- c(.1,.5,.9); PCT <- c(10,50,90)

one <- function(slug, label) {
  b <- readRDS(file.path(BL, sprintf("bundle_%s.rds", slug))); sd <- b$stan_data; a0 <- sd$a0; lH <- sd$log_H
  ch <- read.csv(file.path(SUMM, sprintf("%s_m3_child.csv", slug)))
  dj <- read.csv(file.path(SUMM, sprintf("%s_m3_psi.csv", slug)))$delta_j
  s <- as.data.frame(readRDS(file.path(SUMM, sprintf("%s_m3.summary.rds", slug)))); g <- function(v) s$median[s$variable==v]
  mu<-g("mu_xi"); dl<-g("delta"); sa<-g("sigma_a"); sb<-g("sigma_b"); rho<-g("rho_ab")

  emp <- tibble(aa=sd$aa, y=sd$y) |> group_by(aa) |> summarise(prop=mean(y), .groups="drop") |> mutate(age=sd$admin_age[aa])
  aw <- quantile(emp$age, c(.03,.97)); ages <- seq(floor(aw[1]), ceiling(aw[2]), by=1.5)
  empq <- function(a,q){ v<-emp$prop[abs(emp$age-a)<=1.6]; if(length(v)<10) NA else quantile(v,q) }
  ## fitted-kappa fan: quantiles over the ACTUAL fitted children, actual item difficulties
  fanfit <- function(xi,kap) t(sapply(ages, function(t){ A<-log(t/a0)
    quantile(vapply(seq_along(xi), function(i) mean(plogis(xi[i]+kap[i]*A - dj + lH)), 0.0), QS) }))
  fit <- fanfit(ch$xi_median, ch$kappa_median)
  ## MVN(sigma_b) parametric fan (the over-disperser)
  Zc <- matrix(rnorm(4000),2000,2) %*% chol(matrix(c(sa^2,rho*sa*sb,rho*sa*sb,sb^2),2))
  mvn <- fanfit(mu+Zc[,1], 1+dl+Zc[,2])
  ## GAMLSS auto-lambda
  assign("GDF", emp |> mutate(yp=pmin(pmax(prop,0.5/sd$J),1-0.5/sd$J)) |> as.data.frame(), envir=.GlobalEnv)
  gm <- gamlss(yp ~ pbm(age), sigma.formula=~pb(age), family=BE, data=GDF, trace=FALSE)
  gc <- centiles.pred(gm, cent=PCT, xname="age", xvalues=ages)[,2:4]

  mk <- function(m,src) data.frame(age=ages, q10=m[,1], q50=m[,2], q90=m[,3], model=src, dataset=label)
  lines <- bind_rows(mk(fit,"M3 (fitted-κ)"), mk(gc,"GAMLSS (auto-λ)"), mk(mvn,"M3 (MVN draw)"))
  empdf <- expand.grid(age=ages, pct=PCT) |> rowwise() |> mutate(v=empq(age, pct/100)) |> ungroup() |> filter(!is.na(v)) |> mutate(dataset=label)
  cat(sprintf("%-10s kappa=%.1f sigma_b=%.2f  fitted-kappa SD=%.2f  skew=%.2f\n", slug, 1+dl, sb, sd(ch$kappa_median),
              {x<-ch$kappa_median-mean(ch$kappa_median); mean(x^3)/mean(x^2)^1.5}))
  list(lines=lines, emp=empdf, kappa=tibble(dataset=label, kappa=ch$kappa_median))
}

R <- Map(one, names(LANGS), LANGS)
overlay <- list(lines = bind_rows(lapply(R,`[[`,"lines")), emp = bind_rows(lapply(R,`[[`,"emp")))
kappa   <- bind_rows(lapply(R,`[[`,"kappa"))
saveRDS(list(overlay=overlay, kappa=kappa, levels=unname(LANGS)), here("paper","cache","si_gamlss.rds"))
cat("wrote paper/cache/si_gamlss.rds\n")
