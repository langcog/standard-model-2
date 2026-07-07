## 03_fan.R -- by-dataset model-implied fan (Fig 1 bottom panel) from the M3 fits.
##
## Draws a synthetic population from the fitted M3 posterior -- (xi, kappa) ~
## MVN(mu_xi/1+delta, Sigma[sigma_a, sigma_b, rho]) and item difficulties from
## N(0, tau_delta) -- computes each child's expected-vocab curve, and overlays
## 10/25/50/75/90 quantiles on the empirical trajectories. Uses only the saved
## scalar posteriors (fits/bayes_long/summaries/<slug>_m3.summary.rds) + the
## bundle (empirical data + log_H/a0). Datasets without an m3 fit are skipped.
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(here)})
set.seed(1)
N_SIM <- 500; N_SPAG <- 150; M_ITEM <- 500; QS <- c(.1,.25,.5,.75,.9)
LANGS <- c(thal="English (Thal)", smith="English (Smith)", marchman="English (Marchman)",
           norwegian="Norwegian", japanese="Japanese")

one <- function(slug, label) {
  sf <- here("fits","bayes_long","summaries", paste0(slug,"_m3.summary.rds"))
  bf <- here("fits","bayes_long", paste0("bundle_",slug,".rds"))
  if (!file.exists(sf)) { cat("skip", slug, "(no m3 yet)\n"); return(NULL) }
  s <- as.data.frame(readRDS(sf)); g <- function(v) s$median[s$variable==v]
  b <- readRDS(bf); sd <- b$stan_data
  mu_xi<-g("mu_xi"); delta<-g("delta"); sa<-g("sigma_a"); sb<-g("sigma_b")
  rho<-g("rho_ab"); tau<-g("tau_delta"); log_H<-sd$log_H; a0<-sd$a0

  ## synthetic population: (a,b) ~ MVN(0, Sigma); xi = mu_xi + a; kappa = 1+delta+b
  Sig <- matrix(c(sa^2, rho*sa*sb, rho*sa*sb, sb^2), 2)
  Z <- matrix(rnorm(N_SIM*2), N_SIM, 2) %*% chol(Sig)
  xi <- mu_xi + Z[,1]; kappa <- 1 + delta + Z[,2]
  dj <- rnorm(M_ITEM, 0, tau)                          # synthetic item difficulties
  base_j <- log_H - dj

  ## empirical trajectories (proportion produced per admin)
  emp <- tibble(aa=sd$aa, y=sd$y) |> group_by(aa) |> summarise(prop=mean(y), .groups="drop") |>
    mutate(age=sd$admin_age[aa], child=sd$admin_to_child[aa])
  aw <- quantile(emp$age, c(.05,.95)); ages <- seq(floor(aw[1]), ceiling(aw[2]), by=0.5)
  kids <- sample(unique(emp$child), min(N_SPAG, n_distinct(emp$child)))
  spag <- emp |> filter(child %in% kids) |> mutate(lang=factor(label, levels=LANGS))

  ## model fan: vocab_i(t) = mean_j plogis(xi_i + kappa_i*log(t/a0) + base_j)
  qt <- lapply(ages, function(t){
    A <- log(t/a0); v <- rowMeans(plogis(outer(xi + kappa*A, base_j, "+")))
    data.frame(age=t, q=QS, vocab=quantile(v, QS, names=FALSE))}) |> bind_rows()
  qt$qf <- factor(qt$q, levels=QS, labels=c("10th","25th","50th","75th","90th"))
  qt$lang <- factor(label, levels=LANGS)
  cat(sprintf("%-10s kappa=%.1f sigma_a=%.2f sigma_b=%.2f rho=%.2f\n", slug, 1+delta, sa, sb, rho))
  list(spag=spag, qt=qt)
}

res <- Filter(Negate(is.null), Map(one, names(LANGS), LANGS))
spag <- bind_rows(lapply(res, `[[`, "spag")); qt <- bind_rows(lapply(res, `[[`, "qt"))

pal <- c("10th"="#2c7fb8","25th"="#7fcdbb","50th"="#238b45","75th"="#fe9929","90th"="#c41e37")
p <- ggplot() +
  geom_line(data=spag, aes(age, prop, group=child), color="grey30", alpha=0.12, linewidth=0.22) +
  geom_line(data=qt, aes(age, vocab, color=qf, linewidth=ifelse(qf=="50th",0.95,0.55))) +
  facet_wrap(~lang, nrow=1) +
  scale_color_manual(values=pal, name="Percentile") + scale_linewidth_identity() +
  scale_y_continuous(limits=c(0,1)) +
  labs(x="Age (months)", y="proportion of items produced",
       title="Bayesian M3 model-implied fan over empirical trajectories (by dataset)") +
  theme_minimal(base_size=10) +
  theme(strip.text=element_text(face="bold"), legend.position="bottom",
        panel.grid.minor=element_blank(), plot.title=element_text(size=11, face="bold"))
out <- here("studies","bayes_long","fan_m3.png")
ggsave(out, p, width=3.1*length(res)+0.5, height=3.6, dpi=150)
cat("wrote", out, "\n")
