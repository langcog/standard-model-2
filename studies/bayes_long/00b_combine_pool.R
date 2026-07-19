## 00b_combine_pool.R -- combine the per-dataset bundles into ONE cross-dataset bundle
## for the mega-model (m_pool.stan). Global child/item/admin indices + dataset-membership
## arrays (child_ds, item_ds). Items and children never cross datasets.
##
## Usage:  Rscript studies/bayes_long/00b_combine_pool.R [suffix]   (default suffix "_a3")
suppressPackageStartupMessages({library(dplyr); library(here)})
SFX <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "_a3"
SLUGS <- c("thal","smith","marchman","norwegian","japanese")
BL <- here("fits","bayes_long")

aa<-jj<-y<-integer(0); a2c<-integer(0); age<-numeric(0)
child_ds<-item_ds<-integer(0)
oI<-oJ<-oA<-0L                                  # running child / item / admin offsets
meta_rows <- list()
for (d in seq_along(SLUGS)) {
  b <- readRDS(file.path(BL, sprintf("bundle_%s%s.rds", SLUGS[d], SFX))); sd <- b$stan_data
  aa  <- c(aa,  sd$aa + oA)
  jj  <- c(jj,  sd$jj + oJ)
  y   <- c(y,   sd$y)
  a2c <- c(a2c, sd$admin_to_child + oI)         # admin -> GLOBAL child
  age <- c(age, sd$admin_age)
  child_ds <- c(child_ds, rep(d, sd$I))
  item_ds  <- c(item_ds,  rep(d, sd$J))
  meta_rows[[d]] <- data.frame(dataset=SLUGS[d], I=sd$I, J=sd$J, A=sd$A, N=sd$N)
  oI <- oI + sd$I; oJ <- oJ + sd$J; oA <- oA + sd$A
}
b1 <- readRDS(file.path(BL, sprintf("bundle_%s%s.rds", SLUGS[1], SFX)))$stan_data

stan_data <- list(
  N=length(y), A=oA, I=oI, J=oJ, D=length(SLUGS),
  aa=aa, jj=jj, y=y, admin_to_child=a2c, admin_age=age,
  child_ds=child_ds, item_ds=item_ds,
  log_H=b1$log_H, a0=b1$a0)

meta <- list(slug=paste0("pool",SFX), datasets=SLUGS, per=bind_rows(meta_rows),
             n_kids=oI, n_items=oJ, n_admins=oA, n_obs=length(y))
cat(sprintf("pool%s: D=%d kids=%d items=%d admins=%d obs=%d\n", SFX, length(SLUGS), oI, oJ, oA, length(y)))
print(meta$per)
saveRDS(list(stan_data=stan_data, meta=meta), file.path(BL, sprintf("bundle_pool%s.rds", SFX)))
cat("wrote", file.path(BL, sprintf("bundle_pool%s.rds", SFX)), "\n")
