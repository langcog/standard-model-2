## Subset + reindex a built joint io-proc bundle down to a tiny verification set
## (default 50 kids / 50 words) so Phase-0 simplification fits run in minutes, not hours.
## Children are prioritized for CHANNEL COVERAGE (input + LWL + multi-timepoint CDI) so
## gamma_in / beta_xi stay checkable; items are stratified across lexical classes.
## Usage: Rscript model/scripts/build_subsample_bundle.R [n_kids] [n_items] [in_bundle]
## Output: fits/<stem>_sub<K>x<J>.rds
suppressPackageStartupMessages({ library(dplyr) })
args   <- commandArgs(trailingOnly = TRUE)
NK     <- as.integer(if (length(args) >= 1) args[1] else 50)
NJ     <- as.integer(if (length(args) >= 2) args[2] else 50)
IN     <- if (length(args) >= 3) args[3] else "fits/joint_io_proc_mm_subset_data.rds"
set.seed(1)

b  <- readRDS(IN); sd <- b$stan_data; I <- sd$I; J <- sd$J
## ---- pick children: channel coverage first ----
has_in  <- (1:I) %in% unique(sd$rec_to_child)
has_lwl <- (1:I) %in% unique(sd$lwl_to_child)
n_adm   <- tabulate(sd$admin_to_child, nbins = I)
cdi_ok  <- n_adm > 0
score   <- has_in + has_lwl + (n_adm >= 2)              # 0..3
keep_i  <- (1:I)[cdi_ok][order(-score[cdi_ok], -n_adm[cdi_ok])][1:NK]
keep_i  <- sort(keep_i)
## ---- pick items: stratified across lexical classes ----
cc <- sd$cc
keep_j <- unlist(lapply(sort(unique(cc)), function(c) {
  idx <- which(cc == c); sample(idx, min(length(idx), ceiling(NJ / length(unique(cc)))))
}))
keep_j <- sort(sample(keep_j, min(NJ, length(keep_j))))

## ---- reindex maps ----
imap <- rep(NA_integer_, I); imap[keep_i] <- seq_along(keep_i)
jmap <- rep(NA_integer_, J); jmap[keep_j] <- seq_along(keep_j)
keep_a <- which(!is.na(imap[sd$admin_to_child]))        # admins of kept kids
amap   <- rep(NA_integer_, sd$A); amap[keep_a] <- seq_along(keep_a)
keep_n <- which(!is.na(amap[sd$aa]) & !is.na(jmap[sd$jj]))   # CDI obs in kept admins & items
keep_v <- which(!is.na(imap[sd$rec_to_child]))          # input recs of kept kids
keep_l <- which(!is.na(imap[sd$lwl_to_child]))          # LWL obs of kept kids

n <- sd                                                  # copy, then overwrite
n$I <- NK; n$J <- length(keep_j); n$A <- length(keep_a)
n$N <- length(keep_n); n$V_obs <- length(keep_v); n$N_lwl <- length(keep_l)
n$aa <- amap[sd$aa[keep_n]];  n$jj <- jmap[sd$jj[keep_n]];  n$y <- sd$y[keep_n]
n$admin_to_child <- imap[sd$admin_to_child[keep_a]];  n$admin_age <- sd$admin_age[keep_a]
n$cc <- sd$cc[keep_j];  n$log_p <- sd$log_p[keep_j]
n$study_of_child <- sd$study_of_child[keep_i]
n$rec_to_child <- imap[sd$rec_to_child[keep_v]];  n$log_input_obs <- sd$log_input_obs[keep_v]
n$instr_of_rec <- sd$instr_of_rec[keep_v]
n$lwl_to_child <- imap[sd$lwl_to_child[keep_l]];  n$lwl_log_age <- sd$lwl_log_age[keep_l]
n$lwl_log_rt   <- sd$lwl_log_rt[keep_l]
THREADS <- 4L; n$grainsize <- max(1L, as.integer(floor(n$N / (2 * THREADS))))

## ---- validate index ranges ----
stopifnot(all(n$aa %in% 1:n$A), all(n$jj %in% 1:n$J), all(n$admin_to_child %in% 1:n$I),
          all(n$rec_to_child %in% 1:n$I), all(n$lwl_to_child %in% 1:n$I),
          all(n$cc %in% 1:n$C), all(n$study_of_child %in% 1:n$S))

b$stan_data <- n
stem <- sub("\\.rds$", "", basename(IN)); stem <- sub("_subset_data$|_mm_subset_data$", "", stem)
out <- sprintf("fits/%s_sub%dx%d.rds", stem, NK, n$J)
saveRDS(b, out)
cat(sprintf("wrote %s\n  I=%d J=%d A=%d N=%d N_lwl=%d V_obs=%d  | studies=%d classes=%d\n",
            out, n$I, n$J, n$A, n$N, n$N_lwl, n$V_obs, length(unique(n$study_of_child)), length(unique(n$cc))))
cat(sprintf("  channel coverage of the %d kids: input=%d  LWL=%d  both=%d  multi-tp-CDI=%d\n",
            NK, sum(has_in[keep_i]), sum(has_lwl[keep_i]),
            sum(has_in[keep_i] & has_lwl[keep_i]), sum(n_adm[keep_i] >= 2)))
