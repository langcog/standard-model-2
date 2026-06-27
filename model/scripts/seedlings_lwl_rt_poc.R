## PROOF OF CONCEPT: derive Fernald-style RT from the SEEDLings hand-tailored
## (Study 1, our n=44) EyeLink reports. Target onset = EL_BUTTON_CRIT_WORD
## message; RT = first TARGET-AOI fixation after onset on DISTRACTOR-initial
## trials, 300-1800 ms window. Reports yield by age to test the 13-18mo overlap
## with our Fernald RT window. RUN LOCALLY.
## Out: figs/data_checks/seedlings_lwl_rt_poc.png + console summary
suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })
DIR <- here("data/raw/seedlings/raw_eyetracking_data/HaT")
RT_LO <- 300; RT_HI <- 1800

trim <- function(x) sub("\\s+$", "", x)
fix <- fread(file.path(DIR, "fix_rep_eighttoeighteenmonth_7-26-2016.xls"), sep = "\t",
             select = c("RECORDING_SESSION_LABEL","TRIAL_INDEX","CURRENT_FIX_START",
                        "CURRENT_FIX_END","CURRENT_FIX_INTEREST_AREA_LABEL"))
setnames(fix, c("sess","trial","fstart","fend","aoi"))
fix[, aoi := trim(aoi)][, `:=`(fstart = as.numeric(fstart), fend = as.numeric(fend))]

mes <- fread(file.path(DIR, "mes_rep_eighttoeighteenmonth_7-26-2016.xls"), sep = "\t",
             select = c("RECORDING_SESSION_LABEL","TRIAL_INDEX","CURRENT_MSG_TEXT",
                        "CURRENT_MSG_TIME","TargetSide","TrialType","AudioTarget"))
setnames(mes, c("sess","trial","msg","mtime","side","ttype","audio"))

## per-trial onset (critical word) + trial attributes (constant within trial)
onset <- mes[msg == "EL_BUTTON_CRIT_WORD", .(onset = min(as.numeric(mtime))), by = .(sess, trial)]
attr_ <- mes[, .(side = side[1], ttype = ttype[1], audio = audio[1]), by = .(sess, trial)]
trials <- merge(onset, attr_, by = c("sess","trial"))
trials <- trials[side %in% c("L","R")]            # test trials only (drop mono warm-ups)
trials[, age := as.integer(sub(".*_0?", "", sub("[a-z]$","",sess)))]

## merge onset onto the fixations, then compute RT per (sess,trial) group
fix2 <- merge(fix, trials[, .(sess, trial, onset, side, ttype, age)],
              by = c("sess","trial"))
res <- fix2[order(sess, trial, fstart), {
  o <- onset[1]
  ini <- aoi[fstart <= o & fend > o]                       # fixation ongoing at onset
  if (!length(ini)) { b <- which(fstart <= o); ini <- if (length(b)) aoi[tail(b,1)] else NA_character_ }
  tgt <- fstart[aoi == "TARGET" & fstart > o]
  .(age = age[1], side = side[1], ttype = ttype[1], ini_aoi = ini[1],
    rt = if (length(tgt)) min(tgt) - o else NA_real_)
}, by = .(sess, trial)]
res[, dinit := ini_aoi == "DISTRACTOR"]
res[, valid := !is.na(dinit) & dinit & !is.na(rt) & rt >= RT_LO & rt <= RT_HI]

## ---- summary ----
cat(sprintf("\nHaT test trials (L/R): %d across %d sessions, %d kids\n",
            nrow(res), uniqueN(res$sess), uniqueN(sub("_.*","",res$sess))))
cat(sprintf("  initial-look: TARGET %d (%.0f%%), DISTRACTOR %d (%.0f%%), away/. %d\n",
            sum(res$ini_aoi=="TARGET",na.rm=T), 100*mean(res$ini_aoi=="TARGET",na.rm=T),
            sum(res$dinit,na.rm=T), 100*mean(res$dinit,na.rm=T),
            sum(!res$ini_aoi %in% c("TARGET","DISTRACTOR"))))
cat(sprintf("  valid D-initial RTs in [%d,%d]: %d (%.0f%% of D-initial)\n",
            RT_LO, RT_HI, sum(res$valid), 100*sum(res$valid)/sum(res$dinit,na.rm=T)))
v <- res[valid==TRUE]
cat(sprintf("  RT ms: median=%.0f IQR=[%.0f,%.0f] mean=%.0f\n",
            median(v$rt), quantile(v$rt,.25), quantile(v$rt,.75), mean(v$rt)))
cat("\nvalid-RT yield by age (mo):\n")
print(v[, .(n_RT = .N, n_kids = uniqueN(sub("_.*","",sess)), med_RT = round(median(rt))), by = age][order(age)])
cat("\n  in-window ages (14/16/18) = overlap with our 13-30 Fernald window:\n")
print(v[age>=14, .(n_RT=.N, n_kids=uniqueN(sub("_.*","",sess))), by=.(age,ttype)][order(age,ttype)])
cat("\n  RT by trial type (generic vs hand-picked high-freq 'unique') -- the freq->processing signal:\n")
print(v[, .(n_RT=.N, med_RT=round(median(rt)), mean_RT=round(mean(rt))), by=ttype][order(ttype)])

## ---- plot: RT distribution by age ----
v[, age_f := factor(age)]
p <- ggplot(v, aes(age_f, rt)) +
  geom_violin(fill="#1f78b4", alpha=.3, color=NA) +
  geom_boxplot(width=.25, outlier.size=.4) +
  geom_hline(yintercept=median(v$rt), linetype=2, color="grey40") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=RT_LO, ymax=RT_HI, alpha=0) +
  labs(x="age (months)", y="derived RT (ms, D-initial, 300-1800)",
       title="SEEDLings LWL: derived Fernald RT by age (POC)",
       subtitle=sprintf("%d valid RTs / %d kids; ages 14-18mo overlap our 13-30mo Fernald window",
                        nrow(v), uniqueN(sub("_.*","",v$sess)))) +
  theme_bw(base_size=11)
ggsave(here("figs/data_checks/seedlings_lwl_rt_poc.png"), p, width=8, height=5, dpi=130)
cat("\nwrote figs/data_checks/seedlings_lwl_rt_poc.png\n")
