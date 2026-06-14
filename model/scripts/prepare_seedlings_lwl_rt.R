## Production extractor: derive per-child Fernald RT from the SEEDLings
## hand-tailored LWL (Zhu et al., Study 1, our n=44) EyeLink reports, with a QC
## report card against the paper's trial-funnel + accuracy targets
## (journal/notes/seedlings_lwl_qc_targets.md). RT enters the io-proc model at the
## CHILD level (per-kid intercept+slope), so we emit per-trial RT rows in the
## `lwl` schema and let sigma_lwl absorb trial noise. RUN LOCALLY.
## Out: data/seedlings/seedlings_lwl_rt.csv  +  console QC report card
suppressPackageStartupMessages({ library(here); library(data.table) })
DIR <- here("data/seedlings/raw_eyetracking_data/HaT")
W_LO <- 367; W_HI <- 3500; W_LEN <- W_HI - W_LO            # accuracy window
RT_LO <- 300; RT_HI <- 1800                                # Fernald RT window
trim <- function(x) sub("\\s+$", "", x)
qc <- function(label, got, target, tol) cat(sprintf("  %-34s got %-7s target %-7s %s\n",
      label, got, target, if (abs(got-target) <= tol) "PASS" else sprintf("off by %+d", got-target)))

## ---- load reports ----
fix <- fread(file.path(DIR,"fix_rep_eighttoeighteenmonth_7-26-2016.xls"), sep="\t",
  select=c("RECORDING_SESSION_LABEL","TRIAL_INDEX","CURRENT_FIX_START","CURRENT_FIX_END",
           "CURRENT_FIX_INTEREST_AREA_LABEL"))
setnames(fix, c("sess","trial","fstart","fend","aoi"))
fix[, `:=`(aoi=trim(aoi), fstart=as.numeric(fstart), fend=as.numeric(fend))]
mes <- fread(file.path(DIR,"mes_rep_eighttoeighteenmonth_7-26-2016.xls"), sep="\t",
  select=c("RECORDING_SESSION_LABEL","TRIAL_INDEX","CURRENT_MSG_TEXT","CURRENT_MSG_TIME",
           "TargetSide","TrialType","Pair","TargetImage","DistractorImage"))
setnames(mes, c("sess","trial","msg","mtime","side","ttype","pair","timg","dimg"))

## ---- trial table: onset + attributes ----
onset <- mes[msg=="EL_BUTTON_CRIT_WORD", .(onset=min(as.numeric(mtime))), by=.(sess,trial)]
att <- mes[, .(side=side[1], ttype=ttype[1], pair=pair[1], timg=timg[1], dimg=dimg[1]), by=.(sess,trial)]
att[, age := as.integer(sub(".*_0?","",sub("[a-z]$","",sess)))]
## zero-pad subj to "01".."46" -- fixes the lone mis-padded session "9_10" -> "09"
att[, subj := sprintf("%02d", as.integer(sub("_.*","",sess)))]
trials <- merge(att, onset, by=c("sess","trial"), all.x=TRUE)        # all attempts; onset NA if no button
test <- trials[side %in% c("L","R")]                                 # drop mono warm-ups

## ---- Tier 1: trial funnel ----
n_fix_per_trial <- fix[, .N, by=.(sess,trial)]
test <- merge(test, n_fix_per_trial, by=c("sess","trial"), all.x=TRUE)
test[is.na(N), N := 0L]
no_onset <- test[is.na(onset), .N]
no_fix   <- test[!is.na(onset) & N==0, .N]
## within-window looking time at target/distractor (clip each fixation to the window)
fw <- merge(fix, test[!is.na(onset), .(sess,trial,onset)], by=c("sess","trial"))
fw[, ov := pmax(0, pmin(fend, onset+W_HI) - pmax(fstart, onset+W_LO))]
look <- fw[aoi %in% c("TARGET","DISTRACTOR"), .(t_td = sum(ov),
            t_tgt=sum(ov[aoi=="TARGET"]), t_dis=sum(ov[aoi=="DISTRACTOR"])), by=.(sess,trial)]
test <- merge(test, look, by=c("sess","trial"), all.x=TRUE)
test[is.na(t_td), `:=`(t_td=0, t_tgt=0, t_dis=0)]
lowlook <- test[!is.na(onset) & N>0 & t_td < W_LEN/3, .N]
usable  <- test[!is.na(onset) & N>0 & t_td >= W_LEN/3]
cat("\n=== Tier 1: trial funnel (target = Zhu et al. Study 1) ===\n")
qc("planned trials (44x6x32)", 44L*6L*32L, 8448L, 0L)
qc("trial attempts (test, L/R)", nrow(test), 8383L, 120L)       # ~tol: warmups/structure
qc("no locatable onset (no button)", no_onset, 78L, 25L)
qc("no fixation data", no_fix, 190L, 60L)
qc("low-look (<1/3 window) excl", lowlook, 1717L, 250L)
qc("usable trials (expect ~+86)", nrow(usable), 6204L, 200L)    # +86 'didnt hear' we cant drop

## ---- Tier 2: accuracy (increase in target looking, pair-salience-corrected) ----
usable[, prop_tgt := t_tgt/(t_tgt+t_dis)]
## per image: prop looking at it when TARGET vs when DISTRACTOR, within a subj+age
img <- rbind(
  usable[, .(subj, age, pair, img=timg, look=prop_tgt,  role="target")],
  usable[, .(subj, age, pair, img=dimg, look=1-prop_tgt, role="distractor")])
inc <- img[, .(look=mean(look)), by=.(subj,age,pair,img,role)]
inc <- dcast(inc, subj+age+pair+img ~ role, value.var="look")[!is.na(target) & !is.na(distractor)]
inc[, increase := target - distractor]
cat("\n=== Tier 2: accuracy 'increase in target looking' ===\n")
cat(sprintf("  grand-mean increase = %.3f  (target ~0.048; PASS if in [0.03,0.07])  %s\n",
            mean(inc$increase), if (mean(inc$increase) %between% c(0.03,0.07)) "PASS" else "CHECK"))
cat("  by age (mean increase; >0 & rising = good):\n")
print(inc[, .(increase=round(mean(increase),3), n=.N), by=age][order(age)])

## ---- RT extraction (the deliverable) ----
fr <- merge(fix, usable[, .(sess,trial,onset,age,subj)], by=c("sess","trial"))
rt <- fr[order(sess,trial,fstart), {
  o <- onset[1]
  ini <- aoi[fstart<=o & fend>o]; if(!length(ini)){b<-which(fstart<=o); ini<-if(length(b)) aoi[tail(b,1)] else NA_character_}
  tgt <- fstart[aoi=="TARGET" & fstart>o]
  .(age=age[1], subj=subj[1], ini=ini[1], rt=if(length(tgt)) min(tgt)-o else NA_real_)
}, by=.(sess,trial)]
rt[, valid := !is.na(ini) & ini=="DISTRACTOR" & !is.na(rt) & rt>=RT_LO & rt<=RT_HI]
v <- rt[valid==TRUE]
cat("\n=== Tier 3: RT (deliverable) ===\n")
cat(sprintf("  valid RTs=%d / %d kids; median=%.0f ms; by age:\n", nrow(v), uniqueN(v$subj), median(v$rt)))
print(v[, .(n=.N, med=round(median(rt))), by=age][order(age)])

## ---- emit lwl-schema RT (per-trial; child-level use) ----
out <- v[, .(dataset_name="seedlings_zhu", lab_subject_id=subj, lwl_age=age,
             lwl_log_rt=log(rt))]
fwrite(out, here("data/seedlings/seedlings_lwl_rt.csv"))
cat(sprintf("\nwrote data/seedlings/seedlings_lwl_rt.csv: %d RT rows, %d kids, ages %d-%d\n",
            nrow(out), uniqueN(out$lab_subject_id), min(out$lwl_age), max(out$lwl_age)))
