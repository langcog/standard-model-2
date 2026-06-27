#!/usr/bin/env Rscript
# Analyze novel-word (nonce) acquisition from the warm-start continuations.
# Reads /tmp/nonce_pull/runs/*/nonce.csv + endpoints.csv (val_bpb,J per source checkpoint).
# Two tests:
#   A (regime): per-arm new-word acquisition vs schedule/J, controlling for starting val_bpb.
#   B (kappa):  within a run, are LATER-introduced cohorts acquired faster? (within-EXPOSURE-WINDOW
#               drop, to remove the forgetting confound that contaminates a naive first->last drop).
suppressWarnings(suppressMessages(library(ggplot2)))
OUT <- "/Users/mcfrank/Projects/standard_model_2/studies/llm/childes_flatlr"
PULL <- "/tmp/nonce_pull"

armmeta <- function(tag){  # cont_flat1e5_s42[_4x] -> sched/lr/seed/budget
  budget <- if (grepl("_4x$", tag)) "4x" else "1x"
  s <- sub("_4x$", "", sub("^cont_", "", tag))
  seed <- if (grepl("_s42$", s)) 42 else 1
  base <- sub("_s(42|1)$", "", s)
  sched <- if (base == "std") "decay" else "flat"
  lr <- if (base %in% c("std","flat1e4")) 1e-4 else if (base=="flat3e5") 3e-5 else 1e-5
  data.frame(tag=tag, sched=sched, lr=lr, seed=seed, budget=budget,
             armlab=if(base=="std") "decay 1e-4" else paste("flat", format(lr,scientific=TRUE)),
             stringsAsFactors=FALSE)
}

dirs <- list.dirs(file.path(PULL,"runs"), recursive=FALSE)
all <- do.call(rbind, lapply(dirs, function(d){
  f <- file.path(d,"nonce.csv"); if(!file.exists(f)) return(NULL)
  x <- read.csv(f); x$tag <- basename(d); x
}))
maxstep <- max(all$step)

# per (tag,nonce): pre = surprisal nearest its segment START; post = nearest segment END;
# within-window acquisition = pre - post. cohort c exposed in step-fraction [c/4,(c+1)/4].
nearest <- function(df, target) df$mean_span_nll[which.min(abs(df$step - target))]
rows <- list()
for(tg in unique(all$tag)) for(nz in unique(all$nonce)){
  d <- all[all$tag==tg & all$nonce==nz,]; if(!nrow(d)) next
  c <- d$cohort[1]; ms <- max(d$step)
  segs <- (c/4)*ms; sege <- ((c+1)/4)*ms
  pre <- nearest(d, segs); post <- nearest(d, sege); floor_ <- min(d$mean_span_nll)
  rows[[length(rows)+1]] <- data.frame(tag=tg, nonce=nz, cohort=c, intro_frac=c/4,
                                       pre=pre, post=post, win_drop=pre-post, floor=floor_)
}
ac <- do.call(rbind, rows)
ac <- merge(ac, do.call(rbind, lapply(unique(ac$tag), armmeta)), by="tag")

# ---- B: kappa test -- within-window acquisition vs cohort intro time ----
byc <- aggregate(win_drop~intro_frac, ac, mean)
pB <- ggplot(ac, aes(factor(intro_frac), win_drop)) +
  geom_boxplot(outlier.size=.6, fill="#cfe8f3") +
  stat_summary(fun=mean, geom="line", aes(group=1), colour="#d73027", linewidth=1) +
  stat_summary(fun=mean, geom="point", colour="#d73027", size=2) +
  labs(title="(B) Within-run kappa test: are later-introduced nonces acquired faster?",
       subtitle=paste0("within-EXPOSURE-WINDOW surprisal drop by cohort introduction time (forgetting-controlled); ",
                       "trend: ", paste(sprintf("%.2f",byc$win_drop), collapse=" -> ")),
       x="cohort introduction time (fraction of run)", y="within-window acquisition (nats)") +
  theme_minimal(base_size=12)
ggsave(file.path(OUT,"fig_nonce_kappa.png"), pB, width=8, height=4.6, dpi=130)

# ---- A: regime test -- per-arm acquisition vs J, controlling for val_bpb (all 12 checkpoints) ----
end <- read.csv(file.path(OUT,"endpoints.csv"))  # tag like std_s42_1x; has val_bpb,J,budget
arm <- aggregate(cbind(win_drop, floor)~tag, ac, mean)         # tag = continuation tag
arm$key <- sub("^cont_", "", arm$tag)
arm$key <- ifelse(grepl("_(1x|4x)$", arm$key), arm$key, paste0(arm$key, "_1x"))
arm <- merge(arm, end[,c("tag","val_bpb","J","sched","lr","seed","budget")],
             by.x="key", by.y="tag")
arm$armlab <- ifelse(arm$sched=="decay","decay 1e-4",
              ifelse(arm$lr==1e-4,"flat 1e-4", ifelse(arm$lr==3e-5,"flat 3e-5","flat 1e-5")))
write.csv(arm[order(arm$budget, arm$lr),], file.path(OUT,"nonce_by_arm.csv"), row.names=FALSE)

# PRIMARY metric = floor (nonce surprisal reached after the fixed ~100 exposures; lower = learned
# better). win_drop is confounded by the pre-level (overfit models start nonces very high -> big drop
# to a still-poor floor), so floor is the clean "how well did it learn the new word" measure.
pal <- c("decay 1e-4"="#d73027","flat 1e-4"="#fc8d59","flat 3e-5"="#91bfdb","flat 1e-5"="#4575b4")
shp <- c("1x"=16, "4x"=17)
pA1 <- ggplot(arm, aes(J, floor)) +
  geom_smooth(method="lm", se=TRUE, colour="grey60", fill="grey90", linewidth=.5) +
  geom_point(aes(colour=armlab, shape=budget), size=3.2) +
  scale_colour_manual(values=pal, name=NULL) + scale_shape_manual(values=shp, name="budget") +
  scale_y_reverse() +
  labs(title="(A) New-word mastery tracks the REGIME (J)",
       x="J of source checkpoint (returns signature)",
       y="nonce floor surprisal  (axis reversed: up = better learned)") + theme_minimal(base_size=12)
pA2 <- ggplot(arm, aes(val_bpb, floor)) +
  geom_smooth(method="lm", se=TRUE, colour="grey60", fill="grey90", linewidth=.5) +
  geom_point(aes(colour=armlab, shape=budget), size=3.2) +
  scale_colour_manual(values=pal, name=NULL) + scale_shape_manual(values=shp, name="budget") +
  scale_y_reverse() +
  labs(title="(A, control) vs starting competence (val_bpb)",
       x="val_bpb of source checkpoint (lower = more competent)",
       y="nonce floor surprisal  (axis reversed: up = better learned)") + theme_minimal(base_size=12)
ggsave(file.path(OUT,"fig_nonce_regime_J.png"), pA1, width=6.5, height=4.4, dpi=130)
ggsave(file.path(OUT,"fig_nonce_regime_valbpb.png"), pA2, width=6.5, height=4.4, dpi=130)

cat("=== per-arm new-word acquisition (mean within-window drop), all 12 checkpoints ===\n")
print(arm[order(-arm$win_drop), c("armlab","seed","budget","val_bpb","J","win_drop","floor")], row.names=FALSE)
cat("\n=== B: within-window acquisition by cohort intro time ===\n"); print(byc, row.names=FALSE)
cat("\n=== PRIMARY: floor (nonce mastery; lower=better). lm(floor ~ J + val_bpb), n=", nrow(arm), " ===\n", sep="")
print(summary(lm(floor~J+val_bpb, arm))$coefficients)
cat("cor(floor, J)=", round(cor(arm$floor, arm$J),3),
    "  cor(floor, val_bpb)=", round(cor(arm$floor, arm$val_bpb),3), "\n", sep="")
cat("(negative J coef = higher-J regime -> lower floor -> better new-word learning, controlling competence)\n")
