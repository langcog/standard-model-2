#!/usr/bin/env Rscript
# Analyze the CHILDES flat-vs-decay LR sweep. Reads /tmp/flatlr_pull/runs/*/grammar.csv
# (+ J from /tmp/flatlr_pull/J.tsv) and writes figures to this script's dir.
suppressWarnings(suppressMessages({library(ggplot2)}))
OUT <- "/Users/mcfrank/Projects/standard_model_2/studies/llm/childes_flatlr"
PULL <- "/tmp/flatlr_pull"

meta <- function(tag) {
  sched <- if (grepl("^std", tag)) "decay" else "flat"
  lr <- if (grepl("^std|flat1e4", tag)) 1e-4 else if (grepl("flat3e5", tag)) 3e-5 else if (grepl("flat1e5", tag)) 1e-5 else NA
  seed <- if (grepl("_s42_", tag)) 42 else 1
  budget <- if (grepl("_4x$", tag)) "4x" else "1x"
  data.frame(tag = tag, sched = sched, lr = lr, seed = seed, budget = budget,
             arm = paste0(sched, "_", format(lr, scientific = TRUE)), stringsAsFactors = FALSE)
}

dirs <- list.dirs(file.path(PULL, "runs"), recursive = FALSE)
dirs <- dirs[!grepl("smoke", dirs)]
traj <- do.call(rbind, lapply(dirs, function(d) {
  tag <- basename(d); g <- read.csv(file.path(d, "grammar.csv"))
  g <- g[!is.na(g$val_bpb), ]; cbind(meta(tag), g)
}))
traj$armlab <- ifelse(traj$sched == "decay", "decay 1e-4",
               ifelse(traj$lr == 1e-4, "flat 1e-4",
               ifelse(traj$lr == 3e-5, "flat 3e-5", "flat 1e-5")))
traj$armlab <- factor(traj$armlab, levels = c("decay 1e-4","flat 1e-4","flat 3e-5","flat 1e-5"))
pal <- c("decay 1e-4"="#d73027","flat 1e-4"="#fc8d59","flat 3e-5"="#91bfdb","flat 1e-5"="#4575b4")

# ---- endpoint table (+J) ----
J <- read.delim(file.path(PULL, "J.tsv")); names(J)[1] <- "tag"
ends <- do.call(rbind, lapply(dirs, function(d){
  tag<-basename(d); g<-read.csv(file.path(d,"grammar.csv")); g<-g[!is.na(g$val_bpb),]
  fin<-g[nrow(g),]; m<-meta(tag)
  data.frame(m, val_bpb=fin$val_bpb, zorro=fin$zorro, blimp=fin$blimp,
             J=J$J[match(tag,J$tag)])
}))
write.csv(ends[order(ends$budget,ends$sched,ends$lr),], file.path(OUT,"endpoints.csv"), row.names=FALSE)

# ---- (1) val_bpb trajectory: 1x, all arms (dose-response of overfitting) ----
t1 <- traj[traj$budget=="1x",]
p1 <- ggplot(t1, aes(epoch, val_bpb, colour=armlab, group=interaction(tag))) +
  geom_line(alpha=.85) + geom_point(size=.7, alpha=.6) +
  scale_x_log10() + scale_colour_manual(values=pal, name=NULL) +
  labs(title="CHILDES val loss vs training (1x = 20 epochs)",
       subtitle="lower flat LR -> lower (better) val_bpb: flat 1e-5 < flat 3e-5 < decay < flat 1e-4 (overfits)",
       x="epoch (log)", y="val bits/token") + theme_minimal(base_size=12)
ggsave(file.path(OUT,"fig_valbpb_1x.png"), p1, width=8, height=4.6, dpi=130)

# ---- (2) the overfitting flip: val_bpb 1x vs 4x, decay vs flat 1e-5 ----
key <- traj[traj$armlab %in% c("decay 1e-4","flat 1e-5"),]
p2 <- ggplot(key, aes(epoch, val_bpb, colour=armlab, linetype=budget, group=tag)) +
  geom_line(linewidth=.8) +
  scale_x_log10() + scale_colour_manual(values=pal, name=NULL) +
  scale_linetype_manual(values=c("1x"="solid","4x"="22"), name="budget") +
  labs(title="Extended training (4x=60ep): standard decay OVERFITS, flat-low is robust",
       subtitle="decay val_bpb 2.06 -> 2.67 with more epochs; flat 1e-5 stays ~2.0-2.1",
       x="epoch (log)", y="val bits/token") + theme_minimal(base_size=12)
ggsave(file.path(OUT,"fig_overfit_4x.png"), p2, width=8, height=4.6, dpi=130)

# ---- (3) Zorro (primary grammar) trajectory: grammar ~flat across schedules ----
p3 <- ggplot(traj, aes(epoch, zorro, colour=armlab, linetype=budget, group=tag)) +
  geom_line(alpha=.85) + geom_hline(yintercept=.5, linetype=3, colour="grey60") +
  scale_x_log10() + scale_colour_manual(values=pal, name=NULL) +
  scale_linetype_manual(values=c("1x"="solid","4x"="22"), name="budget") +
  labs(title="Zorro (CHILDES-vocab grammar) vs training — primary grammar metric",
       subtitle="grammar plateaus ~0.68-0.70 across ALL schedules & budgets: decoupled from the LR/loss story",
       x="epoch (log)", y="Zorro accuracy") + theme_minimal(base_size=12) +
  coord_cartesian(ylim=c(.48,.78))
ggsave(file.path(OUT,"fig_zorro.png"), p3, width=8, height=4.6, dpi=130)

# ---- (4) summary: J dose-response (1x) + endpoint val_bpb/J decay-vs-flatlow 1x/4x ----
e1 <- ends[ends$budget=="1x",]; e1$armlab <- factor(
  ifelse(e1$sched=="decay","decay 1e-4", ifelse(e1$lr==1e-4,"flat 1e-4", ifelse(e1$lr==3e-5,"flat 3e-5","flat 1e-5"))),
  levels=c("decay 1e-4","flat 1e-4","flat 3e-5","flat 1e-5"))
ag <- aggregate(J~armlab, e1, mean)
p4 <- ggplot(ag, aes(armlab, J, fill=armlab)) + geom_col(width=.7) +
  geom_hline(yintercept=0, linetype=2, colour="grey50") +
  scale_fill_manual(values=pal, guide="none") +
  labs(title="CDI returns-to-experience J (1x): the flat-low-LR signature replicates on CHILDES",
       subtitle="lower flat LR -> higher J (more increasing returns); decay & flat-full are lowest",
       x=NULL, y="J = late - early gain-rate") + theme_minimal(base_size=12)
ggsave(file.path(OUT,"fig_J_dose.png"), p4, width=8, height=4.2, dpi=130)

cat("figures + endpoints.csv written to", OUT, "\n")
print(ends[order(ends$budget,ends$sched,ends$lr), c("tag","sched","lr","seed","budget","val_bpb","zorro","blimp","J")], row.names=FALSE)
