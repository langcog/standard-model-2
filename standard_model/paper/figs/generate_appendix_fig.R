require(ggplot2)
require(ggpubr)

# ability distro
ggplot(NULL, aes(c(-3,3))) +
  geom_area(stat = "function", fun = dnorm, fill = "grey80", xlim = c(-3, 3),
            alpha=.5) +
  geom_area(stat = "function", fun = dnorm, fill = "black", xlim = c(.5, .55)) +
  labs(x = "", y = "") + 
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = NULL) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))


# ODS and CDS rate distributions
p1 <- ggplot(NULL, aes(c(-3,3))) +
  geom_area(stat = "function", fun = dnorm, fill = "blue", xlim = c(-3, 3),
            alpha=.7) +
  geom_area(stat = "function", fun = dnorm, fill = "black", xlim = c(1.5, 1.55)) +
  labs(x = "", y = "") + ggtitle("ODS Rate Distribution") + 
  scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) 
#geom_vline(aes(xintercept=0), linetype="dashed")

p2 <- ggplot(NULL, aes(c(-3,3))) +
  geom_area(stat = "function", fun = dnorm, fill = "red", xlim = c(-3, 3),
            alpha=.7) +
  geom_area(stat = "function", fun = dnorm, fill = "black", xlim = c(-1.0, -1.05)) +
  labs(x = "", y = "") + ggtitle("CDS Rate Distribution") + 
  scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) 

ggarrange(p1, p2, nrow=2)
ggsave("CDS_ODS_distros_standard_model.pdf", width=3, height=5)


# per word pie charts

# Basic piechart
make_pie <- function(propCDS, propODS, title="") {
  prop_unlearned = 100 - propCDS - propODS
  dat <- data.frame(
    group=LETTERS[1:3],
    value=c(propCDS,propODS,prop_unlearned)
  )
  p <- ggplot(dat, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=.1, color="white") +
    coord_polar("y", start=0) + ggtitle(title) +
    scale_fill_manual(values=c("red", "blue", "grey80")) +
    theme_void() + # remove background, grid, numeric labels
    theme(legend.position="none")
}

# CHILDES / SUBTLEX frequency per million tokens:
# you - 44502 / 41857
# go - 6631 / 3793
# ball - 5774 / 105
# have - 5330 / 6161
# dog - 4042 / 193
# table - 3746 / 106
# tomorrow - 240 / 336

ggarrange(make_pie(70,20),
          make_pie(20,60),
          make_pie(33,33),
          make_pie(20,20),
          make_pie(10,10),
          nrow=5)

ggarrange(make_pie(42.5,40, title="you"), # 1.06x as frequent in CDS
          make_pie(35,20, title="go"), # 1.75x as freq in CDS
          make_pie(82.5,1.5, title="ball"), # 55x as freq in CDS
          make_pie(17.4,20, title="have"), # .87x as freq in CDS (1.15x as freq in adult speech)
          make_pie(63,3, title="dog"), # 21x as freq in CDS
          make_pie(7.1,10, title="tomorrow"), # .71x as freq in CDS (1.4x as freq in adult)
          nrow=6)
ggsave("standard_model_words_example.pdf", height=6, width=1)