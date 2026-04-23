require(childesr)
require(tidyverse)
require(ggrepel)
require(here)

wf <- read_csv(here("scripts/data/childes_english.csv")) %>%
  select(-final_count, -solo_count, -mean_sent_length) %>% # 45598
  filter(word_count>1, !is.na(word)) %>%# 27736 (21441 > twice, 14376 >5)
  arrange(desc(word_count)) %>%
  filter(!is.element(word, c("www", "w", "r", "e", "g", "y"))) # frequent non-words
# 10236 words

# problem: Mommy and mommy are different; should combine some of these counts
combine_word_counts <- function(wf, sources, target) {
  targ_idx = which(wf$word==target)
  if(length(targ_idx)==0) {
    print(paste("Target (",target,") not found in dataframe"))
    return(wf)
  }
  for(source in sources) {
    targ_idx = which(wf$word==target) # need to do this again in case it changes after removing source word
    source_idx = which(wf$word==source)
    if(length(source_idx)==0) {
      print(paste("Source (",source,") not found in dataframe"))
    } else {
      print(paste("source count:",wf[source_idx,]$word_count, 
                  "target count:",wf[targ_idx,]$word_count))
      wf[targ_idx,]$word_count = wf[targ_idx,]$word_count + wf[source_idx,]$word_count
      wf = wf[-source_idx,]
    }
  }
  return(wf)
}

wf = combine_word_counts(wf, c("Mommy","mom","Mom","mama","momma","Mama","Momma","Mommy's","mother"), "mommy") # 
wf = combine_word_counts(wf, c("Daddy", "Dad", "dad", "dada", "papa","Papa","Dada","father", "Daddy's"), "daddy")
wf = combine_word_counts(wf, c("Grandpa","grandfather","Grandfather","granddaddy"), "grandpa")
wf = combine_word_counts(wf, c("Grandma","grandmother","Grandmother","granny"), "grandma")
wf = combine_word_counts(wf, "have_to", "hafta")
wf = combine_word_counts(wf, "all_gone", "gone")
wf = combine_word_counts(wf, "in", "inside")
wf = combine_word_counts(wf, "kleenex", "tissue")
wf = combine_word_counts(wf, "got_to", "gotta")
wf = combine_word_counts(wf, "going", "gonna")
wf = combine_word_counts(wf, c("lots_of","a_lot","a_lot_of"), "lot")
wf = combine_word_counts(wf, "let", "lemme")
wf = combine_word_counts(wf, c("shush","hush"), "shh")
wf = combine_word_counts(wf, "thank_you", "thanks")
wf = combine_word_counts(wf, "night_night", "night")
wf = combine_word_counts(wf, "quack_quack", "quack")
wf = combine_word_counts(wf, "woof_woof", "woof")
wf = combine_word_counts(wf, "choo_choo", "choo")
wf = rbind(wf, cbind(word="child's own name", word_count=0)) # add this
wf$word_count = as.numeric(wf$word_count)
wf = combine_word_counts(wf, c("Matthew","Adam","Alex","Abe","Sarah","William","Naima",
                               "Eric","Melissa","Michael","Mark","Laura","Ross","Henry"), 
                         "child's own name") # "Child" ?

sum(wf$word_count) # 7,576,872

# maybe convert word_count to count_per_million tokens?
wf <- wf %>% mutate(word_count_norm = 1000000 * word_count / sum(word_count))

# 1 wordbank item is very infrequent: TV appears only 3x

# also want to get CHILDES freq per wordbank item
load("data/wordbank_eng_ws_wg_webcdi31-36mos.Rds")

relabel_col <- function(d_mat, old_name, new_name) {
  colnames(d_mat)[which(colnames(d_mat)==old_name)] = new_name
  return(d_mat)
}

# 594 of the 680 CDI items easily identified in CHILDES, but we'll recode others:

# for polysemous and multiword CDI phrases, we'll replace with the lower frequency word
d_mat = relabel_col(d_mat, "baa baa", "baa")
d_mat = relabel_col(d_mat, "woof woof", "woof")
d_mat = relabel_col(d_mat, "night night", "night")
d_mat = relabel_col(d_mat, "choo choo", "choo")
d_mat = relabel_col(d_mat, "yum yum", "yum")
d_mat = relabel_col(d_mat, "quack quack", "quack")

d_mat = relabel_col(d_mat, "fish (animal)", "fish")
d_mat = relabel_col(d_mat, "fish (food)", "fish")
d_mat = relabel_col(d_mat, "clean (action)", "clean")
d_mat = relabel_col(d_mat, "clean (description)", "clean")
d_mat = relabel_col(d_mat, "play dough", "playdough")
d_mat = relabel_col(d_mat, "potato chip", "chip") # chip 176, potato 424
d_mat = relabel_col(d_mat, "watch (action)", "watch")
d_mat = relabel_col(d_mat, "watch (object)", "watch")
d_mat = relabel_col(d_mat, "water (beverage)", "water")
d_mat = relabel_col(d_mat, "water (not beverage)", "water")
d_mat = relabel_col(d_mat, "work (action)", "work")
d_mat = relabel_col(d_mat, "work (place)", "work")
d_mat = relabel_col(d_mat, "little (description)", "little")
d_mat = relabel_col(d_mat, "lemme/let me", "lemme")
d_mat = relabel_col(d_mat, "need/need to", "need")
d_mat = relabel_col(d_mat, "inside/in", "inside")
d_mat = relabel_col(d_mat, "hafta/have to", "hafta")
d_mat = relabel_col(d_mat, "gonna/going to", "gonna")
d_mat = relabel_col(d_mat, "washing machine", "washing")
d_mat = relabel_col(d_mat, "wanna/want to", "wanna")
d_mat = relabel_col(d_mat, "all gone", "gone")
d_mat = relabel_col(d_mat, "try/try to", "try")
d_mat = relabel_col(d_mat, "dry (description)", "dry")
d_mat = relabel_col(d_mat, "dry (action)", "dry")
d_mat = relabel_col(d_mat, "orange (food)", "orange")
d_mat = relabel_col(d_mat, "orange (description)", "orange")
d_mat = relabel_col(d_mat, "tissue/kleenex", "tissue")
d_mat = relabel_col(d_mat, "on top of", "top")
d_mat = relabel_col(d_mat, "uh-oh", "uhoh")
d_mat = relabel_col(d_mat, "slide (action)", "slide")
d_mat = relabel_col(d_mat, "slide (object)", "slide")
d_mat = relabel_col(d_mat, "turn around", "around")
d_mat = relabel_col(d_mat, "dress (object)", "dress")
d_mat = relabel_col(d_mat, "can (object)", "can")
d_mat = relabel_col(d_mat, "can (auxiliary)", "can")
d_mat = relabel_col(d_mat, "drink (action)", "drink")
d_mat = relabel_col(d_mat, "drink (beverage)", "drink")
d_mat = relabel_col(d_mat, "gotta/got to", "gotta")
d_mat = relabel_col(d_mat, "give me five!", "five") # 
d_mat = relabel_col(d_mat, "uh oh", "uhoh")
d_mat = relabel_col(d_mat, "call (on phone)", "call")
d_mat = relabel_col(d_mat, "chicken (animal)", "chicken")
d_mat = relabel_col(d_mat, "chicken (food)", "chicken")
d_mat = relabel_col(d_mat, "green beans", "beans")
d_mat = relabel_col(d_mat, "swing (object)", "swing")
d_mat = relabel_col(d_mat, "swing (action)", "swing")
d_mat = relabel_col(d_mat, "so big!", "big")
d_mat = relabel_col(d_mat, "toy (object)", "toy")
d_mat = relabel_col(d_mat, "this little piggy", "piggy")
d_mat = relabel_col(d_mat, "rocking chair", "rocking")
d_mat = relabel_col(d_mat, "gas station", "station")
d_mat = relabel_col(d_mat, "lawn mower", "mower")
d_mat = relabel_col(d_mat, "a lot", "lot")
d_mat = relabel_col(d_mat, "did/did ya", "did")
d_mat = relabel_col(d_mat, "high chair", "chair")
d_mat = relabel_col(d_mat, "french fries", "fries")
d_mat = relabel_col(d_mat, "owie/boo boo", "owie")
d_mat = relabel_col(d_mat, "go potty", "potty")
d_mat = relabel_col(d_mat, "soda/pop", "soda") # soda: 318, pop: 1284
d_mat = relabel_col(d_mat, "next to", "next")
d_mat = relabel_col(d_mat, "shh/shush/hush", "shh")
d_mat = relabel_col(d_mat, "thank you", "thanks")
d_mat = relabel_col(d_mat, "ice cream", "cream") # cream: 1468, ice: 1529
d_mat = relabel_col(d_mat, "living room", "room") # living: 677 room: 3901
d_mat = relabel_col(d_mat, "belly button", "belly") # belly: 590 button: 1158
d_mat = relabel_col(d_mat, "buttocks/bottom*", "bottom") # buttocks: 0
d_mat = relabel_col(d_mat, "play pen", "pen") # pen: 918 play: 13999
d_mat = relabel_col(d_mat, "TV", "television")
d_mat = relabel_col(d_mat, "teddybear", "teddy") # teddy: 652, bear: 3766

d_mat = relabel_col(d_mat, "grandpa*", "grandpa")
d_mat = relabel_col(d_mat, "grandma*", "grandma")
d_mat = relabel_col(d_mat, "mommy*", "mommy")
d_mat = relabel_col(d_mat, "daddy*", "daddy")
d_mat = relabel_col(d_mat, "vagina*", "vagina")
d_mat = relabel_col(d_mat, "penis*", "penis")
d_mat = relabel_col(d_mat, "church*", "church")

d_mat = relabel_col(d_mat, "peanut butter", "peanut") # peanut: 262 butter: 1115
d_mat = relabel_col(d_mat, "gonna get you!", "get") # slightly unfair to replace 3-word phrase with 1 HF word

#d_mat = relabel_col(d_mat, "child's own name", "xxx") # yes?

words = intersect(wf$word, colnames(d_mat)) # now we're up to 655 / 680

length(unique(colnames(d_mat))) # 659 (4 duplicate?)

relabel_these = setdiff(colnames(d_mat), words)
# "babysitter's name" "child's own name"  
# "cockadoodledoo" "pet's name"

#  make wordbank deciles for each word and save those

# save CHILDES wordfreqs
wf$on_cdi = ifelse(is.element(wf$word, words), 1, 0)

wf <- wf %>% arrange(desc(word_count))
wf$prob = wf$word_count / sum(wf$word_count)
write.csv(wf, file=here("scripts/data/childes_english_word_freq_cleaned_noHapaxes.csv"))


plot_childes_word_freq_range <- function(wf, idx) {
  wf$label = ""
  wf$label[idx] = wf$word[idx]
  wf_plot = wf[min(idx):max(idx),]
  prop_CDI = paste("Proportion on CDI:", sum(wf_plot$on_cdi) / nrow(wf_plot))
  g <- ggplot(wf_plot, aes(x=min(idx):max(idx), y=word_count_norm)) + 
    geom_point(color = ifelse(wf_plot$on_cdi == 1, "red", "black")) + 
    geom_label_repel(aes(label=label), ylim=c(min(wf_plot$word_count_norm),NA)) + 
    xlab("Word Frequency Rank") + ylab("Word Frequency in CHILDES (per million tokens)") +
    annotate("text", x = median(idx), y = max(wf_plot$word_count_norm), label = prop_CDI, color="red")
  print(prop_CDI)
  ggsave(paste0("figures/CHILDES_word_freq_", min(idx), "-", max(idx), ".pdf"), width=7.5, height=5)
}

#plot_childes_word_freq_range(wf, c(1,2,3,4,5,6, seq(8,100,4))) # .79
#plot_childes_word_freq_range(wf, 100+c(1,seq(4,100,4))) # .49
#plot_childes_word_freq_range(wf, 200+c(1,2,3,seq(4,100,4))) # .61
#plot_childes_word_freq_range(wf, 300+c(1,3,seq(4,100,4))) # .44
#plot_childes_word_freq_range(wf, 400+c(1,3,seq(4,100,4))) # .45
#plot_childes_word_freq_range(wf, 500+c(1,3,seq(4,100,4))) # .35
#plot_childes_word_freq_range(wf, 600+c(1,3,seq(4,100,4))) # .35
#plot_childes_word_freq_range(wf, 700+c(1,3,seq(4,100,4))) # .3
#plot_childes_word_freq_range(wf, 800+c(1,3,seq(4,100,4))) # .33
#plot_childes_word_freq_range(wf, 900+c(1,3,seq(4,100,4))) # .34
