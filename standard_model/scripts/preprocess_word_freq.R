wf = read.csv("data/childes_english_word_freq_cleaned_noHapaxes.csv")
load("data/eng_ws_raw_data.Rds")

d <- eng_ws %>%
  mutate(produces = value == "produces") %>%
  filter(!is.na(category)) %>% 
  mutate(produces1 = if_else(produces==TRUE,1,0)) %>%
  select(data_id, produces1, age, production, sex, definition, lexical_class) %>%
  rename(produces = produces1) %>%
  rename(item = definition) %>%
  rename(person = data_id)

subjs_missing_data = unique(d[which(is.na(d$produces)),]$person)
# get rid of 28 subjects with missing responses
d <- subset(d, !is.element(person, subjs_missing_data))

# relabel CDI items to match CHILDES vocab
wf$word = as.character(wf$word)
wf$word[wf$word == "chicken"] <- "chicken (animal)"
chicken <- c(10183,"chicken2", 1893, 247.5778114, 1) #merge wf
wf<-rbind(wf,chicken)
wf$word[wf$word == "chicken2"] <- "chicken (food)"
wf$word[wf$word == "drink"] <- "drink (beverage)"
wf$word[wf$word == "fish"] <- "fish (animal)"
fish <- c(10183,"fish2", 2548, 333.2426115, 1) #merge wf
wf<-rbind(wf,fish)
wf$word[wf$word == "fish2"] <- "fish (food)"
wf$word[wf$word == "toy"] <- "toy (object)"
wf$word[wf$word == "ice+cream" ] <- "ice cream"	
wf$word[wf$word == "potato+chip" ] <- "potato chip"	
wf$word[wf$word == "dress" ] <- "dress (object)"	
wf$word[wf$word == "can" ] <- "can (object)"	
owie <- c(10183,"owie/boo boo", 245, 32.04256, 1) #merge wf
wf<-rbind(wf,owie)
highchair <- c(10182, "high chair", 191, 24.98012, 1) #merge wf
wf<-rbind(highchair, wf)
wf$word[wf$word == "rocking+chair" ] <- "rocking chair"
wf$word[wf$word == "lawnmower" ] <- "lawn mower"
wf$word[wf$word == "water" ] <- "water (not beverage)"
wf$word[wf$word == "teddy+bear" ] <- "teddybear"
wf$word[wf$word == "fries" ] <- "french fries"
wf$word[wf$word == "orange" ] <- "orange (food)"
soda <- c(10180, "soda/pop", 1602, 209.5191, 1) #merge wf
wf<-rbind(soda, wf)
wf$word[wf$word == "belly+button" ] <- "belly button"
wf$word[wf$word == "penis" ] <- "penis*"
tissue <- c(10181, "tissue/kleenex", 497, 65.00062, 1) #merge wf
wf<-rbind(tissue, wf)
wf$word[wf$word == "living+room" ] <- "living room"
wf$word[wf$word == "slide" ] <- "slide (object)"
wf$word[wf$word == "play+dough" ] <- "play dough"
wf$word[wf$word == "peanut+butter" ] <- "peanut butter"
wf$word[wf$word == "water" ] <- "water (beverage)"
water <- c(10181, "water2", 5874, 768.2366954, 1) #merge wf
wf<-rbind(water, wf)
wf$word[wf$word == "water2" ] <- "water (beverage)"
wf$word[wf$word == "water" ] <- "water (beverage)"
wf$word[wf$word == "bottom" ] <- "buttocks/bottom*" 
wf$word[wf$word == "vagina" ] <- "vagina*"
wf$word[wf$word == "watch" ] <- "watch (object)"
wf$word[wf$word == "play+pen" ] <- "play pen"
wf$word[wf$word == "baa" ] <- "baa baa"
choo <- c(10180, "choo choo", 1113, 145.5648, 1) #merge wf
wf<-rbind(choo, wf)
quacks <- c(10180, "quack quack", 1145, 149.7499, 1)
wf<-rbind(quacks, wf)
wf$word[wf$word == "uhoh" ] <- "uh oh"
woof <- c(10180, "woof woof", 730, 95.47375, 1) #merge wf
wf<-rbind(woof, wf)
yum <- c(10180, "yum yum", 2005, 144.5185, 1) #merge wf
wf<-rbind(yum, wf)
wf$word[wf$word == "church" ] <- "church*"
wf$word[wf$word == "daddy" ] <- "daddy*"
wf$word[wf$word == "grandma" ] <- "grandma*"
wf$word[wf$word == "grandpa" ] <- "grandpa*"
wf$word[wf$word == "mommy" ] <- "mommy*"
wf$word[wf$word == "call" ] <- "call (on phone)"
wf$word[wf$word == "night+night" ] <- "night night"
wf$word[wf$word == "shh" ] <- "shh/shush/hush"
wf$word[wf$word == "thankyou" ] <- "thank you"
wf$word[wf$word == "clean" ] <- "clean (action)"
drink <- c(10183,"drink2", 3353, 438.5253047, 1) #merge wf
wf<-rbind(wf,drink)
wf$word[wf$word == "drink2" ] <- "drink (action)"
wf$word[wf$word == "dry" ] <- "dry (action)"
slide <- c(10183,"slide2", 1240, 162.1745833, 1) #merge wf
wf<-rbind(wf,slide)
wf$word[wf$word == "slide2" ] <- "slide (action)"
wf$word[wf$word == "swing" ] <- "swing (action)"
watch <- c(10183,"watch2", 5337, 698.0046379, 1) #merge wf
wf<-rbind(wf,watch)
wf$word[wf$word == "watch2" ] <- "watch (action)"
wf$word[wf$word == "work" ] <- "work (action)"
clean <- c(10183,"clean2", 2614, 341.8744845, 1) #merge wf
wf<-rbind(wf,clean)
wf$word[wf$word == "clean2" ] <- "clean (description)"
dry <- c(10183,"dry2", 984, 128.6933790, 1) #merge wf
wf<-rbind(wf,dry)
wf$word[wf$word == "dry2" ] <- "dry (description)"
wf$word[wf$word == "little" ] <- "little (description)"
orange <- c(10183,"orange2", 2297, 300.4153370, 1) #merge wf
wf<-rbind(wf,orange)
wf$word[wf$word == "orange2" ] <- "orange (description)"
wf$word[wf$word == "inside" ] <- "inside/in"
wf$word[wf$word == "next_to" ] <- "next to"
can <- c(10183,"can2", 50448, 6597.8898209, 1) #merge wf
wf<-rbind(wf,can)
wf$word[wf$word == "can2" ] <- "can (auxiliary)"
wf$word[wf$word == "gotta" ] <- "gotta/got to"
wf$word[wf$word == "need" ] <- "need/need to"
wf$word[wf$word == "did" ] <- "did/did ya"
wf$word[wf$word == "hafta" ] <- "hafta/have to"
wf$word[wf$word == "try" ] <- "try/try to"
wf$word[wf$word == "gonna" ] <- "gonna/going to"
wf$word[wf$word == "lemme" ] <- "lemme/let me"
wf$word[wf$word == "wanna" ] <- "wanna/want to"
swing <- c(10183,"swing2", 637, 83.3106529, 1) #merge wf
wf<-rbind(wf,swing)
wf$word[wf$word == "swing2" ] <- "swing (object)"
work <- c(10183,"work2", 3762, 492.0167600, 1) #merge wf
wf<-rbind(wf,work)
wf$word[wf$word == "work2" ] <- "work (place)"

# check for inconsistencies in the word
# a list of words with no equivalent on wf: gas station; babysitter's name; pet's name; washing machine; give me five!; gonna get you!; go potty; so big!; this little piggy; turn around; "all gone";on top of; a lot; green beans

wf$word_count = as.numeric(wf$word_count)
wf$word_count_norm = as.numeric(wf$word_count_norm)

wf <- wf %>%
  rename(item = word) %>%
  mutate(prob = word_count_norm / sum(word_count_norm))

d_wf = d %>%
  left_join(wf, by="item")

wf_na = d_wf %>%
  filter(is.na(word_count))
unique(wf_na$item)

# assign minimum probability to all of the CDI items with missing frequencies

d_wf[which(is.na(d_wf$prob)),]$prob = min(d_wf$prob, na.rm=T)
# d_wf$prob = d_wf$prob/sum(d_wf$prob)
d_wf$produces = as.integer(d_wf$produces)

d_unique <- d_wf %>%
  select(item, prob)
d_unique<-d_unique[!duplicated(d_unique$item),]
sum(d_unique$prob)

save(d_wf, file="data/engWS_preprocessed.Rdata")