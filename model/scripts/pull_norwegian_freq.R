## Compute Norwegian CHILDES word frequencies for joining to Wordbank
## Norwegian WS items.
##
## Uses childesr to pull all tokens from Norwegian (nor) transcripts,
## counts each lowercased token, normalizes by total token count to get
## prob per word, then tries to match Wordbank item_definitions.
##
## Writes: model/fits/norwegian_word_freq.rds with columns (item, count,
## prob).

suppressPackageStartupMessages({
  library(childesr)
  library(dplyr)
})

source("model/R/config.R")

message("Pulling Norwegian tokens from CHILDES...")
# Only adult speakers to match child-directed speech
tokens <- get_tokens(language = "nor", role_exclude = "Target_Child",
                     token = "*")
message(sprintf("  got %d tokens", nrow(tokens)))

freq <- tokens %>%
  mutate(w = tolower(gloss)) %>%
  filter(!is.na(w), w != "", !grepl("^[[:punct:]]+$", w)) %>%
  count(w, name = "count") %>%
  arrange(desc(count))

total <- sum(freq$count)
freq <- freq %>% mutate(prob = count / total)
message(sprintf("  %d unique word types, total tokens %d", nrow(freq), total))

saveRDS(freq, file.path(PATHS$fits_dir, "norwegian_word_freq.rds"))
cat("\nSaved model/fits/norwegian_word_freq.rds\n")
cat("Top 10 words:\n")
print(head(freq, 10))
