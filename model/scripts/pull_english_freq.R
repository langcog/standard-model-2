## Compute English (American) CHILDES word frequencies. Parallel
## to pull_norwegian_freq.R: pulls all adult-speaker tokens from
## CHILDES English transcripts via childesr, counts each lowercased
## token, normalizes by total token count to get prob per word.
##
## Replaces the legacy `engWS_preprocessed.Rdata` source from the
## predecessor standard_model project, so English and Norwegian
## frequencies now come from the same reproducible pipeline.
##
## Writes: fits/english_word_freq.rds with columns (w, count, prob).
## Downstream: prepare_longitudinal_data.R joins this to Wordbank items
## via the same normalize-and-match logic as the Norwegian pipeline.

suppressPackageStartupMessages({
  library(childesr)
  library(dplyr)
})

source("model/R/config.R")

# Stick to North American English to match Wordbank's "English (American)"
# instrument. CHILDES `corpus_country` is the standard filter; the
# alternative `collection` ("Eng-NA" vs "Eng-UK") is less reliable
# across corpora.
message("Pulling English (American) tokens from CHILDES...")
message("  (large query -- this can take several minutes)")
tokens <- get_tokens(language = "eng",
                     collection = "Eng-NA",
                     role_exclude = "Target_Child",
                     token = "*")
message(sprintf("  got %d tokens across %d corpora",
                nrow(tokens),
                length(unique(tokens$corpus_id))))

freq <- tokens %>%
  mutate(w = tolower(gloss)) %>%
  filter(!is.na(w), w != "", !grepl("^[[:punct:]]+$", w)) %>%
  count(w, name = "count") %>%
  arrange(desc(count))

total <- sum(freq$count)
freq <- freq %>% mutate(prob = count / total)
message(sprintf("  %d unique word types, total tokens %d",
                nrow(freq), total))

out <- file.path(PATHS$fits_dir, "english_word_freq.rds")
saveRDS(freq, out)
cat(sprintf("\nSaved %s (%.1f MB)\n", out,
            file.info(out)$size / 1e6))

cat("\nTop 20 words:\n")
print(head(freq, 20))

cat("\nSanity check: prob for some common items (lowercased):\n")
anchors <- c("the", "and", "a", "ball", "dog", "mommy", "milk",
             "want", "more", "no", "yes", "bye")
print(freq %>% filter(w %in% anchors) %>% arrange(desc(prob)))
