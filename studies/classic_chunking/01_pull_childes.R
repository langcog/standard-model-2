suppressPackageStartupMessages({library(childesr); library(dplyr); library(readr)})
x <- get_utterances(collection = "Eng-NA", role_exclude = "Target_Child")
x <- x %>% filter(!is.na(gloss), gloss != "", !is.na(target_child_age),
                  target_child_age >= 10, target_child_age <= 60) %>%
  select(id, transcript_id, corpus_name, utterance_order, speaker_role,
         target_child_age, gloss, num_tokens)
cat(nrow(x), "utterances;", n_distinct(x$corpus_name), "corpora\n")
write_csv(x, file.path(Sys.getenv("OUTDIR"), "childes_eng_na_utts.csv.gz"))
cat("done\n")
