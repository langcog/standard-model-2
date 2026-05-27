## Extract longitudinal data for ALL 7 candidate languages.
## Just calls 01_extract_one.R in a loop. Each pull takes
## 10-60 sec depending on language size.

source("model/R/config.R")

LANGS <- c("English (American)", "Norwegian", "Finnish",
            "French (Quebecois)", "Japanese", "Spanish (Mexican)",
            "French (French)")

for (lang in LANGS) {
  cat("\n###########################################\n")
  cat(sprintf("### Extracting: %s\n", lang))
  cat("###########################################\n")
  cmd <- sprintf("Rscript model/scripts/glmer_ladder/01_extract_one.R \"%s\"", lang)
  status <- system(cmd)
  if (status != 0) {
    warning(sprintf("FAILED for %s (exit %d)", lang, status))
  }
}
cat("\n=== Done extracting all languages ===\n")
system("ls -la fits/glmer_ladder/data_*.rds")
