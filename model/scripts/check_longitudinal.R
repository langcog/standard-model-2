## Quick check: does the Wordbank English preprocessed data contain any
## children measured at multiple ages? Lists per-person age counts and
## reports how many children would be usable for longitudinal fits.

source("model/R/config.R")
source("model/R/helpers.R")

df <- load_wordbank_data()

cat(sprintf("Total observations: %d\n", nrow(df)))
cat(sprintf("Unique children: %d\n", length(unique(df$person))))
cat(sprintf("Unique items: %d\n", length(unique(df$item))))
cat(sprintf("Age range: %.0f–%.0f months\n",
            min(df$age), max(df$age)))

# How many unique (person, age) combinations per person
per_person <- df %>%
  distinct(person, age) %>%
  count(person, name = "n_ages")

tbl <- per_person %>%
  count(n_ages, name = "n_children") %>%
  arrange(n_ages)
cat("\nDistribution of #ages per child:\n")
print(tbl, n = Inf)

n_long <- sum(per_person$n_ages >= 2)
cat(sprintf("\nChildren with >=2 age measurements: %d (%.1f%%)\n",
            n_long, 100 * n_long / nrow(per_person)))

if (n_long > 0) {
  # Show some examples
  long_ids <- per_person$person[per_person$n_ages >= 2]
  ex <- df %>% filter(person %in% long_ids[1:min(5, length(long_ids))]) %>%
    distinct(person, age) %>% arrange(person, age)
  cat("\nExample longitudinal children (first 5):\n")
  print(ex, n = Inf)

  # Age gaps between observations
  gaps <- df %>% filter(person %in% long_ids) %>%
    distinct(person, age) %>%
    arrange(person, age) %>%
    group_by(person) %>%
    summarise(min_gap = min(diff(age)),
              max_gap = max(diff(age)),
              span = max(age) - min(age),
              .groups = "drop")
  cat(sprintf("\nAge gaps among longitudinal kids:\n"))
  cat(sprintf("  median min gap: %.1f mo\n", median(gaps$min_gap)))
  cat(sprintf("  median span:    %.1f mo\n", median(gaps$span)))
  cat(sprintf("  max span:       %.1f mo\n", max(gaps$span)))
}
