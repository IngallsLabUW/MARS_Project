library(tidyverse)

# Concept:
# Arrange compounds by mz and assign experiment ID + number to form primary_key.
# Standardize columns to reflect MARS protocolsj (all lowercase, underscore spacing, abbreviated).
untargeted.metco <- read.csv("data_extra/SuppTable3.csv", stringsAsFactors = FALSE) %>%
  select(1:6) %>%
  arrange(mz) %>%
  separate(Fraction, into = c("column", "charge")) %>%
  mutate(primary_key = paste("MetCo", 1:nrow(.), sep = ""),
         z = ifelse(charge == "Negative", -1, 1)) %>%
  rename(previous_match = BestMatch,
         confidence_rank = Confidence) %>%
  select(primary_key, mz, rt, column, z, previous_match, confidence_rank, -charge, -MF)


