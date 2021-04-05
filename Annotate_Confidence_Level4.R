library(fuzzyjoin)
library(tidyverse)

# Confidence Level 4: Everything else!

Confidence.Level.3 <- read.csv("data_processed/confidence_level3.csv")

Confidence.Level.4 <- Confidence.Level.3 %>%
  mutate(confidence_rank = ifelse(is.na(confidence_rank), 4, confidence_rank))