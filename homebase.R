library(fuzzyjoin)
library(tidyverse)


## Raafay: Don't worry about the last section of this script, or about the previous_rank/previous_match columns. 
# Start with homebase production.

# --------------
# Concept:
# Arrange compounds by mz and assign experiment ID + number to form primary_key.
# The primary_key column is a frankenstein of the paper title where the data comes from. For data not in papers, we may need 
# to use the cruise, the culture name, etc. It will end up being identified in a TBD metadata column. 
# Standardize columns to reflect MARS protocols (all lowercase, underscore spacing, abbreviated).

# Potential Questions:
# I have an unknown mass feature. Where else has it been seen in the lab? How similar are they? What is the metadata of both features?
# What is the overall quality of this feature? How are the peaks? 
# Is it a common contaminant? Does it show up on the CC list?

untargeted.metco <- read.csv("data_extra/MetCo.csv") %>% # From Katherine's "Metabolomic Consequences" paper
  select(MF_Frac, BestMatch, Confidence, mz, rt, MS2dat1) %>%
  arrange(mz) %>%
  separate(MF_Frac, into = c("MF", "columncharge"), sep = "_") %>%
  mutate(column = ifelse(str_detect(columncharge, "HILIC"), "HILIC", "RP"),
         z = ifelse(str_detect(columncharge, "Neg"), -1, 1),
         primary_key = paste("MetCo", 1:nrow(.), sep = "")) %>%
  rename(previous_match = BestMatch,
         previous_rank = Confidence) %>%
  select(primary_key, mz, rt, column, z, MS2dat1, previous_match, previous_rank) %>%
  rename(MS2 = MS2dat1) %>%
  mutate(Meta = "From Katherine Heal's 2019 'Metabolomic Consequences' paper, culture work") # or something like this

untargeted.sulfnet <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") %>% # From Katherine's "Sulfonate-Based Networks" paper
  separate(MassFeature_Column, into = c("MassFeature", "drop"), sep = "_") %>%
  rename(column = Column,
         previous_match = Identification,
         previous_rank = Confidence) %>%
  mutate(primary_key = paste("SulfNet", 1:nrow(.), sep = "")) %>%
  select(primary_key, mz, rt, column, z, MS2, previous_match, previous_rank) %>%
  mutate(Meta = "From Katherine Heal's 2019 'Sulfonate-based Networks' paper, collected in the NPSG") # and other metadata here or in another column

# Combine to potential homebase
homebase <- untargeted.metco %>%
  rbind(untargeted.sulfnet)


## Raafay you can ignore this for now, haven't yet thought about how to do comparisons. 
test <- read.csv("data_from_lab_members/final_peaks_Will.csv") %>%
  filter(feature == "FT0089") %>%
  select(feature, mz, rt) %>%
  mutate(compound.ID = 1)

t <- homebase %>%
  difference_left_join(test, by = c("mz"), max_dist = 0.02) %>%
  filter(!is.na(compound.ID))