library(fuzzyjoin)
library(tidyverse)
options(scipen = 999)
options(digits = 6)

source("Functions.R")

# Notes from this step

# WHY ARE ROWS GETTING DROPPED IN FUZZYJOIN
# Timestamped/SHA-marked/Groundhog standards & MS2 sheets need to be used
# Can make a with/without MS2 total similarity score
# Why the rt bajillion zeros?
# Extra rows in mission.accomplished: why?
# What about multiple matches?

mz.flexibility <- 0.02
rt.flexibility <- 0.02 # seconds

# Experimental Values -----------------------------------------------------
# Must be in the following format, including capitalization:
# MassFeature (character), mz (numeric), rt (numeric), column (numeric), z (numeric), MS2 (in concatenated format, character)

Experimental.Values <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") %>%
  separate(MassFeature_Column, into = c("MassFeature", "drop"), sep = "_") %>%
  rename(column = Column,
         KRH_identification = Identification) %>%
  select(MassFeature, KRH_identification, mz, rt, column, z, MS2)

# Known Values -----------------------------------------------------

# MS2s
Cyano.MS2 <- read.csv("data_extra/Standards_MS2/Cyano_stds_withMS2s.csv") %>%
  mutate(column = "RP",
         z = NA)
HILICPos.MS2 <- read.csv("data_extra/Standards_MS2/HILICPos_stds_withMS2s.csv") %>%
  mutate(column = "HILIC",
         z = 1)
HILICNeg.MS2 <- read.csv("data_extra/Standards_MS2/HILICNeg_stds_withMS2s.csv") %>%
  mutate(column = "HILIC",
         z = -1)

# MS2s combined
Standards.MS2 <- Cyano.MS2 %>%
  rbind(HILICPos.MS2, HILICNeg.MS2) %>%
  rename(MS2 = MS2s) %>%
  select(-rtmin, -rtmax, -rt, -mz_obs, -mz)

# Standards csv
Standards <- read.csv("data_extra/Ingalls_Lab_Standards_Feb1.csv") %>%
  mutate(rt = (RT..min. * 60)) %>%
  rename(compound = Compound.Name_old,
         mz = m.z,
         column = Column) %>%
  select(compound, mz, rt, column, z)

# All known variables
Knowns <- Standards %>% 
  left_join(Standards.MS2, by = c("compound", "column", "z")) %>%
  transform(mz = as.numeric(mz))

# All unknown variables
Unknowns <- Experimental.Values %>% 
  mutate(compound_unknown = 1:nrow(Experimental.Values)) %>%
  select(compound_unknown, KRH_identification, mz, rt, column, z, MS2) 

## Confidence Level 1 Matching ----------------------------------------
My.Fuzzy.Join <- Knowns %>%
  difference_left_join(Unknowns, by = c("mz"), max_dist = mz.flexibility) %>% 
  rename(compound_known = compound,
         mz_known = mz.x,
         rt_seconds_known = rt.x,
         column_known = column.x,
         z_known = z.x,
         MS2_known = MS2.x,
         mz_unknown = mz.y,
         rt_seconds_unknown = rt.y,
         column_unknown = column.y,
         z_unknown = z.y,
         MS2_unknown = MS2.y) %>%
  select(compound_unknown, KRH_identification, compound_known, mz_unknown, mz_known, rt_seconds_unknown, rt_seconds_known, 
         column_unknown, column_known, z_unknown, z_known, MS2_unknown, MS2_known)  %>%
  arrange(compound_unknown)

## Confidence Level 1 ----------------------------------------  
# DHPS as example
Confidence.Level.1 <- My.Fuzzy.Join %>%
  filter(z_unknown == z_known,
         column_unknown == column_known) %>%
  mutate(mz_similarity_score = exp(-0.5 * (((mz_unknown - mz_known) / mz.flexibility) ^ 2)),
         rt_similarity_score = exp(-0.5 * (((rt_seconds_unknown - rt_seconds_known) / rt.flexibility) ^ 2))) %>%
  rowwise() %>%
  mutate(MS2_cosine_similarity = ifelse(is.na(MS2_unknown) | is.na(MS2_known), 
                                        NA, MS2CosineSimilarity(MakeScantable(MS2_unknown), MakeScantable(MS2_known)))) %>%
  mutate(total_similarity_score = ifelse(is.na(MS2_cosine_similarity), 
                                         ((mz_similarity_score + rt_similarity_score) / 2) * 100,
                                            ((MS2_cosine_similarity + mz_similarity_score + rt_similarity_score) / 3) * 100))

# Sanitycheck: add here
No.Fuzzy.Match <- Unknowns %>%
  select(compound_unknown) %>%
  filter(compound_unknown %in% setdiff(1:nrow(Unknowns), My.Fuzzy.Join$compound_unknown)) %>%
  pull()
CL1 <- sort(unique(Confidence.Level.1$compound_unknown))
everything.else <- setdiff(1:nrow(Unknowns), sort(c(CL1, No.Fuzzy.Match)))
all.unknowns <- sort(c(No.Fuzzy.Match, CL1, everything.else))

everything.else.df <- Unknowns %>%
  filter(compound_unknown %in% everything.else)
No.Fuzzy.Match.df <- Unknowns %>%
  filter(compound_unknown %in% No.Fuzzy.Match)

Mission.Accomplished <- My.Fuzzy.Join %>%
  left_join(Confidence.Level.1) %>%
  mutate(confidence_rank = ifelse(mz_similarity_score == 1 & rt_similarity_score == 1, 1, NA),
         confidence_source = ifelse(!is.na(confidence_rank), "Ingalls_Standards", NA)) %>%
  left_join(everything.else.df) %>%
  bind_rows(No.Fuzzy.Match.df)
