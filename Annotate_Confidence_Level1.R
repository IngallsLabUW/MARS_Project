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
Confidence.Level.1 <- My.Fuzzy.Join %>%
  filter(z_unknown == z_known,
         column_unknown == column_known) %>%
  mutate(mz_similarity_score = exp(-0.5 * (((mz_unknown - mz_known) / mz.flexibility) ^ 2)),
         rt_similarity_score = exp(-0.5 * (((rt_seconds_unknown - rt_seconds_known) / rt.flexibility) ^ 2))) %>%
##################
# test <- Confidence.Level.1 %>%
#   select(compound_unknown:rt_seconds_known, mz_similarity_score, rt_similarity_score)
##################
  rowwise() %>%
  mutate(MS2_cosine_similarity = ifelse(is.na(MS2_unknown) | is.na(MS2_known), 
                                        NA, MS2CosineSimilarity(MakeScantable(MS2_unknown), MakeScantable(MS2_known)))) %>%
  mutate(total_similarity_score = ifelse(is.na(MS2_cosine_similarity), 
                                         ((mz_similarity_score + rt_similarity_score) / 2) * 100,
                                            ((MS2_cosine_similarity + mz_similarity_score + rt_similarity_score) / 3) * 100))

mission.accomplished <- My.Fuzzy.Join %>%
  left_join(Confidence.Level.1) %>%
  mutate(confidence_rank = ifelse(mz_similarity_score == 1 & rt_similarity_score == 1, 1, NA),
         confidence_source = ifelse(!is.na(confidence.rank), "Ingalls_Standards", NA))



## Confidence Level A3 ----------------------------------------

within10duplicate <- function(df, column) {
  if (column > 8) # replace with "more than one unique match" rather than a number
    # Also "if" statement is probably unnecessary
    df2 <- df %>%
      mutate(RT.diff = abs(RT.seconds_Unknowns - RT.seconds_Standards)) %>%
      group_by(Unknown.Compound) %>%
      mutate(Closest.rt.Match = min(RT.diff)) %>%
      mutate(Closest.rt.Match2 = ifelse(Closest.rt.Match == RT.diff, TRUE, FALSE))
  else df
}

A3Confidence_MS2s <- My.Fuzzy.Join %>% # mz and 10 RT
  filter(!Unknown.Compound %in% unique(Confidence.Level.1$Unknown.Compound),
         !Unknown.Compound %in% unique(A1Confidence$Unknown.Compound),
         !Unknown.Compound %in% unique(A2Confidence_MS2s$Unknown.Compound),
         !Unknown.Compound %in% unique(A2Confidence$Unknown.Compound)) %>%
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  group_by(Unknown.Compound) %>%
  add_tally() %>%
  within10duplicate(column = "n") %>%
  ungroup() %>%
  mutate(mz_Similarity = exp(-0.5 * (((mz_Unknowns - mz_Standards) / 0.02) ^ 2)),
         RT_Similarity = exp(-0.5 * (((RT.seconds_Unknowns - RT.seconds_Standards) / 0.04) ^ 2))) %>% # Why all the 0s?
  filter_at(vars(MS2_Unknowns, MS2_Standards),all_vars(!is.na(.))) %>%
  rowwise() %>% 
  mutate(MS2cosinesim = MS2CosineSimilarity(MakeScantable(MS2_Unknowns), MakeScantable(MS2_Standards))) %>%
  mutate(Total.Similarity.Score = ((MS2cosinesim + mz_Similarity + RT_Similarity) / 3) * 100)
