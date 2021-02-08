library(fuzzyjoin)
library(tidyverse)
options(scipen = 999)
options(digits = 6)

source("Functions.R")

# TODO: see where KRH and MARS identification differ

# Experimental values
Experimental.Values <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") %>%
  select(-contains("cluster"))

# Lab MS2 values
Cyano.MS2 <- read.csv("data_extra/Cyano_stds_withMS2s.csv") %>%
  mutate(Column = "RP",
         z = NA)
HILICPos.MS2 <- read.csv("data_extra/HILICPos_stds_withMS2s.csv") %>%
  mutate(Column = "HILIC",
         z = 1)
HILICNeg.MS2 <- read.csv("data_extra/HILICNeg_stds_withMS2s.csv") %>%
  mutate(Column = "HILIC",
         z = -1)

Standards.MS2 <- Cyano.MS2 %>%
  rbind(HILICPos.MS2, HILICNeg.MS2) %>%
  mutate(rt = rt * 60) %>%
  rename(MS2 = MS2s) %>%
  select(-rtmin, -rtmax, -mz_obs)

Standards <- read.csv("data_extra/Ingalls_Lab_Standards_Feb1.csv") %>%
  mutate(rt = (RT..min. * 60),
         m.z = as.numeric(m.z)) %>%
  rename(mz = m.z,
         compound = Compound.Name_old) %>%
  select(compound, Column, rt, ionization_form, z) # no appreciable difference between mzs, easier to have just one

## Confidence Level 1 Matching ----------------------------------------
Knowns <- Standards %>% # our known standards
  left_join(Standards.MS2, by = c("compound", "Column", "z")) %>%
  rename(RT.seconds_Standards = rt.x,
         RT.seconds_MS2s = rt.y)

Unknowns <- Experimental.Values %>% # Experimental
  mutate(Unknown.Compound = paste("Compound_", 1:nrow(Experimental.Values), sep = "")) %>%
  rename(KRH.Identification = Identification) %>%
  select(Unknown.Compound, KRH.Identification, mz, rt, Column, z, MS2) 

MyFuzzyJoin <- Unknowns %>%
  difference_left_join(Knowns, by = c("mz"), max_dist = 0.02) %>% # needs to be swapped out for variable
  rename(Compound_Standards = compound,
         mz_Unknowns = mz.x,
         Column_Unknowns = Column.x, # can switch these out for variables eventually?
         z_Unknowns = z.x,
         RT.seconds_Unknowns = rt,
         MS2_Unknowns = MS2.x,
         mz_Standards = mz.y,
         Column_Standards = Column.y,
         z_Standards = z.y,
         MS2_Standards = MS2.y) %>%
  select(Unknown.Compound, KRH.Identification, Compound_Standards, mz_Unknowns, mz_Standards, RT.seconds_Unknowns, RT.seconds_Standards, RT.seconds_MS2s, 
         Column_Unknowns, Column_Standards, z_Unknowns, z_Standards, MS2_Unknowns, MS2_Standards) 

## Confidence Level A1 ----------------------------------------
A1Confidence_MS2s <- MyFuzzyJoin %>% # mz and 0.02 RT, 0.02 will change to "cutoff 1" or something
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  group_by(Unknown.Compound) %>%
  filter(!RT.seconds_Standards < (RT.seconds_Unknowns - 0.02) &
           !RT.seconds_Standards > (RT.seconds_Unknowns + 0.02)) %>% 
  ungroup() %>%
  mutate(mz_Similarity = exp(-0.5 * (((mz_Unknowns - mz_Standards) / 0.02) ^ 2)),
         RT_Similarity = exp(-0.5 * (((RT.seconds_Unknowns - RT.seconds_Standards) / 0.04) ^ 2))) %>% # Check cutoff and similarity of 1 always?
  filter_at(vars(MS2_Unknowns, MS2_Standards), all_vars(!is.na(.))) %>% 
  rowwise() %>% 
  mutate(MS2cosinesim = MS2CosineSimilarity(MakeScantable(MS2_Unknowns), MakeScantable(MS2_Standards))) %>%
  mutate(Total.Similarity.Score = ((MS2cosinesim + mz_Similarity + RT_Similarity) / 3) * 100)

A1Confidence <- MyFuzzyJoin %>% 
  filter(!Unknown.Compound %in% A1Confidence_MS2s$Unknown.Compound) %>%     
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  group_by(Unknown.Compound) %>%
  filter(!RT.seconds_Standards < (RT.seconds_Unknowns - 0.02) &
           !RT.seconds_Standards > (RT.seconds_Unknowns + 0.02)) %>% 
  ungroup() %>%
  mutate(mz_Similarity = exp(-0.5 * (((mz_Unknowns - mz_Standards) / 0.02) ^ 2)),
         RT_Similarity = exp(-0.5 * (((RT.seconds_Unknowns - RT.seconds_Standards) / 0.04) ^ 2))) %>% 
  rowwise() %>% 
  mutate(Total.Similarity.Score = ((mz_Similarity + RT_Similarity) / 2) * 100) ### See if this works


## Confidence Level A2 ---------------------------------------- # Having issues with a bajillion zeroes in RT similarity
A2Confidence_MS2s <- MyFuzzyJoin %>% # mz and 10 RT
  filter(!Unknown.Compound %in% unique(A1Confidence_MS2s$Unknown.Compound),
         !Unknown.Compound %in% unique(A1Confidence$Unknown.Compound)) %>%
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  filter(!RT.seconds_Standards < (RT.seconds_Unknowns - 10) &
           !RT.seconds_Standards > (RT.seconds_Unknowns + 10)) %>%
  ungroup() %>%
  mutate(mz_Similarity = exp(-0.5 * (((mz_Unknowns - mz_Standards) / 0.02) ^ 2)),
         RT_Similarity = exp(-0.5 * (((RT.seconds_Unknowns - RT.seconds_Standards) / 0.04) ^ 2))) %>% 
  filter_at(vars(MS2_Unknowns, MS2_Standards),all_vars(!is.na(.))) %>%
  rowwise() %>% 
  mutate(MS2cosinesim = MS2CosineSimilarity(MakeScantable(MS2_Unknowns), MakeScantable(MS2_Standards))) %>%
  mutate(Total.Similarity.Score = ((MS2cosinesim + mz_Similarity + RT_Similarity) / 3) * 100)

A2Confidence <- MyFuzzyJoin %>% 
  filter(!Unknown.Compound %in% unique(A1Confidence_MS2s$Unknown.Compound),
         !Unknown.Compound %in% unique(A1Confidence$Unknown.Compound),
         !Unknown.Compound %in% unique(A2Confidence_MS2s$Unknown.Compound)) %>%
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  filter(!RT.seconds_Standards < (RT.seconds_Unknowns - 10) &
           !RT.seconds_Standards > (RT.seconds_Unknowns + 10)) %>%
  ungroup() %>%
  mutate(mz_Similarity = exp(-0.5 * (((mz_Unknowns - mz_Standards) / 0.02) ^ 2)),
         RT_Similarity = exp(-0.5 * (((RT.seconds_Unknowns - RT.seconds_Standards) / 0.04) ^ 2))) %>% 
  rowwise() %>% 
  mutate(Total.Similarity.Score = ((mz_Similarity + RT_Similarity) / 2) * 100)

## Confidence Level A3 ----------------------------------------

within10duplicate <- function(df, column) {
  if (column > 8) # replace with "more than one unique match" rather than a number
    # Also "if" statement is probably unnecessary
    df2 <- df %>%
      mutate(RT.diff = abs(RT.seconds_Unknowns - RT.seconds_Standards)) %>%
      group_by(Unknown.Compound) %>%
      mutate(Closest.Match = min(RT.diff)) %>%
      mutate(Closest.Match2 = ifelse(Closest.Match == RT.diff, TRUE, FALSE))
  else df
}

A3Confidence_MS2s <- MyFuzzyJoin %>% # mz and 10 RT
  filter(!Unknown.Compound %in% unique(A1Confidence_MS2s$Unknown.Compound),
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

A3Confidence <- MyFuzzyJoin %>% # mz and 10 RT
  filter(!Unknown.Compound %in% unique(A1Confidence_MS2s$Unknown.Compound),
         !Unknown.Compound %in% unique(A1Confidence$Unknown.Compound),
         !Unknown.Compound %in% unique(A2Confidence_MS2s$Unknown.Compound),
         !Unknown.Compound %in% unique(A2Confidence$Unknown.Compound),
         !Unknown.Compound %in% unique(A3Confidence_MS2s$Unknown.Compound)) %>%
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  group_by(Unknown.Compound) %>%
  add_tally() %>%
  within10duplicate(column = "n") %>%
  ungroup() %>%
  mutate(mz_Similarity = exp(-0.5 * (((mz_Unknowns - mz_Standards) / 0.02) ^ 2)),
         RT_Similarity = exp(-0.5 * (((RT.seconds_Unknowns - RT.seconds_Standards) / 0.04) ^ 2))) %>% # Why all the 0s?
  rowwise() %>% 
  mutate(Total.Similarity.Score = ((mz_Similarity + RT_Similarity) / 2) * 100)


everything.else <- MyFuzzyJoin %>%
  filter(!Unknown.Compound %in% unique(A1Confidence_MS2s$Unknown.Compound),
         !Unknown.Compound %in% unique(A1Confidence$Unknown.Compound),
         !Unknown.Compound %in% unique(A2Confidence_MS2s$Unknown.Compound),
         !Unknown.Compound %in% unique(A2Confidence$Unknown.Compound),
         !Unknown.Compound %in% unique(A3Confidence_MS2s$Unknown.Compound),
         !Unknown.Compound %in% unique(A3Confidence$Unknown.Compound))

# At this point, add a column to say where it was identified and to what level of satisfaction.
