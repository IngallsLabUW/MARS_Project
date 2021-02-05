library(fuzzyjoin)
library(tidyverse)

source("Functions.R")


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
A1Confidence <- MyFuzzyJoin %>% # mz and 0.02 RT
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  group_by(Unknown.Compound) %>%
  filter(!RT.seconds_Standards < (RT.seconds_Unknowns - 0.02) &
           !RT.seconds_Standards > (RT.seconds_Unknowns + 0.02)) %>% ## 0.02 will change to "cutoff 1" or something
  ungroup() %>%
  mutate(mz_Similarity = exp(-0.5 * ((mz_Unknowns - mz_Standards) / 0.02)),
         RT_Similarity = exp(-0.5 * ((RT.seconds_Unknowns - RT.seconds_Standards) / 0.04))) %>% # Check cutoff and similarity of 1 always?
  filter_at(vars(MS2_Unknowns, MS2_Standards),all_vars(!is.na(.))) %>%
  rowwise() %>% 
  mutate(MS2cosinesim = MS2CosineSimilarity(MakeScantable(MS2_Unknowns), MakeScantable(MS2_Standards))) %>%
  mutate(Total.Similarity.Score = ((MS2cosinesim + mz_Similarity + RT_Similarity) / 3) * 100)


## On to A2 confidence
within10duplicate <- function(df, column) {
  if (column > 1)
  df2 <- df %>% 
    slice(which.min(abs(RT.seconds_Unknowns - Rt.seconds_Standards)))
  else df
}

A2Confidence <- MyFuzzyJoin %>% # mz and 10 RT
  filter(!Compound_Standards %in% A1Confidence$Compound_Standards,
         !KRH.Identification %in% A1Confidence$KRH.Identification) %>%
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  filter(!Rt.seconds_Standards < (RT.seconds_Unknowns - 10) & !Rt.seconds_Standards > (RT.seconds_Unknowns + 10)) %>% ## 10 will change to "cutoff #2"
  group_by(Unknown.Compound) %>%
  add_tally() %>%
  within10duplicate(column = "n")
  

A3Confidence <- MyFuzzyJoin %>% # mz and closest RT
  filter(!Compound_Standards %in% A1Confidence$Compound_Standards,
         !KRH.Identification %in% A1Confidence$KRH.Identification) %>%
  filter(!Compound_Standards %in% A2Confidence$Compound_Standards,
         !KRH.Identification %in% A2Confidence$KRH.Identification) %>%
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  group_by(KRH.Identification) %>%
  slice(which.min(abs(RT.seconds_Unknowns - Rt.seconds_Standards)))

everything.else <- MyFuzzyJoin %>%
  filter(!Standards.Compound.Name %in% A1Confidence$Standards.Compound.Name,
         !KRH.Identification %in% A1Confidence$KRH.Identification) %>%
  filter(!Standards.Compound.Name %in% A2Confidence$Standards.Compound.Name,
         !KRH.Identification %in% A2Confidence$KRH.Identification) %>%
  filter(!Standards.Compound.Name %in% A3Confidence$Standards.Compound.Name,
         !KRH.Identification %in% A3Confidence$KRH.Identification)
  
