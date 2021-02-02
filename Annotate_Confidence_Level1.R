library(fuzzyjoin)
library(tidyverse)


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
  select(compound, Column, rt, ionization_form, z)

## Confidence Level 1 Matching ----------------------------------------
Knowns <- Standards %>% # our known standards
  left_join(Standards.MS2, by = c("compound", "Column", "z")) %>%
  rename(RT.seconds_standards = rt.x,
         RT.seconds_MS2s = rt.y)

Unknowns <- Experimental.Values %>% # Experimental
  mutate(Unknown.Compound = paste("Compound_", 1:nrow(Experimental.Values), sep = "")) %>%
  rename(KRH.Identification = Identification) %>%
  select(Unknown.Compound, KRH.Identification, mz, rt, Column, z, MS2) 

MyFuzzyJoin <- Unknowns %>%
  difference_left_join(Knowns, by = c("mz"), max_dist = 0.02) %>% # needs to be swapped out for variable
  rename(mz_Unknowns = mz.x,
         Column_Unknowns = Column.x, # can switch these out for variables eventually?
         z_Unknowns = z.x,
         RT.seconds_Unknowns = rt,
         MS2_Unknowns = MS2.x,
         mz_Standards = mz.y,
         Column_Standards = Column.y,
         z_Standards = z.y,
         MS2_Standards = MS2.y) %>%
  select(Unknown.Compound, KRH.Identification, compound, mz_Unknowns, mz_Unknowns, RT.seconds_Unknowns, RT.seconds_standards, RT.seconds_MS2s, 
         Column_Unknowns, Column_Standards, z_Unknowns, z_Standards, MS2_Unknowns, MS2_Standards) 

A1Confidence <- MyFuzzyJoin %>% # mz and 0.02 RT
  filter(z.Unknowns == z.Standards,
         Column.Unknowns == Column.Standards) %>%
  group_by(Unknown.Compound) %>%
  filter(!RT.seconds.Standards < (RT.seconds.Unknowns - 0.02) &
           !RT.seconds.Standards > (RT.seconds.Unknowns + 0.02)) ## 0.02 will change to "cutoff 1" or something

## Confidence Level 1 MS2 matching
betaine_exp <- Experimental.Values %>%
  filter(Identification == "Betaine") %>%
  select(-contains("cluster"))

betaine_std <- Cyano.MS2 %>% 
  rbind(HILICPos.MS2, HILICNeg.MS2) %>%
  filter(compound == "Betaine",
         Column == "HILICPos")




## On to A2 confidence
within10duplicate <- function(df, column) {
  if (column > 1)
  df2 <- df %>% 
    slice(which.min(abs(RT.seconds.Unknowns - RT.seconds.Standards)))
  else df
}

A2Confidence <- MyFuzzyJoin %>% # mz and 10 RT
  filter(!Standards.Compound.Name %in% A1Confidence$Standards.Compound.Name,
         !KRH.Identification %in% A1Confidence$KRH.Identification) %>%
  filter(z.Unknowns == z.Standards,
         Column.Unknowns == Column.Standards) %>%
  filter(!RT.seconds.Standards < (RT.seconds.Unknowns - 10) & !RT.seconds.Standards > (RT.seconds.Unknowns + 10)) %>% ## 10 will change to "cutoff 2"
  filter(Unknown.Compound != "Compound_135") %>%
  group_by(Unknown.Compound) %>%
  add_tally() %>%
  within10duplicate(column = "n")
  

A3Confidence <- MyFuzzyJoin %>% # mz and closest RT
  filter(!Standards.Compound.Name %in% A1Confidence$Standards.Compound.Name,
         !KRH.Identification %in% A1Confidence$KRH.Identification) %>%
  filter(!Standards.Compound.Name %in% A2Confidence$Standards.Compound.Name,
         !KRH.Identification %in% A2Confidence$KRH.Identification) %>%
  filter(z.Unknowns == z.Standards,
         Column.Unknowns == Column.Standards) %>%
  group_by(KRH.Identification) %>%
  slice(which.min(abs(RT.seconds.Unknowns - RT.seconds.Standards)))

everything.else <- MyFuzzyJoin %>%
  filter(!Standards.Compound.Name %in% A1Confidence$Standards.Compound.Name,
         !KRH.Identification %in% A1Confidence$KRH.Identification) %>%
  filter(!Standards.Compound.Name %in% A2Confidence$Standards.Compound.Name,
         !KRH.Identification %in% A2Confidence$KRH.Identification) %>%
  filter(!Standards.Compound.Name %in% A3Confidence$Standards.Compound.Name,
         !KRH.Identification %in% A3Confidence$KRH.Identification)
  
