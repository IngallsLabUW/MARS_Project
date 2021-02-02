library(fuzzyjoin)
library(tidyverse)


# Experimental values
Experimental.Values <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") %>%
  select(-contains("cluster"))

# Lab MS2 values
Cyano.MS2 <- read.csv("data_extra/Cyano_stds_withMS2s.csv") %>%
  mutate(Column = "Cyano")
HILICPos.MS2 <- read.csv("data_extra/HILICPos_stds_withMS2s.csv") %>%
  mutate(Column = "HILICPos")
HILICNeg.MS2 <- read.csv("data_extra/HILICNeg_stds_withMS2s.csv") %>%
  mutate(Column = "HILICNeg")

Standards.MS2 <- Cyano.MS2 %>%
  rbind(HILICPos.MS2, HILICNeg.MS2) %>%
  mutate(rt = rt * 60) %>%
  rename(rt_ms2 = rt)
  
Standards <- read.csv("data_extra/Ingalls_Lab_Standards_Feb1.csv") %>%
  mutate(rt = (RT..min. * 60)) %>%
  select(Compound.Name_old, Column, m.z, rt, ionization_form, z)


## Confidence Level 1 Matching ----------------------------------------
Unknowns <- Experimental.Values %>% # Experimental
  mutate(Unknown.Compound = paste("Compound_", 1:nrow(Experimental.Values), sep = "")) %>%
  rename(m.z = mz,
         KRH.Identification = Identification) %>%
  select(Unknown.Compound, KRH.Identification, m.z, rt, Column, z, MS2) 

Knowns <- Standards %>% # our known standards
  rename(compound = Compound.Name_old,
         mz_known = m.z) %>%
  left_join(Standards.MS2, by = "compound")

  select(Standards.Compound.Name, m.z, RT.seconds, rt, Column, z.x, z.y, file_name, MS2s)


## Temporary test for mz_obs diffs
# There are no differences larger than 0.02 between mz and mz_obs.
# There are duplicate entries, from either a missed MS2 or a compound detected in multiple columns.
# Standards_mz <- Knowns %>% select(Standards.Compound.Name, m.z) %>% unique()
# MS2_mz <- Cyano.MS2 %>% # no differences larger than 0.02 between mz and mz_obs 
#   rbind(HILICPos.MS2, HILICNeg.MS2) %>% 
#   select(1:3, Column) %>% 
#   unique() %>% # There are duplicates! Sometimes a missed ms2, sometimes detected in both
#   rename(Standards.Compound.Name = compound) %>%
#   group_by(Standards.Compound.Name) %>%
#   add_tally()
# Difftest <- Standards_mz %>% # 
#   left_join(MS2_mz) %>%
#   mutate(diff = abs(m.z - mz)) %>%
#   filter(diff > 0.02)


MyFuzzyJoin <- Unknowns %>%
  difference_inner_join(Knowns, by = c("m.z"), max_dist = 0.02) %>% # needs to be swapped out for variable
  rename(m.z.Unknowns = m.z.x,
         Column.Unknowns = Column.x, # can switch these out for variables eventually?
         z.Unknowns = z.x,
         RT.seconds.Unknowns = rt.x,
         m.z.Standards = m.z.y,
         Column.Standards = Column.y,
         z.Standards = z.y,
         RT.seconds.Standards = RT.seconds,
         MS2.Unknowns = MS2) %>%
  select(Unknown.Compound, KRH.Identification, Standards.Compound.Name,
         m.z.Unknowns, m.z.Standards, RT.seconds.Unknowns, RT.seconds.Standards,
         Column.Unknowns, Column.Standards, z.Unknowns, z.Standards, MS2.Unknowns) 

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
  
