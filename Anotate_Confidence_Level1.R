library(fuzzyjoin)
library(tidyverse)
options(digits = 6)


## Katherine csv: fully curated
MF_ClusterAssignments_Katherine <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") 

Standards <- read.csv("data_extra/Ingalls_Standards_Jan11.csv", stringsAsFactors = FALSE) %>%
  mutate(m.z = as.numeric(m.z)) %>%
  mutate(RT.seconds = RT..min. * 60)

# Filtered_KRH_ClusterAssignments <- MF_ClusterAssignments_Katherine %>%
#   select(!contains("cluster")) %>%
#   filter(Confidence == 1)

## Confidence Level 1 Matching ----------------------------------------
Unknowns_mz <- MF_ClusterAssignments_Katherine %>%
  mutate(Unknown.Compound = paste("Compound_", 1:nrow(MF_ClusterAssignments_Katherine), sep = "")) %>%
  rename(m.z = mz,
         KRH.Identification = Identification) %>%
  select(Unknown.Compound, KRH.Identification, m.z, rt, Column, z) 

Knowns_mz <- Standards %>%
  select(Compound.Name_old, m.z, RT.seconds, Column, z) %>%
  rename(Standards.Compound.Name = Compound.Name_old) 


MyFuzzyJoin <- Unknowns_mz %>%
  difference_inner_join(Knowns_mz, by = c("m.z"), max_dist = 0.05) %>%
  rename(m.z.Unknowns = m.z.x,
         Column.Unknowns = Column.x, # can switch these out for variables eventually
         z.Unknowns = z.x,
         RT.seconds.Unknowns = rt,
         m.z.Standards = m.z.y,
         Column.Standards = Column.y,
         z.Standards = z.y,
         RT.seconds.Standards = RT.seconds) %>%
  select(Unknown.Compound, KRH.Identification, Standards.Compound.Name,
         m.z.Unknowns, m.z.Standards, RT.seconds.Unknowns, RT.seconds.Standards,
         Column.Unknowns, Column.Standards, z.Unknowns, z.Standards) 

A1Confidence <- MyFuzzyJoin %>% # mz and 0.02 RT
  filter(z.Unknowns == z.Standards,
         Column.Unknowns == Column.Standards) %>%
  group_by(Unknown.Compound) %>%
  filter(!RT.seconds.Standards < (RT.seconds.Unknowns - 0.02) & !RT.seconds.Standards > (RT.seconds.Unknowns + 0.02)) ## 0.02 will change to "cutoff 1" or something


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
  
