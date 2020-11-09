library(fuzzyjoin)
library(tidyverse)
options(digits = 6)


## Katherine csv: fully curated
MF_ClusterAssignments_Katherine <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") 

## Filter by Confidence Level 1 ----------------------------------------
Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                               stringsAsFactors = FALSE, header = TRUE) %>%
  mutate(RT.seconds = RT..min. * 60) %>%
  select(-RT..min.)

Filtered_KRH_ClusterAssignments <- MF_ClusterAssignments_Katherine %>%
  select(!contains("cluster")) %>%
  filter(Confidence == 1)

## Check matching equality All the same except that KRH ID'd iodine, while that isn't in the standards list.
IDd_KRH <- as.data.frame(sort(unique(Filtered_KRH_ClusterAssignments$Identification))) %>%
  rename(Compound = 1)
IDd_stds <- as.data.frame(sort(unique(Standards$Compound.Name_old))) %>%
  rename(Compound = 1)
Difference = setdiff(IDd_KRH, IDd_stds)

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
         Column.Unknowns, Column.Standards, z.Unknowns, z.Standards) %>%
  filter(z.Unknowns == z.Standards,
         Column.Unknowns == Column.Standards)
###

RT_Testing <- MyFuzzyJoin %>%
  select(Unknown.Compound:Standards.Compound.Name, RT.seconds.Unknowns, RT.seconds.Standards) %>%
  group_by(Unknown.Compound) %>%
  add_tally() %>%
  mutate(testing = ifelse((RT.seconds.Unknowns > (RT.seconds.Standards + 0.2) | 
                             RT.seconds.Unknowns < (RT.seconds.Standards - 0.2)), 
                          FALSE, TRUE)) 
