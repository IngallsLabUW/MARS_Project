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

Filtered_Standards <- Standards %>%
  select(Column, Compound.Name_old, Emperical.Formula:z, Fraction1, Fraction2) %>%
  filter(Compound.Name_old %in% MF_ClusterAssignments_Katherine$Identification) 

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

MyFuzzyJoin_mz <- Unknowns_mz %>%
  difference_inner_join(Knowns_mz, by = c("m.z"), max_dist = 0.05) %>%
  rename(Standards.m.z = m.z.y,
         Standards.Column = Column.y,
         Standards.z = z.y)

