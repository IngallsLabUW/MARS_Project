library(fuzzyjoin)
library(tidyverse)


# Katherine csv: fully curated
MF_ClusterAssignments_Katherine <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") 

# Filter by Confidence Level 1----------------------------------------
Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                               stringsAsFactors = FALSE, header = TRUE) %>%
  select(Column, Compound.Name_old, Emperical.Formula:z, Fraction1, Fraction2) %>%
  filter(Compound.Name_old %in% MF_ClusterAssignments_Katherine$Identification)

Filtered_KRH_ClusterAssignments <- MF_ClusterAssignments_Katherine %>%
  select(!contains("cluster")) %>%
  filter(Confidence == 1)

## Check equality. All the same except that KRH ID'd iodine, while that isn't in the standards list.
IDd_KRH <- as.data.frame(sort(unique(Filtered_KRH_ClusterAssignments$Identification))) %>%
  rename(Compound = 1)
IDd_stds <- as.data.frame(sort(unique(Standards$Compound.Name_old))) %>%
  rename(Compound = 1)
Difference = setdiff(IDd_KRH, IDd_stds)

# Test run of ID'ing only Stds by mz ----------------------------------------
Unknowns <- MF_ClusterAssignments_Katherine %>%
  mutate(Unknown.Compound = paste("Compound_", 1:nrow(MF_ClusterAssignments_Katherine), sep = "")) %>%
  rename(m.z = mz) %>%
  select(Unknown.Compound, m.z) 

Knowns <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
           stringsAsFactors = FALSE, header = TRUE) %>%
  select(Compound.Name_old, m.z)

MyFuzzyJoin <- Unknowns %>%
  difference_inner_join(Knowns, by = c("m.z"), max_dist = 0.5)
  