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
  rename(Rt.seconds_Standards = rt.x,
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
  select(Unknown.Compound, KRH.Identification, Compound_Standards, mz_Unknowns, mz_Standards, RT.seconds_Unknowns, Rt.seconds_Standards, RT.seconds_MS2s, 
         Column_Unknowns, Column_Standards, z_Unknowns, z_Standards, MS2_Unknowns, MS2_Standards) 

## Confidence Level A1 ----------------------------------------
A1Confidence <- MyFuzzyJoin %>% # mz and 0.02 RT
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  group_by(Unknown.Compound) %>%
  filter(!Rt.seconds_Standards < (RT.seconds_Unknowns - 0.02) &
           !Rt.seconds_Standards > (RT.seconds_Unknowns + 0.02)) %>%
  ungroup() %>% ## 0.02 will change to "cutoff 1" or something
  select(-RT.seconds_MS2s, -MS2_Unknowns, -MS2_Standards) %>%
  unique()


A1MS2 <- MyFuzzyJoin %>%
  filter(Unknown.Compound %in% A1Confidence$Unknown.Compound) %>%
  select(Unknown.Compound, Compound_Standards, MS2_Unknowns, MS2_Standards) 

tosplit <- A1MS2 %>%
  group_by(Unknown.Compound) %>% 
  group_split()

scantable_Unknowns <- lapply(tosplit, function(df) {
  lapply(unique(df[, 3]), MakeScantable)
})

scantable_Standards <- lapply(tosplit, function(df){
  apply(df[4], 1, MakeScantable)
})

####
# Filter out NAs and put elsewhere since it doesn't matter anyway
test <- A1MS2 %>%
  ungroup() %>%
  slice(1:56) %>%
  rowwise() %>% 
  mutate(cosinesim= MS2CosineSimilarity(MakeScantable(MS2_Unknowns), MakeScantable(MS2_Standards)))







## Confidence Level 1 MS2 matching
Uracil <- MyFuzzyJoin %>%
  filter(Compound_Standards == "Uracil") %>%
  select(Compound_Standards, MS2_Unknowns, MS2_Standards)

Uracil.Standards <- Uracil %>%
  select(MS2_Standards) %>%
  slice(1) %>%
  MakeScantable()

Uracil.Experimental <- Uracil %>%
  select(MS2_Unknowns) %>%
  slice(1) %>%
  MakeScantable()

uracil.MS2cosine.sim <- MS2CosineSimilarity(scan1 = Uracil.Experimental, scan2 = Uracil.Standards)

## MS1 Similarity --------------------------------
uracil.krh <- Experimental.Values %>%
  filter(Identification == "Uracil") %>%
  select(Identification, mz, rt)

uracil.standard <- read.csv("data_extra/Ingalls_Lab_Standards_Feb1.csv") %>%
  mutate(rt = (RT..min. * 60),
         mz = as.numeric(m.z)) %>%
  rename(compound = Compound.Name_old) %>%
  filter(compound == "Uracil") %>%
  select(compound, mz, rt)

exp.value.mz = uracil.krh$mz
theor.value.mz = uracil.standard$mz

exp.value.rt = uracil.krh$rt
theor.value.rt = uracil.standard$rt

MS1.mz.similarity <- exp(-0.5 * (((exp.value.mz - theor.value.mz) / 0.02) ^ 2))
MS1.rt.similarity <- exp(-0.5 * (((exp.value.rt - theor.value.rt) / 0.4) ^ 2))

# Total Similarity Score
# TS = ((MS2 Similarity + MS1 Similarity) / 2) * 100
Uracil.Total.Similarity_AllVariables <- ((uracil.MS2cosine.sim + MS1.mz.similarity + MS1.rt.similarity) / 3) * 100
Uracil.Total.Similarity_NoMS2 <- ((MS1.mz.similarity + MS1.rt.similarity) / 2) * 100

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
  
