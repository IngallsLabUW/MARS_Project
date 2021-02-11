library(fuzzyjoin)
library(parallel)
library(tidyverse)
source("Functions.R")
source("KRH_Functions.R")

Cosine.Score.Cutoff <- 0.5 
MassBank.ppm.Cutoff <- 5

# Four relational spreadsheets from the MoNA download
MoNA.Spectra <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Spectra.csv")
MoNA.Names <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Names.csv")
MoNA.MetaData <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_MetaData.csv")
MoNA.CmpInfo <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_CmpInfo.csv")

# Experimental spectra from lab 
filter.variables <- read.csv("data_underway/Filtered_For_CL1.csv")
All.Unknowns <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv")

Experimental.Spectra <- All.Unknowns %>%
  mutate(Unknown.Compound = paste("Compound_", 1:nrow(All.Unknowns), sep = "")) %>%
  filter(Unknown.Compound %in% filter.variables$Unknown.Compound) %>%
  rename(MF.Fraction = MassFeature_Column) %>%
  select(Unknown.Compound, MF.Fraction, mz, MS2) %>%
  as.data.frame()

# Reassign variables from Experimental.Spectra
mz <- as.numeric(unlist(Experimental.Spectra["mz"])) 
MS2 <- as.character(Experimental.Spectra["MS2"]) 
MF.Fraction <- as.data.frame(Experimental.Spectra["MF.Fraction"]) 

MS2.Test <- Experimental.Spectra %>% 
  rename(MH_mass = mz) %>%
  mutate(Has.MS2 = ifelse(is.na(MS2), FALSE, TRUE))

Experimental.Spectra.ForJoin <- MS2.Test %>%
  filter(Has.MS2 == TRUE) %>% # Matching only on MS2 cosine similarity
  select(-Has.MS2) %>%
  ###
  select(-Unknown.Compound)
  ###

Experimental.Spectra.NoMS2 <- MS2.Test %>%
  filter(Has.MS2 == FALSE) %>%
  select(-Has.MS2)

# Subtract hydrogen for reference database
MoNA.Spectra.MHMass <- MoNA.Spectra %>%
  mutate(MH_mass = M_mass - 1.0072766) %>%
  select(ID, Names, spectrum_KRHform_filtered, MH_mass) 
MoNA.Spectra.MHMass[MoNA.Spectra.MHMass == ""] <- NA
MoNA.Spectra.MHMass <- MoNA.Spectra.MHMass %>%
  drop_na()

numCores <- detectCores()
numCores

system.time(
  outputdf_parallel <- mclapply(MoNA.Spectra.MHMass["MH_mass"], IsolateMoNACandidates, mc.cores = numCores)
)

outputdf_parallel <- bind_rows(outputdf_parallel) %>%
  unique()

