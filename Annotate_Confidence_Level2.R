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
Experimental.Spectra <- read.csv("data_underway/Filtered_For_CL1.csv") 

# Reassign variables from Experimental.Spectra
mz <- as.numeric(unlist(Experimental.Spectra["mz_Unknowns"])) 
MS2 <- as.character(Experimental.Spectra["MS2_Unknowns"]) 
MF.Fraction <- as.data.frame(Experimental.Spectra["Unknown.Compound"]) 

MS2.Test <- Experimental.Spectra %>% 
  rename(MH_mass = mz_Unknowns) %>%
  mutate(Has.MS2 = ifelse(is.na(MS2_Unknowns), FALSE, TRUE))

Experimental.Spectra.ForJoin <- MS2.Test %>%
  filter(Has.MS2 == TRUE) %>% # Matching only on MS2 cosine similarity
  select(-Has.MS2)

Experimental.Spectra.NoMS2 <- MS2.Test %>%
  filter(Has.MS2 == FALSE) %>%
  select(-Has.MS2)

# Subtract hydrogen for reference database
MoNA.Spectra.MHMass <- MoNA.Spectra %>%
  mutate(MH_mass = M_mass - 1.0072766) %>%
  select(ID, Names, spectrum_KRHform_filtered, MH_mass) %>%
  drop_na()

numCores <- detectCores()
numCores

system.time(
  outputdf_parallel <- mclapply(MoNA.Spectra.MHMass["MH_mass"], IsolateMoNACandidates, mc.cores = numCores)
)

outputdf_parallel <- bind_rows(outputdf_parallel) %>%
  unique()

