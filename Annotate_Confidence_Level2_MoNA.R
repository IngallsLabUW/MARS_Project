library(fuzzyjoin)
library(parallel)
library(tidyverse)
source("Functions.R")

cosine.score.cutoff <- 0.5 
massbank.ppm.cutoff <- 5

mz.flexibility <- 0.02
rt.flexibility <- 30 # seconds

# Are you comparing to positive or negative ion mode?
ion.mode <- "negative"

# Scraped data from MoNA --------------------------------------------------
# Four relational spreadsheets for each ion mode

# Subtract hydrogen for reference database
if (ion.mode == "negative") {
  MoNA.Spectra <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Spectra.csv") %>%
    mutate(MH_mass = M_mass - 1.0072766)
} else {
  MoNA.Spectra <- read.csv("data_extra/MoNA_RelationalSpreadsheets/POS_Spectra.csv") %>%
    mutate(MH_mass = M_mass + 1.0072766)
}

# Tidy theoretical spectra, dropping NA MS2s
MoNA.Spectra <- MoNA.Spectra %>%
  select(ID, Names, spectrum_KRHform_filtered, MH_mass) %>%
  mutate_all(., list(~na_if(.,""))) %>%
  drop_na()
  

MoNA.CmpInfo.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_CmpInfo.csv") %>%
  select(SpectraID, name, value) %>%
  filter(name == "molecular formula")
 
# MoNA.CmpInfo.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/POS_CmpInfo.csv") %>%
#   select(SpectraID, name, value) %>%
#   filter(name == "molecular formula")


# Experimental Spectra ----------------------------------------------------
Confidence.Level.1 <- read.csv("data_processed/confidence_level1.csv")
Experimental.Spectra <- Confidence.Level.1 %>%
  filter(!is.na(MS2_experimental),
         z_experimental == 1) %>%
  select(compound_experimental, mz_experimental, MS2_experimental) %>%
  rename(MH_mass = mz_experimental) %>%
  unique()

# Reassign variables from Experimental.Spectra
mz <- as.numeric(unlist(Experimental.Spectra["MH_mass"])) 
MS2 <- as.character(Experimental.Spectra["MS2_experimental"]) 
Mass.Feature <- as.data.frame(Experimental.Spectra["compound_experimental"]) 

# Assign variables for paralellized comparison
experimental.df <- Experimental.Spectra
numCores <- detectCores()
numCores

# Compare
MoNA.Matched <- mclapply(MoNA.Spectra["MH_mass"], IsolateMoNACandidates, experimental.df = experimental.df,  mc.cores = numCores) %>%
  bind_rows() %>%
  full_join(Confidence.Level.1) %>%
  select(compound_experimental, KRH_identification, compound_theoretical, massbank_match, ID, mz_experimental, mz_theoretical, mz_massbank, 
         rt_sec_experimental, rt_sec_theoretical, column_experimental, column_theoretical, z_experimental, z_theoretical, 
         MS2_experimental, MS2_theoretical, MS2_massbank, ppm_mass_error, massbank_ppm, mz_similarity_score, rt_similarity_score, 
         MS2_cosine_similarity, total_similarity_score, massbank_cosine_similarity, confidence_rank, confidence_source) %>%
  arrange(compound_experimental)

# Combine Confidence Level 2 with Confidence Level 1 ----------------------
Confidence.Level.2 <- MoNA.Matched %>%
  mutate(confidence_rank = ifelse(!is.na(massbank_match), 
                                  paste(confidence_rank, "2", sep = "; "), confidence_rank),
         confidence_source = ifelse(!is.na(massbank_match), 
                                    paste(confidence_source, "MoNA", sep = "; "), confidence_source)) %>%
  mutate(across(starts_with("confidence"), ~ReplaceNA(.x)))
