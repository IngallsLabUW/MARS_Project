library(fuzzyjoin)
library(parallel)
library(tidyverse)
source("Functions.R")

Cosine.Score.Cutoff <- 0.5 
MassBank.ppm.Cutoff <- 5

# Notes from this step
# Would like to extract the other online database keys and work that in to check matches there.
# Last Auto-Curation, looks useful but need to figure out what it is/what it means exactly. How does this relate to "date"?
# There is a LOT of metadata per compound in MoNA. Definitely open to incorporating more past what is already here! 
# Retention time issues on the metadata sheet. Will need custom attention for seconds, minutes, random words, etc...
# Varied rt from some other parameter. That will need to be ID'd for what's causing the range of values. 

# Scraped data from MoNA --------------------------------------------------
# Four relational spreadsheets
MoNA.Spectra.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Spectra.csv")
MoNA.Names.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Names.csv") %>%
  select(SpectraID, name)
MoNA.MetaData.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_MetaData.csv")
MoNA.CmpInfo.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_CmpInfo.csv")

relevant.parameters <- c("exact mass", "retention time")

# Unknowns (experimental) -------------------------------------------------
# MH_mass is the compound mass minus the weight of a proton
Unknowns <- read.csv("data_processed/confidence_level1.csv") %>%
  select(-X) %>%
  filter(z_unknown == -1) %>%
  rename(MH_mass = mz_unknown)

# Knowns (MoNA) -----------------------------------------------------------
# exact_mass is what we should be matching to
Knowns.MoNA <- MoNA.MetaData.Neg %>%
  filter(name %in% relevant.parameters) %>%
  mutate(name = str_replace(name, " ", "_"),
         name = str_replace(name, "/", "")) %>%
  pivot_wider(., names_from = "name", values_from = "value") %>%
  filter(str_detect(retention_time, "min"),
         !str_detect(retention_time, "N/A")) %>%
  mutate(retention_time = gsub(" .*", "", retention_time),
         rt_seconds_MoNA = as.numeric(retention_time) * 60) %>%
  left_join(MoNA.Names.Neg, by = "SpectraID") %>%
  rename(name_MoNA = name) %>%
  mutate(MH_mass = as.numeric(exact_mass) - 1.0072766) %>%
  select(SpectraID, MH_mass, name_MoNA, rt_seconds_MoNA) 

# Fuzzy Join --------------------------------------------------------------
MoNA.Fuzzy.Join <- Knowns.MoNA %>%
  difference_left_join(Unknowns, by = c("MH_mass"), max_dist = 0.02) %>%
  rename(MH_mass_MoNA = MH_mass.x,
         MH_mass_unknown = MH_mass.y,
         compound_standards = compound_known,
         mz_standards = mz_known,
         rt_seconds_standards = rt_seconds_known, 
         column_standards = column_known, 
         z_standards = z_known, 
         MS2_standards = MS2_known,
         mz_similarity_score_stds = mz_similarity_score,
         rt_similarity_score_stds = rt_similarity_score,
         MS2_cosine_similarity_stds = MS2_cosine_similarity,
         total_similarity_score_stds = total_similarity_score) %>%
  select(compound_unknown, KRH_identification, compound_standards, SpectraID, name_MoNA,
         MH_mass_unknown, MH_mass_MoNA, mz_standards,
         rt_seconds_unknown, rt_seconds_MoNA, rt_seconds_standards,
         #MS2_standards, MS2_unknown,
         mz_similarity_score_stds, rt_similarity_score_stds, MS2_cosine_similarity_stds,
         total_similarity_score_stds, confidence_rank, confidence_source) %>%
  # select(compound_unknown, KRH_identification, compound_standards, SpectraID, name_MoNA,
  #        MH_mass_unknown, MH_mass_MoNA, mz_standards,
  #        rt_seconds_unknown, rt_seconds_MoNA, rt_seconds_standards,
  #        column_unknown, column_standards,
  #        z_unknown, z_standards,
  #        MS2_standards, MS2_unknown,
  #        mz_similarity_score_stds, rt_similarity_score_stds, MS2_cosine_similarity_stds,
  #        total_similarity_score_stds, confidence_rank, confidence_source) %>%
  filter(!is.na(compound_unknown)) %>%
  unique()

Confidence.Level.2 <- MoNA.Fuzzy.Join %>%
  rowwise() %>%
  mutate(mz_similarity_score_MoNA = exp(-0.5 * (((MH_mass_unknown - MH_mass_MoNA) / mz.flexibility) ^ 2)),
         rt_similarity_score_MoNA = exp(-0.5 * (((rt_seconds_unknown - rt_seconds_MoNA) / rt.flexibility) ^ 2)))


# Experimental Spectra ----------------------------------------------------
Experimental.Spectra <- read.csv("data_processed/confidence_level1.csv") %>%
  select(compound_unknown, mz_unknown, MS2_unknown) %>%
  unique()

# Reassign variables from Experimental.Spectra
mz <- as.numeric(unlist(Experimental.Spectra["mz_unknown"])) 
MS2 <- as.character(Experimental.Spectra["MS2_unknown"]) 
MF.Fraction <- as.data.frame(Experimental.Spectra["compound_unknown"]) 

MS2.Test <- Experimental.Spectra %>% 
  rename(MH_mass = mz_unknown) %>%
  mutate(Has.MS2 = ifelse(is.na(MS2_unknown), FALSE, TRUE))

Experimental.Spectra.ForJoin <- MS2.Test %>%
  filter(Has.MS2 == TRUE) %>% # Matching only on MS2 cosine similarity
  select(-Has.MS2) 

Experimental.Spectra.NoMS2 <- MS2.Test %>%
  filter(Has.MS2 == FALSE) %>%
  select(-Has.MS2)

IsolateMoNACandidates2 <- function(MoNA.Mass) {
  potential.candidates <- MoNA.Spectra.MHMass %>% 
    filter(MH_mass > MoNA.Mass - 0.020,  
           MH_mass < MoNA.Mass + 0.020) %>% 
    difference_inner_join(Experimental.Spectra.ForJoin, by = "MH_mass", max_dist = 0.02) %>% 
    mutate(scan1 = spectrum_KRHform_filtered, # scan1 is MS2 from MoNA
           scan2 = MS2_unknown,               # scan2 is MS2 from the experimental data
           mass1 = MH_mass.x,                 # mass1 is the primary mass from MoNA
           mass2 = MH_mass.y) %>%                # mass2 is the primary mass from experimental data
    rename(MH_mass_MoNA = MH_mass.x,
           MH_mass_unknown = MH_mass.y)
  
  if (length(potential.candidates$ID) == 0) {
    print("There are no potential candidates.")
    No.Match.Return <- MF.Fraction %>%
      mutate(MassBankMatch = NA,
             MassBankppm = NA,
             MassBankCosine1 = NA)
    
    return(No.Match.Return)
  }
  
  # Add cosine similarity scores
  print("Making potential candidates")
  
  potential.candidates$Cosine1 <- apply(potential.candidates, 1, FUN = function(x) MakeMS2CosineDataframe(x)) 
  
  Candidates.Filtered.Cosine <- potential.candidates %>%
    filter(Cosine1 > Cosine.Score.Cutoff) %>%
    arrange(desc(Cosine1))
  
  final.candidates <- Candidates.Filtered.Cosine %>%
    mutate(MassBankMatch = paste(Names, ID, sep = " ID:"),
           MassBankppm = abs(mass1 - mass2) / mass2 * 10^6,
           MassBankCosine1 = Cosine1) %>%
    unique() %>%
    filter(MassBankppm < MassBank.ppm.Cutoff) 
  
  return(final.candidates)
}

# Subtract hydrogen for reference database
MoNA.Spectra.MHMass <- MoNA.Spectra.Neg %>%
  mutate(MH_mass = M_mass - 1.0072766) %>%
  select(ID, Names, spectrum_KRHform_filtered, MH_mass) 
MoNA.Spectra.MHMass[MoNA.Spectra.MHMass == ""] <- NA
MoNA.Spectra.MHMass <- MoNA.Spectra.MHMass %>%
  drop_na()

numCores <- detectCores()
numCores

system.time(
  outputdf_parallel <- mclapply(MoNA.Spectra.MHMass["MH_mass"], IsolateMoNACandidates2, mc.cores = numCores)
)

outputdf_parallel <- bind_rows(outputdf_parallel) %>%
  unique()

everything <- Confidence.Level.2 %>%
  left_join(outputdf_parallel, by = "compound_unknown") %>%
  unique()

# What got lost
No.MoNA.Match <- Unknowns %>%
  select(compound_unknown) %>%
  filter(compound_unknown %in% setdiff(1:nrow(Unknowns), MoNA.Fuzzy.Join$compound_unknown)) %>%
  pull() %>% 
  unique()