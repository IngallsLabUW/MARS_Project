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
# eV vs V in collision energy?
# Check exact mass on mona, does it include the proton?

# Four relational spreadsheets from the MoNA download
MoNA.Spectra.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Spectra.csv")
MoNA.Names.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Names.csv")
MoNA.MetaData.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_MetaData.csv")
MoNA.CmpInfo.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_CmpInfo.csv")

relevant <- c("exact mass", "retention time", "precursor m/z", "collision energy")

MoNA.Known <- MoNA.MetaData.Neg %>%
  select(name, value, SpectraID) %>%
  filter(name %in% relevant) %>%
  mutate(name = str_replace(name, " ", "_")) %>%
  pivot_wider(., names_from = "name", values_from = "value") %>%
  filter(str_detect(retention_time, "min"),
         !str_detect(retention_time, "N/A")) %>%
  mutate(retention_time = gsub(" .*", "", retention_time)) %>%
  mutate(rt_seconds_known = as.numeric(retention_time) * 60) %>%
  select(-retention_time) %>%
  mutate(collision_energy = ifelse(!str_detect(collision_energy, "V"), NA, collision_energy))
  

# For initial MoNa matching, w/out ms2
Unknowns4Mona.neg <- read.csv("data_processed/confidence_level1.csv") %>%
  select(-X) %>%
  filter(z_unknown == -1) %>%
  mutate(MH_mass = mz_unknown - 1.0072766) 

KnownsfromMona.neg <- MoNA.Known %>%
  rename(MH_mass = exact_mass) %>%
  mutate(MH_mass = as.numeric(MH_mass))


MoNA.Fuzzy.Join <- KnownsfromMona.neg %>%
  difference_left_join(Unknowns4Mona.neg, by = c("MH_mass"), max_dist = 0.02) 


# Experimental spectra from lab, Confidence Level 1
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
  select(-Has.MS2) #%>%
  #select(-compound_unknown)

Experimental.Spectra.NoMS2 <- MS2.Test %>%
  filter(Has.MS2 == FALSE) %>%
  select(-Has.MS2)

IsolateMoNACandidates2 <- function(MoNA.Mass) {
  potential.candidates <- MoNA.Spectra.MHMass %>% 
    filter(MH_mass > MoNA.Mass - 0.020,  
           MH_mass < MoNA.Mass + 0.020) %>% 
    difference_inner_join(Experimental.Spectra.ForJoin, by = "MH_mass", max_dist = 0.02) %>% 
    mutate(scan1 = spectrum_KRHform_filtered, # scan1 is MS2 from MoNA
           scan2 = MS2_unknown,                       # scan2 is MS2 from the experimental data
           mass1 = MH_mass.x,                 # mass1 is the primary mass from MoNA
           mass2 = MH_mass.y)                 # mass2 is the primary mass from experimental data
  
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
