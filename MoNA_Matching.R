library(fuzzyjoin)
library(tidyverse)
source("Functions.R")
source("KRH_Functions.R")

# TODO
# Need to filter out the standards from this "check". No point in assigning things we already have!

# All code originally from the massbankMS2MatchNEGATIVE function. Edited by RML.
Cosine.Score.Cutoff <- 0.5 
MassBank.ppm.Cutoff <- 5

# Imports -----------------------------------------------------------------

# One of four relational spreadsheets from the MoNA download
MoNA.Spectra <- read.csv("data_extra/NEG_Spectra.csv") 

Experimental.Spectra <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") %>% 
  rename(MF.Fraction = MassFeature_Column) %>%
  select(MF.Fraction, mz, MS2) %>%
  as.data.frame()

# Reassign variables from Experimental.Spectra
mz <- as.numeric(unlist(Experimental.Spectra["mz"])) 
MS2 <- as.character(Experimental.Spectra["MS2"]) 
MF.Fraction <- as.data.frame(Experimental.Spectra["MF.Fraction"]) 

Experimental.Spectra.ForJoin <- Experimental.Spectra %>% 
  rename(MH_mass = mz) 

# Subtract hydrogen
MoNA.Spectra.MHMass <- MoNA.Spectra %>%
  mutate(MH_mass = M_mass - 1.0072766) 

# Create candidate and no matches dataframe ------------------------------------------------------

testmass <- 116.07093 # Test mass from the MoNA database

Candidates <- MoNA.Spectra.MHMass %>% 
  filter(MH_mass > testmass - 0.020,  
         MH_mass < testmass + 0.020) %>% 
  # Join with Experimental.Spectra to see which experimental mzs found near matches in the MoNA Spectra df.
  difference_inner_join(Experimental.Spectra.ForJoin, by = "MH_mass", max_dist = 0.02) %>%  
  mutate(scan1 = spectrum_KRHform_filtered, # scan1 is MS2 from MoNA 
         scan2 = MS2, # scan2 is MS2 from the experimental data 
         mass1 = MH_mass.x, # mass1 is the primary mass from MoNA 
         mass2 = MH_mass.y) # mass2 is the primary mass from experimental data

No.Match.Return <- MF.Fraction %>%
  mutate(MassBankMatch = NA, 
         MassBankppm = NA,
         MassBankCosine1 = NA)


# Add cosine similarity scores --------------------------------------------

if (length(Candidates$ID > 1)) { 
  Candidates$Cosine1 <- apply(Candidates, 1, FUN=function(x) MSMScosine1_df(x)) 
  
  Candidates.Filtered.Cosine <- Candidates %>%
    filter(Cosine1 > Cosine.Score.Cutoff) %>% 
    arrange(desc(Cosine1)) 
  
  Final.Candidates <- Candidates.Filtered.Cosine %>% 
    mutate(MassBankMatch = paste(Names, ID, sep = " ID:"), 
           MassBankppm = abs(mass1 - mass2) / mass2 * 10^6,
           MassBankCosine1 = Cosine1) %>%
    #filter(MassBankppm < MassBank.ppm.Cutoff) %>% # In the current data, this step removes everything.
    select(MF.Fraction, MassBankMatch:MassBankCosine1)

} else if (length(Candidates$MF.Fraction) == 0) {
  
    Final.Candidates <- No.Match.Return

}

# Test run for adply function ---------------------------------------------
# Here’s how I normally implemented the massbankMS2MatchPOSITIVE function I gave you earlier:
# Make sure whatever you try can be run easily with parallel processing, I think that was why I went with it.
# Where PosDat is a dataframe with mz and MS2 (in KRH form)

# plyr performs function on row of dataframe. has to be parallelizable. 

# parallel is differnet on mac/pcs. 

PositiveMBMatches <- adply(PosDat, 1, .fun=function(x) massbankMS2MatchPOSITIVE(x), .parallel = TRUE) %>%
  filter(!is.na(MassBankMatch))



