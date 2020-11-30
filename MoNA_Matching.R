library(fuzzyjoin)
library(tidyverse)
source("Functions.R")
source("KRH_Functions.R")


# Specific Questions ------------------------------------------------------

# - The massbankMS2MatchNEGATIVE function doesn't like when the experimental data comes with a NA in the MS2 column. 
#   For now I am just filtering those rows out (line 67), but would it be preferable to hang on to that data somewhere?
## I think that's fine for now - for a 'level 2' confidence match you definately need an MS2, so for those mass featuers without MS2s, we can't get any 'level 2' confidence match.

# - Currently we are filtering out any matches that have a cosine similarity < 0.8 (line 74). 
#   Is it ok to totally throw those out or better to keep them as a lower-confidence match? 
## I think it'd be better to have the cosine similarity score cut off as an input and probably lower it to maybe 0.5 since we will then incorporate it into a 'total score'

#   Ditto on the MassBankppm < 5 filter (line 83).
# - There is a step in the Candidates df production that combines the experimental "MF_Frac" column 
#   with the narrowed-down candidates from MoNA (lines 36 and 76). In the original code this inserts a character vector 
#   of the whole MF_Frac column into the dataframe cell. Better to match back to the relevant MF_Frac identification?
## I think it'd be better to have the pmme cut off as an input and then incorporate it into a 'total score'


# TODO
# Need to filter out the standards from this "check". No point in assigning things we already have!

# Imports -----------------------------------------------------------------
Spectra <- read.csv("NEG_Spectra.csv") # One of four relational spreadsheets

ShortestDat <- read.csv("MFCluster_Assignments_Katherine.csv") %>% 
  ## ShortestDat is the unknowns, aka experimental values
  rename(MF_Frac = MassFeature_Column) %>%
  select(MF_Frac, mz, MS2) %>%
  as.data.frame()

# Originally from the massbankMS2MatchNEGATIVE function
# Edited by RML

mz <- as.numeric(unlist(ShortestDat["mz"])) # Had to add the unlist() here, otherwise threw a conversion error
MS2 <- as.character(ShortestDat["MS2"]) 
#MF_Frac <- as.data.frame(ShortestDat["MF_Frac"]) # Changed from character vector to dataframe for making the NoMatchReturn later.
# Works for line 62 but not 77.
MF_Frac <- as.character(ShortestDat["MF_Frac"]) # original. Works for line 77 but not for line 62.
Spectra.MH.Mass <- Spectra %>%
  mutate(MH_mass = M_mass - 1.0072766) 

ShortestDatForJoin <- ShortestDat %>% 
  rename(MH_mass = mz) %>%
  select(MH_mass, MS2)


# Create candidate dataframe ----------------------------------------------
testmass <- 116.07093 #Mass from the ShortestDat
Test <- Spectra.MH.Mass %>% filter(MH_mass > testmass-.020) %>% 
  filter(MH_mass < testmass+.020)  # All of these spectra should be in your 'candidates', not sure why they are not in the following lines


Candidates <- data.frame(near(Spectra.MH.Mass$MH_mass, ShortestDat$mz, tol = 0.02)) %>% # Create T/F matching dataframe
  # Join with Spectra df
  cbind(Spectra.MH.Mass) %>% 
  # Rename and filter for only those compounds that matched
  rename(mz.match = 1) %>%
  filter(mz.match == TRUE) %>% 
  # Join with ShortestDat to see which experimental mzs found near matches in the MoNA Spectra df.
  # Like the doing the near() filter in reverse.
  difference_inner_join(ShortestDatForJoin, by = "MH_mass", max_dist = 0.02) %>% 
  mutate(scan1 = spectrum_KRHform_filtered,
         scan2 = MS2, # from 
         mass1 = MH_mass.x, # this should be from Spectra
         mass2 = MH_mass.y) # this should be from experimental data

NoMatchReturn <- MF_Frac %>%
  mutate(MassBankMatch = NA,
         MassBankppm = NA,
         MassBankCosine1 = NA)

# Currently filtering out NAs until this is cleaned up
Candidates <- Candidates[c(1:5, 7), ]
  
if (length(Candidates$ID > 1)) { 
  Candidates$Cosine1 <- apply(Candidates, 1, FUN=function(x) MSMScosine1_df(x)) 
  Candidates$Cosine2 <- apply(Candidates, 1, FUN=function(x) MSMScosine2_df(x))
  
  Candidates2 <- Candidates %>%
    filter(Cosine1 > 0.8) %>% 
    arrange(desc(Cosine1)) %>% 
    mutate(MF_Frac = MF_Frac)
  
  Candidates3 <- Candidates2 %>% 
    mutate(MassBankMatch = paste(Names, ID, sep= " ID:"), 
           MassBankppm = abs(mass1-mass2)/mass2 *10^6,
           MassBankCosine1 = Cosine1) %>%
    #filter(MassBankppm < 5) %>% # In this data, this step removes everything.
    select(MF_Frac, MassBankMatch:MassBankCosine1)
}
if (length(Candidates$MF_Frac) == 0) {
  Candidates <- NoMatchReturn
}