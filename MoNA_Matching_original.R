library(fuzzyjoin)
library(tidyverse)
source("Functions.R")
source("MyKRHFunctions.R")

# Need to filter out the standards from this "check". No point in assigning things we already have!

# Imports -----------------------------------------------------------------
Spectra <- read.csv("data_extra/NEG_Spectra.csv") # One of four relational spreadsheets

ShortestDat <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") %>% 
  ## ShortestDat is the unknowns, aka experimental values
  rename(MF_Frac = MassFeature_Column) %>%
  select(MF_Frac, mz, MS2) %>%
  as.data.frame()

# Originally from the massbankMS2MatchNEGATIVE function
# Edited by RML

mz <- as.numeric(unlist(ShortestDat["mz"])) # Had to add the unlist() here, otherwise threw a conversion error
MS2 <- as.character(ShortestDat["MS2"]) 
MF_Frac <- as.data.frame(ShortestDat["MF_Frac"]) # Changed from character vector to dataframe for making the NoMatchReturn later.
MF_Frac <- as.character(ShortestDat["MF_Frac"])
Spectra.MH.Mass <- Spectra %>%
  mutate(MH_mass = M_mass - 1.0072766) 

ShortestDatForJoin <- ShortestDat %>% 
  rename(MH_mass = mz) %>%
  select(MH_mass, MS2)


# Create candidate dataframe ----------------------------------------------
  
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
   # filter(MassBankppm < 5) %>% # what is this step? 
    select(MF_Frac, MassBankMatch:MassBankCosine1)
}
if (length(Candidates$MF_Frac) == 0) {
  Candidates <- NoMatchReturn
}


#####################################################
# ORIGINAL CANDIDATES + NOMATCHRETURN FRAME

#filter(near(MH_mass, mz, tol = 0.02)) %>%
# mutate(scan1 = spectrum_KRHform_filtered,
#        scan2 = MS2,
#        mass1 = MH_mass, # from Spectra
#        mass2 = mz) %>%
# filter(!is.na(scan1))

# NoMatchReturn <- Spectra %>%
#   mutate(MassBankMatch = NA,
#          MassBankppm = NA,
#          MassBankCosine1 = NA,
#          MF_Frac = MF_Frac) %>% # Doesn't work as character or as df right now
#   select(MF_Frac, MassBankMatch:MassBankCosine1) %>% head(1)
#####################################################

