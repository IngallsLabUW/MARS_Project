library(fuzzyjoin)
library(parallel)
library(tidyverse)
source("Functions.R")

cosine.score.cutoff <- 0.5 
massbank.ppm.cutoff <- 5

mz.flexibility <- 0.02
rt.flexibility <- 0.02 # seconds

# Scraped data from MoNA --------------------------------------------------
# Four relational spreadsheets
MoNA.Spectra.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Spectra.csv")

# MoNA.Names.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Names.csv") %>%
#   select(SpectraID, name)
# MoNA.MetaData.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_MetaData.csv")
# MoNA.CmpInfo.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_CmpInfo.csv")

MoNA.Spectra.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/POS_Spectra.csv")
# MoNA.Names.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/POS_Names.csv") %>%
#   select(SpectraID, name)
# MoNA.MetaData.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/POS_MetaData.csv")
# MoNA.CmpInfo.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/POS_CmpInfo.csv")


# Experimental Spectra ----------------------------------------------------
Experimental.Spectra <- read.csv("data_processed/confidence_level1.csv") %>%
  filter(!is.na(MS2_experimental),
         z_experimental == 1) %>%
  select(compound_experimental, mz_experimental, MS2_experimental) %>%
  rename(MH_mass = mz_experimental) %>%
  unique()

# Reassign variables from Experimental.Spectra
mz <- as.numeric(unlist(Experimental.Spectra["MH_mass"])) 
MS2 <- as.character(Experimental.Spectra["MS2_experimental"]) 
Mass.Feature <- as.data.frame(Experimental.Spectra["compound_experimental"]) 

# Subtract hydrogen for reference database
MoNA.Spectra.MHMass <- MoNA.Spectra.Pos %>%
  mutate(MH_mass = M_mass - 1.0072766) %>%
  select(ID, Names, spectrum_KRHform_filtered, MH_mass) %>%
  mutate_all(., list(~na_if(.,""))) %>%
  drop_na()

IsolateMoNACandidates <- function(MoNA.Mass, experimental.df) {
  potential.candidates <- MoNA.Spectra.MHMass %>% 
    filter(MH_mass > MoNA.Mass - 0.020,  
           MH_mass < MoNA.Mass + 0.020) %>% 
    difference_inner_join(experimental.df, by = "MH_mass", max_dist = 0.02) %>% 
    rename(scan1 = spectrum_KRHform_filtered, # scan1 is MS2 from MoNA
           scan2 = MS2_experimental,          # scan2 is MS2 from the experimental data
           mass1 = MH_mass.x,                 # mass1 is the mass from MoNA
           mass2 = MH_mass.y)                 # mass2 is the mass from experimental data
  
  if (length(potential.candidates$ID) == 0) {
    print("There are no potential candidates.")
    No.Match.Return <- Mass.Feature %>%
      mutate(massbank_match = NA,
             massbank_ppm = NA,
             massbank_cosine1 = NA)

    return(No.Match.Return)
  }

  # Add cosine similarity scores
  print("Making potential candidates")

  potential.candidates$massbank_cosine1 <- apply(potential.candidates, 1, FUN = function(x) MakeMS2CosineDataframe(x))

  final.candidates <- potential.candidates %>%
    mutate(massbank_match = paste(Names, ID, sep = " ID:"),
           massbank_ppm = abs(mass1 - mass2) / mass2 * 10^6) %>%
    unique() %>%
    filter(massbank_ppm < massbank.ppm.cutoff,
           massbank_cosine1 > cosine.score.cutoff) %>%
    arrange(desc(massbank_cosine1))

  return(final.candidates)
}
experimental.df <- Experimental.Spectra

numCores <- detectCores()
numCores

system.time(
  outputdf_parallel <- mclapply(MoNA.Spectra.MHMass["MH_mass"], IsolateMoNACandidates, experimental.df = experimental.df,  mc.cores = numCores)
)

outputdf_parallel <- bind_rows(outputdf_parallel) %>%
  unique()
