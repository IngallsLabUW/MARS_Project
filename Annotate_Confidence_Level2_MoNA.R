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

IsolateMoNACandidates <- function(MoNA.Mass, experimental.df) {
  potential.candidates <- MoNA.Spectra %>% 
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
             massbank_cosine_similarity = NA)

    return(No.Match.Return)
  }

  # Add cosine similarity scores
  print("Making potential candidates")

  potential.candidates$massbank_cosine_similarity <- apply(potential.candidates, 1, FUN = function(x) MakeMS2CosineDataframe(x))

  final.candidates <- potential.candidates %>%
    mutate(massbank_match = paste(Names, ID, sep = " ID:"),
           massbank_ppm = abs(mass2 - mass1) / mass1 * 10^6) %>% 
    rename(MS2_massbank = scan1,
           mz_massbank = mass1, # should maybe name it something else because of mh?...
           MS2_experimental = scan2,
           mz_experimental = mass2) %>%
    unique() %>%
    filter(massbank_ppm < massbank.ppm.cutoff,
           massbank_cosine_similarity > cosine.score.cutoff) %>%
    arrange(desc(massbank_cosine_similarity))

  return(final.candidates)
}

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
