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
MoNA.Names.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Names.csv") %>%
  select(SpectraID, name)
MoNA.MetaData.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_MetaData.csv")
MoNA.CmpInfo.Neg <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_CmpInfo.csv")

# MoNA.Spectra.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/POS_Spectra.csv")
# MoNA.Names.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/POS_Names.csv") %>%
#   select(SpectraID, name)
# MoNA.MetaData.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/POS_MetaData.csv")
# MoNA.CmpInfo.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/POS_CmpInfo.csv")
# 
# 
# # Experimental Values-------------------------------------------------
# # MH_mass is the compound mass minus the weight of a proton
# Experimental <- read.csv("data_processed/confidence_level1.csv") %>%
#   filter(z_experimental == -1) %>%
#   rename(MH_mass = mz_experimental) %>%
#   select(-contains("column"), -z_experimental, -z_theoretical)
# 
# # Theoretical (MoNA) -----------------------------------------------------------
# # exact_mass is what we should be matching to
# Theoretical.MoNA <- MoNA.MetaData.Neg %>%
#   filter(name == "exact mass") %>%
#   mutate(name = str_replace(name, " ", "_"),
#          name = str_replace(name, "/", "")) %>%
#   pivot_wider(., names_from = "name", values_from = "value") %>%
#   filter(str_detect(retention_time, "min"),
#          !str_detect(retention_time, "N/A")) %>%
#   mutate(retention_time = gsub(" .*", "", retention_time),
#          rt_seconds_MoNA = as.numeric(retention_time) * 60) %>%
#   left_join(MoNA.Names.Neg, by = "SpectraID") %>%
#   rename(name_MoNA = name) %>%
#   mutate(MH_mass = as.numeric(exact_mass) - 1.0072766) %>%
#   select(SpectraID, MH_mass, name_MoNA, rt_seconds_MoNA) %>%
#   unique()
# 
# # Fuzzy Join --------------------------------------------------------------
# # Please note: the below selections/renamings are subject to change. Adjust
# # namings and selections according to your needs; be aware that adding more
# # columns may cause join explosions.
# MoNA.Fuzzy.Join <- Theoretical.MoNA %>%
#   difference_left_join(Experimental, by = c("MH_mass"), max_dist = 0.02) %>%
#   rename(MH_mass_MoNA = MH_mass.x,
#          MH_mass_unknown = MH_mass.y,
#          compound_standards = compound_known,
#          mz_standards = mz_known,
#          rt_seconds_standards = rt_seconds_known, 
#          mz_similarity_score_stds = mz_similarity_score,
#          rt_similarity_score_stds = rt_similarity_score,
#          total_similarity_score_stds = total_similarity_score) %>%
#   select(compound_unknown, KRH_identification, compound_standards, SpectraID, name_MoNA,
#          MH_mass_unknown, MH_mass_MoNA, mz_standards,
#          rt_seconds_unknown, rt_seconds_MoNA, rt_seconds_standards,
#          mz_similarity_score_stds, rt_similarity_score_stds, 
#          total_similarity_score_stds, confidence_rank, confidence_source) %>%
#   filter(!is.na(compound_unknown)) %>%
#   unique()
# 
# Confidence.Level.2 <- MoNA.Fuzzy.Join %>%
#   rowwise() %>%
#   mutate(mz_similarity_score_MoNA = exp(-0.5 * (((MH_mass_unknown - MH_mass_MoNA) / mz.flexibility) ^ 2)),
#          rt_similarity_score_MoNA = exp(-0.5 * (((rt_seconds_unknown - rt_seconds_MoNA) / rt.flexibility) ^ 2)))

# Experimental Spectra ----------------------------------------------------
Experimental.Spectra <- read.csv("data_processed/confidence_level1.csv") %>%
  select(compound_experimental, mz_experimental, MS2_experimental) %>%
  unique()

# Reassign variables from Experimental.Spectra
mz <- as.numeric(unlist(Experimental.Spectra["mz_experimental"])) 
MS2 <- as.character(Experimental.Spectra["MS2_experimental"]) 
Mass.Feature <- as.data.frame(Experimental.Spectra["compound_experimental"]) 

Experimental.Spectra.ForJoin <- Experimental.Spectra %>%
  filter(!is.na(MS2_experimental)) %>%
  rename(MH_mass = mz_experimental)

# Subtract hydrogen for reference database
MoNA.Spectra.MHMass <- MoNA.Spectra.Neg %>%
  mutate(MH_mass = M_mass - 1.0072766) %>%
  select(ID, Names, spectrum_KRHform_filtered, MH_mass)

# Tidy blank values within the dataframe
MoNA.Spectra.MHMass[MoNA.Spectra.MHMass == ""] <- NA
MoNA.Spectra.MHMass <- MoNA.Spectra.MHMass %>%
  drop_na()

IsolateMoNACandidates <- function(MoNA.Mass, experimental.df) {
  potential.candidates <- MoNA.Spectra.MHMass %>% 
    filter(MH_mass > MoNA.Mass - 0.020,  
           MH_mass < MoNA.Mass + 0.020) %>% 
    difference_inner_join(experimental.df, by = "MH_mass", max_dist = 0.02) %>% 
    rename(scan1 = spectrum_KRHform_filtered, # scan1 is MS2 from MoNA
           scan2 = MS2_experimental,          # scan2 is MS2 from the experimental data
           mass1 = MH_mass.x,                 # mass1 is the primary mass from MoNA
           mass2 = MH_mass.y) #%>%             # mass2 is the primary mass from experimental data
    # rename(MH_mass_MoNA = MH_mass.x,
    #        MH_mass_experimental = MH_mass.y)
  
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

  potential.candidates$Cosine1 <- apply(potential.candidates, 1, FUN = function(x) MakeMS2CosineDataframe(x))

  filtered.cosine.candidates <- potential.candidates %>%
    filter(Cosine1 > cosine.score.cutoff) %>%
    arrange(desc(Cosine1))

  final.candidates <- filtered.cosine.candidates %>%
    mutate(massbank_match = paste(Names, ID, sep = " ID:"),
           massbank_ppm = abs(mass1 - mass2) / mass2 * 10^6,
           massbank_cosine1 = Cosine1) %>%
    unique() %>%
    filter(massbank_ppm < massbank.ppm.cutoff)

  return(final.candidates)
}
experimental.df <- Experimental.Spectra.ForJoin 


numCores <- detectCores()
numCores

system.time(
  outputdf_parallel <- mclapply(MoNA.Spectra.MHMass["MH_mass"], IsolateMoNACandidates, experimental.df = experimental.df,  mc.cores = numCores)
)

outputdf_parallel <- bind_rows(outputdf_parallel) %>%
  unique()
