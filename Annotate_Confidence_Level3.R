library(fuzzyjoin)
library(tidyverse)

# [3]: Possible matches: “Putatively characterized compound classes (e.g. based upon characteristic physicochemical properties of a chemical class of compounds,
# or by spectral similarity to known compounds of a chemical class).” Level 3 is essentially an educated guess, a step above the complete unknown of level 4. 
# Currently this has been done manually. An option is MS1 matches from MoNA/Metlin/KEGG. PubMed and HMDB are in the same boat.

# MS1 match within 0.02 to MoNA/Metlin/KEGG, z match, and column match.

Experimental.Spectra <- read.csv("data_processed/confidence_level2.csv") %>%
  select(compound_experimental:KRH_identification, mz_experimental, z_experimental, confidence_rank, confidence_source) %>%
  unique() %>%
  group_by(compound_experimental) %>% 
  add_tally() %>%
  mutate(testcol = ifelse(n == 2 & is.na(confidence_rank), TRUE, FALSE)) %>%
  filter(testcol != TRUE) %>%
  rename(MH_mass = mz_experimental) %>%
  select(-n, -testcol)

MoNA.Spectra.NEG <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Spectra.csv") %>%
  mutate(z_MoNA = -1)
MoNA.Spectra.Pos <- read.csv("data_extra/MoNA_RelationalSpreadsheets/NEG_Spectra.csv") %>%
  mutate(z_MoNA = 1)

MoNA.Spectra <- MoNA.Spectra.NEG %>%
  rbind(MoNA.Spectra.Pos) %>%
  mutate(MH_mass = M_mass - 1.0072766) %>%
  select(ID, Names, MH_mass, z_MoNA) %>%
  mutate_all(., list(~na_if(.,""))) %>%
  drop_na() 

# Compare (remember some get filtered out)
CL3 <- MoNA.Spectra %>%
  difference_left_join(Experimental.Spectra, by = c("MH_mass"), max_dist = 0.02) %>%
  filter(is.na(confidence_source),
         z_MoNA == z_experimental) %>% # options filtered out here too
  rename(MH_mass_MoNA = MH_mass.x,
         MH_mass_experimental = MH_mass.y) %>%
  mutate(mz_similarity_score = exp(-0.5 * (((MH_mass_experimental - MH_mass_MoNA) / 0.02) ^ 2))) %>%
  arrange(compound_experimental)

# Sanity check -------------------------------------------------------------
No.CL3.Match <- setdiff(1:nrow(Experimental.Spectra),
                        sort(unique(CL3$compound_experimental)))
cl.match <- as.integer(unique(CL3$compound_experimental))

# Have any compounds been lost? Check for a TRUE output
all.experimentals <- sort(c(No.CL3.Match, cl.match))
length(all.experimentals) == length(Experimental.Spectra$compound_experimental)

# Make "no match" dataframes for comparison
No.CL3.Match.df <- Experimental.Spectra %>%
  filter(compound_experimental %in% No.CL3.Match)


final <- CL3 %>%
  bind_rows(No.CL3.Match.df) %>%
  mutate(confidence_rank3 = ifelse(mz_similarity_score > 0.9, 3, NA),
         confidence_source = ifelse(!is.na(confidence_rank), "MoNA", NA),
         confidence_rank = ifelse(is.na(confidence_rank) & !is.na(confidence_rank3), confidence_rank3, confidence_rank)) %>%
  select(compound_experimental, ID, Names, KRH_identification, MH_mass, MH_mass_MoNA, MH_mass_experimental,
         z_MoNA, z_experimental, mz_similarity_score, confidence_rank, confidence_source, mz_similarity_score)

