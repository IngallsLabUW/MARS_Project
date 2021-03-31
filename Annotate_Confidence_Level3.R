library(fuzzyjoin)
library(tidyverse)

# MS1 match within 0.02 to MoNA/Metlin/KEGG, z match, and column match.

# Experimental ------------------------------------------------------------
Confidence.Level.2 <- read.csv("data_processed/confidence_level2.csv")

Experimental.Values <- Confidence.Level.2 %>%
  select(compound_experimental:KRH_identification, mz_experimental, z_experimental, confidence_rank, confidence_source) %>%
  unique() %>%
  group_by(compound_experimental) %>%
  add_tally() %>%
  mutate(temp = ifelse(n == 2 & is.na(confidence_rank), TRUE, FALSE)) %>%
  filter(temp != TRUE) %>%
  rename(MH_mass = mz_experimental) %>%
  select(-n, -temp)


# MoNA  -------------------------------------------------------------------
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

# Compare MS1 -------------------------------------------------------------
My.Fuzzy.Join <- MoNA.Spectra %>%
  difference_left_join(Experimental.Values, by = c("MH_mass"), max_dist = 0.02) %>%
  filter(is.na(confidence_source),
         z_MoNA == z_experimental) %>% 
  rename(MH_mass_MoNA = MH_mass.x,
         MH_mass_experimental = MH_mass.y) %>%
  mutate(mz_similarity_score = exp(-0.5 * (((MH_mass_experimental - MH_mass_MoNA) / 0.02) ^ 2))) %>%
  arrange(compound_experimental)


# Sanity check -------------------------------------------------------------
No.Fuzzy.Match <- setdiff(1:nrow(Experimental.Values),
                        sort(unique(My.Fuzzy.Join$compound_experimental)))
CL3.Match <- as.integer(unique(My.Fuzzy.Join$compound_experimental))

# Have any compounds been lost? Check for a TRUE output
all.experimentals <- sort(c(No.Fuzzy.Match, CL3.Match))
length(all.experimentals) == length(Experimental.Values$compound_experimental)

# Make "no match" dataframes for comparison
No.Fuzzy.Match.df <- Experimental.Values %>%
  filter(compound_experimental %in% No.Fuzzy.Match) %>%
  rename(MH_mass_experimental = MH_mass)


# Go to MARS ------------------------------------------------
final <- My.Fuzzy.Join %>%
  bind_rows(No.Fuzzy.Match.df) %>%
  mutate(confidence_rank3 = ifelse(mz_similarity_score > 0.9, 3, NA),
         confidence_source = ifelse(!is.na(confidence_rank), "MoNA", NA),
         confidence_rank = ifelse(is.na(confidence_rank) & !is.na(confidence_rank3), confidence_rank3, confidence_rank)) %>%
  select(compound_experimental, "MoNA_ID" = ID, "MoNA_Names" = Names, KRH_identification, MH_mass_MoNA, MH_mass_experimental,
         z_MoNA, z_experimental, "mz_similarity_score_MoNa" = mz_similarity_score, confidence_rank, confidence_source) 

##ISSUES HAPPENING WITH CONFIDENCE SCORE
Confidence.Level.3 <- final %>%
  left_join(Confidence.Level.2, by = c("compound_experimental", "KRH_identification", "z_experimental",
                                       "confidence_rank", "confidence_source")) %>%
  arrange(compound_experimental)
