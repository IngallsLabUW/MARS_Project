library(fuzzyjoin)
library(tidyverse)
options(scipen = 999)
options(digits = 6)

source("Functions.R")

mz.flexibility <- 0.02
rt.flexibility <- 30 # seconds

# Theoretical Values ----------------------------------------------------------------------
# Gather and summarize theoretical values from the standards sheet and the Ingalls MS2 data
# MS2s
Cyano.MS2 <- read.csv("data_extra/Standards_MS2/Cyano_stds_withMS2s.csv") %>%
  mutate(column = "RP",
         z = NA)
HILICPos.MS2 <- read.csv("data_extra/Standards_MS2/HILICPos_stds_withMS2s.csv") %>%
  mutate(column = "HILIC",
         z = 1)
HILICNeg.MS2 <- read.csv("data_extra/Standards_MS2/HILICNeg_stds_withMS2s.csv") %>%
  mutate(column = "HILIC",
         z = -1)

# MS2s combined
Standards.MS2 <- Cyano.MS2 %>%
  rbind(HILICPos.MS2, HILICNeg.MS2) %>%
  rename(MS2 = MS2s) %>%
  select(-rtmin, -rtmax, -rt, -mz_obs, -mz)

# Standards csv
Ingalls.Standards <- read.csv("data_extra/Ingalls_Lab_Standards_Feb1.csv") %>%
  mutate(rt = (RT..min. * 60)) %>%
  rename(compound = Compound.Name_old,
         mz = m.z,
         column = Column) %>%
  select(compound, mz, rt, column, z)

# All theoretical values 
Theoretical.Values <- Ingalls.Standards %>% 
  left_join(Standards.MS2, by = c("compound", "column", "z")) %>%
  transform(mz = as.numeric(mz))

# Unknown (Experimental) Values -----------------------------------------------------
# Must be in the following format, including capitalization:
# mz (numeric), rt (numeric), column (numeric), z (numeric), MS2 (in concatenated format, character)

# Experimental.Values <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") %>%
#   separate(MassFeature_Column, into = c("MassFeature", "drop"), sep = "_") %>%
#   rename(column = Column,
#          KRH_identification = Identification) %>%
#   mutate(compound_experimental = 1:nrow(.)) %>%
#   select(compound_experimental, KRH_identification, mz, rt, column, z, MS2) 

Experimental.Values <- untargeted.metco %>%
  mutate(primary_key = 1:nrow(.)) %>%
  rename(MS2 = MS2dat1)
  
## Confidence Level 1 Matching ----------------------------------------
# My.Fuzzy.Join <- Theoretical.Values %>%
#   difference_left_join(Experimental.Values, by = c("mz"), max_dist = mz.flexibility) %>%
#   rename(compound_theoretical = compound,
#          mz_theoretical = mz.x,
#          rt_sec_theoretical = rt.x,
#          column_theoretical = column.x,
#          z_theoretical = z.x,
#          MS2_theoretical = MS2.x,
#          mz_experimental = mz.y,
#          rt_sec_experimental = rt.y,
#          column_experimental = column.y,
#          z_experimental = z.y,
#          MS2_experimental = MS2.y) %>%
#   select(compound_experimental, KRH_identification, compound_theoretical, mz_experimental, mz_theoretical, rt_sec_experimental, rt_sec_theoretical,
#          column_experimental, column_theoretical, z_experimental, z_theoretical, MS2_experimental, MS2_theoretical)  %>%
#   arrange(compound_experimental)

My.Fuzzy.Join <- Theoretical.Values %>%
  difference_left_join(Experimental.Values, by = c("mz"), max_dist = mz.flexibility) %>%
  rename(compound_theoretical = compound,
         mz_theoretical = mz.x,
         rt_sec_theoretical = rt.x,
         column_theoretical = column.x,
         z_theoretical = z.x,
         MS2_theoretical = MS2.x,
         mz_experimental = mz.y,
         rt_sec_experimental = rt.y,
         column_experimental = column.y,
         z_experimental = z.y,
         MS2_experimental = MS2.y) %>%
  select(primary_key, compound_theoretical, mz_experimental, mz_theoretical, rt_sec_experimental, rt_sec_theoretical,
         column_experimental, column_theoretical, z_experimental, z_theoretical, MS2_experimental, MS2_theoretical,
         previous_match, previous_rank)  %>%
  arrange(primary_key)


# # Confidence Level 1 ----------------------------------------  
# Confidence.Level.1 <- My.Fuzzy.Join %>%
#   filter(z_experimental == z_theoretical,
#          column_experimental == column_theoretical) %>%
#   mutate(mz_similarity_score = exp(-0.5 * (((mz_experimental - mz_theoretical) / mz.flexibility) ^ 2)),
#          rt_similarity_score = exp(-0.5 * (((rt_sec_experimental - rt_sec_theoretical) / rt.flexibility) ^ 2))) %>%
#   rowwise() %>%
#   mutate(ppm_mass_error = ((abs(mz_experimental - mz_theoretical)) / mz_theoretical) * 10^6) %>%
#   mutate(MS2_cosine_similarity = ifelse(is.na(MS2_experimental) | is.na(MS2_theoretical), 
#                                         NA, MS2CosineSimilarity(MakeScantable(MS2_experimental), MakeScantable(MS2_theoretical)))) %>%
#   mutate(total_similarity_score = ifelse(is.na(MS2_cosine_similarity), 
#                                          ((mz_similarity_score + rt_similarity_score) / 2) * 100,
#                                             ((MS2_cosine_similarity + mz_similarity_score + rt_similarity_score) / 3) * 100))

# Confidence Level 1 ----------------------------------------  
Confidence.Level.1 <- My.Fuzzy.Join %>%
  filter(z_experimental == z_theoretical,
         column_experimental == column_theoretical) %>%
  mutate(mz_similarity_score = exp(-0.5 * (((mz_experimental - mz_theoretical) / mz.flexibility) ^ 2)),
         rt_similarity_score = exp(-0.5 * (((rt_sec_experimental - rt_sec_theoretical) / rt.flexibility) ^ 2))) %>%
  rowwise() %>%
  mutate(ppm_mass_error = ((abs(mz_experimental - mz_theoretical)) / mz_theoretical) * 10^6) %>%
  mutate(MS2_cosine_similarity = ifelse(is.na(MS2_experimental) | is.na(MS2_theoretical),
                                        NA, MS2CosineSimilarity(MakeScantable(MS2_experimental), MakeScantable(MS2_theoretical)))) %>%
  mutate(total_similarity_score = ifelse(is.na(MS2_cosine_similarity),
                                         ((mz_similarity_score + rt_similarity_score) / 2) * 100,
                                         ((MS2_cosine_similarity + mz_similarity_score + rt_similarity_score) / 3) * 100))


# Sanity check -------------------------------------------------------------
# No fuzzy match (no mz within 0.02 daltons)
# No.Fuzzy.Match <- Experimental.Values %>%
#   select(compound_experimental) %>%
#   filter(compound_experimental %in% setdiff(1:nrow(.), My.Fuzzy.Join$compound_experimental)) %>%
#   pull()

No.Fuzzy.Match <- Experimental.Values %>%
  select(primary_key) %>%
  filter(primary_key %in% setdiff(1:nrow(.), My.Fuzzy.Join$primary_key)) %>%
  pull()

# # Fuzzy match, but wrong z/column
# No.CL1.Match <- setdiff(1:nrow(Experimental.Values), 
#                         sort(c(unique(Confidence.Level.1$compound_experimental), 
#                                No.Fuzzy.Match)))

No.CL1.Match <- setdiff(1:nrow(Experimental.Values),
                        sort(c(unique(Confidence.Level.1$primary_key),
                               No.Fuzzy.Match)))


# # Have any compounds been lost? Check for a TRUE output
# all.experimentals <- sort(c(No.Fuzzy.Match, No.CL1.Match, unique(Confidence.Level.1$compound_experimental)))
# length(all.experimentals) == length(Experimental.Values$compound_experimental)

# Have any compounds been lost? Check for a TRUE output
all.experimentals <- sort(c(No.Fuzzy.Match, No.CL1.Match, unique(Confidence.Level.1$primary_key)))
length(all.experimentals) == length(Experimental.Values$primary_key)

# Make "no match" dataframes for comparison
# No.CL1.Match.df <- Experimental.Values %>%
#   filter(compound_experimental %in% No.CL1.Match)
# No.Fuzzy.Match.df <- Experimental.Values %>%
#   filter(compound_experimental %in% No.Fuzzy.Match)

No.CL1.Match.df <- Experimental.Values %>%
  filter(primary_key %in% No.CL1.Match)
No.Fuzzy.Match.df <- Experimental.Values %>%
  filter(primary_key %in% No.Fuzzy.Match)

# Let's land on MARS ------------------------------------------------------
# Mission.Accomplished <- Confidence.Level.1 %>%
#   bind_rows(No.CL1.Match.df) %>%
#   bind_rows(No.Fuzzy.Match.df) %>%
#   mutate(confidence_rank = ifelse(mz_similarity_score > 0.9 & rt_similarity_score > 0.75 & ppm_mass_error < 7, 1, NA),
#          confidence_source = ifelse(!is.na(confidence_rank), "Ingalls_Standards", NA)) %>%
#   mutate(mz_experimental = ifelse(is.na(mz_experimental) & !is.na(mz), mz, mz_experimental),
#          rt_sec_experimental = ifelse(is.na(rt_sec_experimental) & !is.na(rt), rt, rt_sec_experimental),
#          column_experimental = ifelse(is.na(column_experimental) & !is.na(column), column, column_experimental),
#          z_experimental = ifelse(is.na(z_experimental) & !is.na(z), z, z_experimental),
#          MS2_experimental = ifelse(is.na(MS2_experimental) & !is.na(MS2), MS2, MS2_experimental)) %>%
#   select(-c("mz", "rt", "column", "z", "MS2")) %>%
#   arrange(compound_experimental)

Mission.Accomplished <- Confidence.Level.1 %>%
  bind_rows(No.CL1.Match.df) %>%
  bind_rows(No.Fuzzy.Match.df) %>%
  mutate(confidence_rank = ifelse(mz_similarity_score > 0.9 & rt_similarity_score > 0.75 & ppm_mass_error < 7, 1, NA),
         confidence_source = ifelse(!is.na(confidence_rank), "Ingalls_Standards", NA)) %>%
  mutate(mz_experimental = ifelse(is.na(mz_experimental) & !is.na(mz), mz, mz_experimental),
         rt_sec_experimental = ifelse(is.na(rt_sec_experimental) & !is.na(rt), rt, rt_sec_experimental),
         column_experimental = ifelse(is.na(column_experimental) & !is.na(column), column, column_experimental),
         z_experimental = ifelse(is.na(z_experimental) & !is.na(z), z, z_experimental),
         MS2_experimental = ifelse(is.na(MS2_experimental) & !is.na(MS2), MS2, MS2_experimental)) %>%
  select(-c("mz", "rt", "column", "z", "MS2")) %>%
  arrange(primary_key)

# Save your csv -----------------------------------------------------------
# Save your Mission.Accomplished dataframe for use in the next level.
#write.csv(Mission.Accomplished, "data_processed/confidence_level1.csv", row.names = FALSE)

write.csv(Mission.Accomplished, "data_processed/MetCo_confidence_level1.csv", row.names = FALSE)