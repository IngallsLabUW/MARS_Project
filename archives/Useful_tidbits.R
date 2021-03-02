## Useful tidbits

A1MS2 <- MyFuzzyJoin %>% # From Annotate_Confidence_Level1
  filter(Unknown.Compound %in% A1Confidence$Unknown.Compound) %>%
  select(Unknown.Compound, Compound_Standards, MS2_Unknowns, MS2_Standards) 

tosplit <- A1MS2 %>%
  group_by(Unknown.Compound) %>% 
  group_split()

scantable_Unknowns <- lapply(tosplit, function(df) {
  lapply(unique(df[, 3]), MakeScantable)
})

scantable_Standards <- lapply(tosplit, function(df){
  apply(df[4], 1, MakeScantable)
})


## Confidence Level 1 MS2 matching
Uracil <- MyFuzzyJoin %>%
  filter(Compound_Standards == "Uracil") %>%
  select(Compound_Standards, MS2_Unknowns, MS2_Standards)

Uracil.Standards <- Uracil %>%
  select(MS2_Standards) %>%
  slice(1) %>%
  MakeScantable()

Uracil.Experimental <- Uracil %>%
  select(MS2_Unknowns) %>%
  slice(1) %>%
  MakeScantable()

uracil.MS2cosine.sim <- MS2CosineSimilarity(scan1 = Uracil.Experimental, scan2 = Uracil.Standards)

#group_by(compound_known) %>%
# filter(!rt_seconds_known < (rt_seconds_unknown - rt.flexibility) & 
#          !rt_seconds_known > (rt_seconds_unknown + rt.flexibility)) %>% 
#ungroup() %>%

# Total Similarity Score
# TS = ((MS2 Similarity + MS1 Similarity) / 2) * 100
Uracil.Total.Similarity_AllVariables <- ((uracil.MS2cosine.sim + MS1.mz.similarity + MS1.rt.similarity) / 3) * 100
Uracil.Total.Similarity_NoMS2 <- ((MS1.mz.similarity + MS1.rt.similarity) / 2) * 100

## Confidence Level A3 ----------------------------------------

within10duplicate <- function(df, column) {
  if (column > 8) # replace with "more than one unique match" rather than a number
    # Also "if" statement is probably unnecessary
    df2 <- df %>%
      mutate(RT.diff = abs(RT.seconds_Unknowns - RT.seconds_Standards)) %>%
      group_by(Unknown.Compound) %>%
      mutate(Closest.rt.Match = min(RT.diff)) %>%
      mutate(Closest.rt.Match2 = ifelse(Closest.rt.Match == RT.diff, TRUE, FALSE))
  else df
}

A3Confidence_MS2s <- My.Fuzzy.Join %>% # mz and 10 RT
  filter(!Unknown.Compound %in% unique(Confidence.Level.1$Unknown.Compound),
         !Unknown.Compound %in% unique(A1Confidence$Unknown.Compound),
         !Unknown.Compound %in% unique(A2Confidence_MS2s$Unknown.Compound),
         !Unknown.Compound %in% unique(A2Confidence$Unknown.Compound)) %>%
  filter(z_Unknowns == z_Standards,
         Column_Unknowns == Column_Standards) %>%
  group_by(Unknown.Compound) %>%
  add_tally() %>%
  within10duplicate(column = "n") %>%
  ungroup() %>%
  mutate(mz_Similarity = exp(-0.5 * (((mz_Unknowns - mz_Standards) / 0.02) ^ 2)),
         RT_Similarity = exp(-0.5 * (((RT.seconds_Unknowns - RT.seconds_Standards) / 0.04) ^ 2))) %>% # Why all the 0s?
  filter_at(vars(MS2_Unknowns, MS2_Standards),all_vars(!is.na(.))) %>%
  rowwise() %>% 
  mutate(MS2cosinesim = MS2CosineSimilarity(MakeScantable(MS2_Unknowns), MakeScantable(MS2_Standards))) %>%
  mutate(Total.Similarity.Score = ((MS2cosinesim + mz_Similarity + RT_Similarity) / 3) * 100)


IsolateMoNACandidates <- function(MoNA.Mass) {
  # Create a dataframe of potential compound matches from the downloaded MoNA database.
  #
  # Args
  #   MoNA.Mass: Individual dataframe value of mass, minus hydrogen, 
  #              taken from the MoNA Spectra relational csvs. 
  #
  # Returns
  #   final.candidates: Dataframe of potential matches for experimental features.
  potential.candidates <- MoNA.Spectra.MHMass %>% 
    filter(MH_mass > MoNA.Mass - 0.020,  
           MH_mass < MoNA.Mass + 0.020) %>% 
    difference_inner_join(Experimental.Spectra.ForJoin, by = "MH_mass", max_dist = 0.02) %>% 
    mutate(scan1 = spectrum_KRHform_filtered, # scan1 is MS2 from MoNA
           scan2 = MS2,                       # scan2 is MS2 from the experimental data
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
