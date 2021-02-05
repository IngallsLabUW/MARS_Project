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

# Total Similarity Score
# TS = ((MS2 Similarity + MS1 Similarity) / 2) * 100
Uracil.Total.Similarity_AllVariables <- ((uracil.MS2cosine.sim + MS1.mz.similarity + MS1.rt.similarity) / 3) * 100
Uracil.Total.Similarity_NoMS2 <- ((MS1.mz.similarity + MS1.rt.similarity) / 2) * 100
