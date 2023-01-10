IsolateMoNACandidates <- function(MoNA.Mass, experimental.df) {
  # Compare experimental mass features to scraped MoNA data,
  # based on MS1 and MS2 information.
  #
  # Args
  #   MoNA.Mass: Single MoNA mass, isolated from scraped MoNA df.
  #   experimental.df: experimental dataframe, islated to contain 
  #                    only MS1 and MS2 data.
  # 
  # Returns
  #   final.candidates: dataframe of experimental data, matched with
  #                     information from MoNA.
  #
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

MakeScantable <- function(scan) {
  # Create a filtered MS2 scantable from a concatenated scanlist of MS2s.
  #
  # Args
  #   scan: Single observation from a dataframe containing MS2 mz and intensity 
  #         spectra, separated by semicolons (;).
  #
  # Returns
  #   scantable: Tiny dataframe, containing columns of mz and intensity. 
  #              Intensity is scaled to 100 and filtered to drop all intensity
  #              values below 0.5.
  #
  scantable <- read.table(text = as.character(scan),
                          col.names = c("mz", "intensity"), fill = TRUE) %>%
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", "")),
           intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
    filter(intensity > 0.5) %>%
    arrange(desc(intensity))
  
  return(scantable)
}

MS2CosineSimilarity <- function(scan1, scan2) {
  # Finds the weighted cosine similarity between two sets of MS2 spectra.
  #
  # Args
  #   scan1 & scan2: Tiny dataframes of MS2. First column is mz, second column is intensity.
  #                  These are outputs of the MakeScantable() function.
  #
  # Returns
  #   cosine.similarity: A weighted similarity score between 0 and 1, indicating the cosine
  #                      relationship of the two vectors.
  #   
  mz.tolerance <- 0.02
  
  weight1 <- (scan1[, 1] ^ 2) * sqrt(scan1[, 2])
  weight2 <- (scan2[, 1] ^ 2) * sqrt(scan2[, 2])
  
  diff.matrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x) 
  same.index <- which(abs(diff.matrix) < mz.tolerance, arr.ind = TRUE)
  cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) / 
    (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))
  
  return(cosine.similarity)
}

MakeMS2CosineDataframe <- function(df) {
  scan1 <- MakeScantable(df["scan1"])
  scan2 <- MakeScantable(df["scan2"])
  mass1 <- df["mass1"]
  mass2 <- df["mass2"]
  
  mz.tolerance <- 0.02
  
  weight1 <- (scan1[, 1] ^ 2) * sqrt(scan1[, 2])
  weight2 <- (scan2[, 1] ^ 2) * sqrt(scan2[, 2])
  
  difference.matrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x)
  same.index <- which(abs(difference.matrix) < mz.tolerance, arr.ind = TRUE)

  similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) /
    (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))
  
  return(similarity)
}

ReplaceNA <- function(column) {
  gsub("NA; ", "", column)
}
