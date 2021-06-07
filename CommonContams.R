
CommonContams <- read.csv("data_extra/CommonContams.csv", header = FALSE) %>%
  slice(-1) %>%
  row_to_names(row = 1) 

#flag for known Contaminants -------
checkContaminants <- function(MFs, ppmtol) {
  if(missing(ppmtol)) {
    ppmtol <- 15
  }
  ppmtol <- 15 # temp line for function run through 
  matchedKnownContaminants <- list()
  knownContaminants <- CommonContams %>% 
    mutate(mz = m.z) %>% 
    mutate(mz = as.numeric(mz)) %>%
    select(Fraction1:Fraction3, mz, Compound) %>% 
    mutate(Flag = "PossibleContamin")
  MFstry <- MFs %>% 
    mutate(MF_Frac2 = MF_Frac) %>% 
    separate(MF_Frac2, c("MF", "Frac"), sep =  "_X_") %>% # modified sep  
    select(-MF)
  matchedKnownContaminants[[1]] <- difference_inner_join(x= knownContaminants, y = MFstry, 
                                                         by = "mz", max_dist = .01,  distance_col = NULL) %>%
    mutate(ppm = (abs(mz.x-mz.y )/mz.x *10^6)) %>%
    filter(ppm < ppmtol, Frac == Fraction1 | Frac == Fraction2 | Frac == Fraction3) %>%
    mutate(Frac = as.factor(Frac), mz = mz.x) %>% 
    select(MF_Frac, Compound, ppm, Flag ) %>% 
    rename(ContamMatches = Compound, Contamppm = ppm, ContamFlag = Flag)
  
  matchedKnownContaminants[[2]] <- matchedKnownContaminants[[1]] %>% 
    group_by(MF_Frac) %>%
    summarise(ContamMatches = as.character(paste(ContamMatches,  collapse="; ")),
              Contamppm = as.character(paste(Contamppm,  collapse="; ")))
  return(matchedKnownContaminants)
}

mytest2 <- checkContaminants(MFs.frame)
