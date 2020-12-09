#Check against the KEGG list of compounds -----


KEGGCompounds_withMasses <- read.csv("data_extra/KEGGCompounds_withMasses.csv", header = TRUE)
CommonContams <- read.csv("data_extra/CommonContams.csv", header = TRUE)

checkKEGG <- function(MFs, ppmtol) { # Guessing MFs is the mass features from the experiment? Modified version of the MF clusters?
  if(missing(ppmtol)) {
    ppmtol <- 5
  }
  matchedKEGGs <- list()
  keggCompounds <- read.csv(text = getURL("https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/KEGGCompounds_withMasses.csv"), header = T) %>% rename(Compound= OtherCmpds)
  keggPos <- keggCompounds %>% select(Compound, PosMZ) %>% 
    mutate(Fraction1 = "HILICPos", Fraction2="CyanoAq", Fraction3= "CyanoDCM") %>% unique() %>% 
    rename(mz = PosMZ)
  keggCompoundsShort <- keggCompounds %>% select(Compound, NegMZ) %>% 
    mutate(Fraction1 = "HILICNeg", Fraction2=NA, Fraction3= NA) %>% unique() %>%
    rename(mz = NegMZ) %>% rbind(keggPos)
  MFstry <- MFs %>% mutate(MF_Frac2 = MF_Frac) %>% separate(MF_Frac2, c("MF", "Frac"), sep =  "_") %>% select(-MF)
  matched <- difference_inner_join(x= keggCompoundsShort, y = MFstry, 
                                   by = "mz", max_dist = .01,  distance_col = NULL) %>%
    mutate(ppm = (abs(mz.x-mz.y )/mz.x *10^6)) %>%
    filter(ppm < ppmtol, Frac == Fraction1 | Frac == Fraction2 | Frac == Fraction3) %>%
    mutate(Frac = as.factor(Frac),mz = mz.x) %>% 
    select(MF_Frac, Compound, ppm ) %>% rename(Keggppm = ppm)
  
  matchedNames <- left_join(matched, keggCompounds) %>% select(Compound, Name) %>%
    group_by(Compound) %>%
    summarise(KEGGMatchesNames = as.character(paste(Name,  collapse=" ")))
  matchedKEGGs[[1]] <- left_join(matched, matchedNames) %>% rename(KeggMatches = Compound)
  matchedKEGGs[[2]] <- matchedKEGGs[[1]] %>%
    arrange(Keggppm) %>%
    group_by(MF_Frac) %>%
    summarise(KeggMatches = as.character(paste(KeggMatches,  collapse="; ")),
              Keggppm = as.character(paste(Keggppm,  collapse="; ")),
              KeggNames = as.character(paste(KEGGMatchesNames,  collapse="; ")))
  return(matchedKEGGs)}



#flag for known Contaminants -------
checkContaminants <- function(MFs, ppmtol) {
  if(missing(ppmtol)) {
    ppmtol <- 15
  }
  matchedKnownContaminants <- list()
  knownContaminants <- read.csv(text = getURL("https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/CommonContams.csv"), comment.char = "#", header = T)%>%
    mutate(mz = m.z) %>% select(Fraction1:Fraction3, mz, Compound) %>% mutate(Flag = "PossibleContamin")
  MFstry <- MFs %>% mutate(MF_Frac2 = MF_Frac) %>% separate(MF_Frac2, c("MF", "Frac"), sep =  "_") %>% select(-MF)
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
