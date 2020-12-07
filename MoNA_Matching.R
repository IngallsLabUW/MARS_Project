library(fuzzyjoin)
library(parallel)
library(tidyverse)
source("Functions.R")
source("KRH_Functions.R")


### Takes forever. Parallelization is not working right now.
### Creating many duplicate rows. Join explosion is happening.
MSMScosine1_df <- function(df) {
  scan1 <- scantable(df["scan1"])
  scan2 <- scantable(df["scan2"])
  mass1 <- df["mass1"]
  mass2 <- df["mass2"]
  print("scan1")
  print(scan1)
  print("scan2")
  print(scan2)
  
  mztolerance<-0.02
  
  w1<-(scan1[,1]^2)*sqrt(scan1[,2])
  w2<-(scan2[,1]^2)*sqrt(scan2[,2])
  
  diffmatrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x)
  print("diffmatrix")
  print(diffmatrix)
  sameindex <- which(abs(diffmatrix)<mztolerance,arr.ind=T)
  
  similarity<-sum(w1[sameindex[,2]]*w2[sameindex[,1]])/(sqrt(sum(w2^2))*sqrt(sum(w1^2)))
  
  return(similarity)
}

scantable <- function(Scan) {
  Try <- read.table(text=as.character(Scan),
                    col.names = c("mz","intensity"), fill = TRUE) %>% 
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", ""))) %>%
    mutate(intensity = round(intensity/max(intensity)*100, digits = 1)) %>%
    filter(intensity > 0.5) %>%
    arrange(desc(intensity))
  
  return(Try)
}  

Cosine.Score.Cutoff <- 0.5 
MassBank.ppm.Cutoff <- 5

# One of four relational spreadsheets from the MoNA download
MoNA.Spectra <- read.csv("data_extra/NEG_Spectra.csv") 

Experimental.Spectra <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") %>% 
  rename(MF.Fraction = MassFeature_Column) %>%
  select(MF.Fraction, mz, MS2) %>%
  as.data.frame()

# Reassign variables from Experimental.Spectra
mz <- as.numeric(unlist(Experimental.Spectra["mz"])) 
MS2 <- as.character(Experimental.Spectra["MS2"]) 
MF.Fraction <- as.data.frame(Experimental.Spectra["MF.Fraction"]) 

Experimental.Spectra.ForJoin <- Experimental.Spectra %>% 
  rename(MH_mass = mz) %>%
  mutate(No.MS2 = ifelse(is.na(MS2), TRUE, FALSE)) %>%
  filter(No.MS2 == FALSE) %>%
  select(-No.MS2)

# Subtract hydrogen for reference database
MoNA.Spectra.MHMass <- MoNA.Spectra %>%
  mutate(MH_mass = M_mass - 1.0072766) %>%
  select(ID, Names, spectrum_KRHform_filtered, MH_mass) 
MoNA.Spectra.MHMass[MoNA.Spectra.MHMass==""]<-NA
MoNA.Spectra.MHMass <- MoNA.Spectra.MHMass %>%
  drop_na()


## Example masses for testing
testmass1 <- 116.07093 
testmass2 <- 74.02475

testdf1 <- MoNA.Spectra.MHMass %>%
  filter(MH_mass > testmass1 - 0.020,  
         MH_mass < testmass1 + 0.020) 
testdf2 <- MoNA.Spectra.MHMass %>%
  filter(MH_mass > testmass2 - 0.020,  
         MH_mass < testmass2 + 0.020) 
testdf <- testdf1 %>% rbind(testdf2) %>%
  arrange(MH_mass)

## Full function
MakeCandidates <- function(MoNA.Mass) {
  Candidates <- MoNA.Spectra.MHMass %>% 
    filter(MH_mass > MoNA.Mass - 0.020,  
           MH_mass < MoNA.Mass + 0.020) %>% 
    difference_inner_join(Experimental.Spectra.ForJoin, by = "MH_mass", max_dist = 0.02) %>% 
    mutate(scan1 = spectrum_KRHform_filtered, # scan1 is MS2 from MoNA
           scan2 = MS2, # scan2 is MS2 from the experimental data
           mass1 = MH_mass.x, # mass1 is the primary mass from MoNA
           mass2 = MH_mass.y) # mass2 is the primary mass from experimental data
  
  if (length(Candidates$ID) == 0) {
    print("There are no potential candidates.")
    No.Match.Return <- MF.Fraction %>%
      mutate(MassBankMatch = NA,
             MassBankppm = NA,
             MassBankCosine1 = NA)
    
    return(No.Match.Return)
  }
  
  # Add cosine similarity scores
  print("Making potential candidates")
  
  Candidates$Cosine1 <- apply(Candidates, 1, FUN=function(x) MSMScosine1_df(x)) 
  
  Candidates.Filtered.Cosine <- Candidates %>%
    filter(Cosine1 > Cosine.Score.Cutoff) %>%
    arrange(desc(Cosine1))
  
  Final.Candidates <- Candidates.Filtered.Cosine %>%
    mutate(MassBankMatch = paste(Names, ID, sep = " ID:"),
           MassBankppm = abs(mass1 - mass2) / mass2 * 10^6,
           MassBankCosine1 = Cosine1) %>%
    unique() %>%
    filter(MassBankppm < MassBank.ppm.Cutoff) #%>% 
  #select(MF.Fraction, MassBankMatch:MassBankCosine1)
  
  return(Final.Candidates)
}

output <- lapply(unique(MoNA.Spectra.MHMass$MH_mass), MakeCandidates)
outputdf <- bind_rows(output)

#test <- mclapply(unique(MoNA.Spectra.MHMass$MH_mass), MakeCandidates, detectCores())

