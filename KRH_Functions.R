library(tidyverse)

## Modified KRH functions 

#Check known compounds
checkKnownCompounds <- function(MFs, ppmtol, rttolRP, rttolHILIC) {
  if(missing(ppmtol)) {
    ppmtol <- 15
  }
  if(missing(rttolRP)) {
    rttolRP <- 0.5
  }
  if(missing(rttolHILIC)) {
    rttolHILIC <- 1
  }
  matchedKnownCompounds <- list()
  knownShortCompounds <- read.csv(text = getURL("https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"), header = T)  %>%
    mutate(mz = m.z, RT = RT..min.) %>% select(Compound.Name, mz, RT, Fraction1, Fraction2) %>% 
    mutate(Fraction1 = as.character(Fraction1),
           Fraction2 = as.character(Fraction2))
  MFstry <- MFs %>% mutate(MF_Frac2 = MF_Frac) %>% separate(MF_Frac2, c("MF", "Frac"), sep =  "_") %>% select(-MF)  
  matchedKnownCompounds[[1]] <- difference_inner_join(x= MFstry, y = knownShortCompounds, 
                                                      by = "mz", max_dist = .01,  distance_col = NULL) %>%
    mutate(MZdiff = abs(mz.x-mz.y ))%>%
    mutate(RTdiff = abs(RT.x-RT.y),
           ppm = (MZdiff/mz.x *10^6)) %>%
    filter(ppm < ppmtol) %>%
    filter(Frac == Fraction1 | Frac == Fraction2) %>%
    mutate(Frac = as.factor(Frac), mz = mz.x) %>%
    filter(!((Frac == "CyanoAq" | Frac == "CyanoDCM") &  RTdiff > rttolRP)) %>%
    filter(!((Frac == "HILICNeg" | Frac == "HILICPos") &  RTdiff > rttolHILIC))   %>%
    select(MF_Frac, mz, rt, Compound.Name, RTdiff, ppm)
  matchedKnownCompounds[[2]] <- matchedKnownCompounds[[1]] %>% 
    arrange(RTdiff) %>%
    group_by(MF_Frac) %>%
    summarise(TargetMatches = as.character(paste(Compound.Name,  collapse="; ")),
              ppmMatches = as.character(paste(ppm,  collapse="; ")),
              RTdiffMatches = as.character(paste(RTdiff,  collapse="; ")))
  return(matchedKnownCompounds)
  
}

#Pull an MS2 spectra from a DDA file according to the a targeted mass and time------
#Sums intensities and filters out fragments like in Tabb2003
getms2spectra <- function(xs=xs, mass=mass, time=time, timetol=15, masstolLR=0.4, masstolHR=0.02){
  #xs is msn2xcmsRaw(xcmsRaw(DDAFILE, includeMSn=TRUE))
  masslow<-mass-masstolLR
  masshigh<-mass+masstolLR
  masslowhr <- mass - masstolHR
  masshighr <- mass + masstolHR
  timelow<-time-timetol
  timehigh<-time+timetol
  scanlist<-which(xs@msnPrecursorMz>masslow & xs@msnPrecursorMz<masshigh & xs@scantime>timelow & xs@scantime<timehigh)
  allMS2list <- list()
  if(length(scanlist)>0){
    for (l in 1:length(scanlist)){
      ms2dat <- as.data.frame(getScan(xs,scanlist[l])) %>% mutate(scannumber = scanlist[l])
      allMS2list[[l]] <- ms2dat} #This loop gets all matching scans with LR
    allMS2_matched <- do.call(rbind, allMS2list) %>%
      filter(mz > masslowhr & mz < masshighr) 
    matchedscanlist <- allMS2_matched$scannumber} else {matchedscanlist <- c()} #Do we find any HR matches?  
  if(length(matchedscanlist > 0)){
    MS2Summed <- do.call(rbind, allMS2list) %>%
      filter(scannumber %in% allMS2_matched$scannumber) %>%
      mutate(mzRND = as.factor(round(mz, digits = 2))) %>%
      group_by(mzRND) %>%
      summarise(mz = mean(mz),
                intensity = sum(intensity),
                sqrtintensity = sqrt(intensity),
                scancount = n())
    bestscan <- MS2Summed %>%
      arrange(desc(intensity)) %>%
      mutate(intensity = round(intensity/max(intensity)*100, digits = 1))%>%
      filter(intensity > 0.5) %>%
      mutate(mz = round (mz, digits = 5),
             mash = paste(mz, intensity, sep = ", " ))
    sortedscanrange <- paste(bestscan$mash, collapse = "; ")
    return(sortedscanrange)}else{return(NA)} #Close if/then for scanlist >1
  
}

#The way that I implement the cosine function is that I only do the comparison after mass1 and mass2. 
# See function massbankMS2MatchPOSITIVE for an example. 
massbankMS2MatchPOSITIVE <- function(ShortestDat) { # ShortestDat is the unknowns, the items that are being annotated. 
  # can be taken from MF cluster assignments, MFfrac is the compound ID
  mz <- as.numeric(ShortestDat["mz"])
  MS2 <- as.character(ShortestDat["MS2"])
  MF_Frac <- as.character(ShortestDat["MF_Frac"])
  
  Candidates <- Spectra %>% ## possible matches in mz space.
    # only do ms2 matches on possible ms1 matches. ok on standards becuase its short, but need a filter for larger dbs. 
    mutate(MH_mass = M_mass + 1.0072766) %>%
    filter(near(MH_mass, mz, tol= 0.02)) %>%
    mutate(scan1 = spectrum_KRHform_filtered, 
           scan2 = MS2, mass1 = MH_mass, mass2 = mz)%>% # Sets up dataframe to do comparison
    filter(!is.na(scan1))
  
  NoMatchReturn <- Spectra %>% # What happens if there are no matches. NAs go back to original dataset
    mutate(MassBankMatch = NA,
           MassBankppm = NA,
           MassBankCosine1 = NA,
           MF_Frac = MF_Frac) %>%
    select(MF_Frac, MassBankMatch:MassBankCosine1) %>% head(1)
  
  # some filters below may be thrown out
  if (length(Candidates$ID > 1)){ # if there is a match, 
    Candidates$Cosine1 <- apply(Candidates, 1, FUN=function(x) MSMSconsine1_df(x)) # applies cosine function
    Candidates$Cosine2 <- apply(Candidates, 1, FUN=function(x) MSMSconsine2_df(x))
    Candidates <- Candidates %>% filter(Cosine1 > 0.8) %>% arrange(desc(Cosine1)) %>% 
      mutate(MF_Frac = MF_Frac) %>% head(1)
    Candidates <- Candidates %>% 
      mutate(MassBankMatch = paste(Names, ID, sep= " ID:"),
             MassBankppm = abs(mass1-mass2)/mass2 *10^6,
             MassBankCosine1 = Cosine1) %>%
      filter(MassBankppm < 5) %>%
      select(MF_Frac, MassBankMatch:MassBankCosine1)
  }
  if (length(Candidates$MF_Frac) == 0){
    Candidates <-NoMatchReturn
  }
  return(Candidates)
}

# Finds the cosine similarity method for comparing two mass spectra.
MSMScosine_1 <- function(scan1, scan2) {
  # Finds the cosine similarity method for comparing two mass spectra.
  #
  # Args
  #   scan1 & scan2: Tiny dataframes of ms2, first column is mz, second column is intensity.
  #
  # Returns
  #   similarity: A similarity score between 0 and 1, indicating the similarity of the two vectors.
  #   
  mztolerance <- 0.02
  
  w1 <- (scan1[, 1] ^ 2) * sqrt(scan1[, 2])
  w2 <- (scan2[, 1] ^ 2) * sqrt(scan2[, 2])
  
  diffmatrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x) 
  sameindex <- which(abs(diffmatrix) < mztolerance, arr.ind = T)
  
  similarity <- sum(w1[sameindex[, 2]] * w2[sameindex[, 1]]) / (sqrt(sum(w2 ^ 2)) * sqrt(sum(w1 ^ 2)))
  
  return(similarity)
}

#Check against MassBank Spectrum (m/z and spectra match for putative ID)
# I'd look at MSMScosine_1_df to see how you can feed in a df of mass1, mass2, scan1 (in my format), scan2 and get an output of df$cosine.
MSMScosine1_df <- function(df) {
  scan1 <- scantable(df["scan1"])
  scan2 <- scantable(df["scan2"])
  mass1 <- df["mass1"]
  mass2 <- df["mass2"]
  
  mztolerance<-0.02
  
  w1<-(scan1[,1]^2)*sqrt(scan1[,2])
  w2<-(scan2[,1]^2)*sqrt(scan2[,2])
  
  diffmatrix<-sapply(scan1[,1], function(x) scan2[,1]-x)
  sameindex<-which(abs(diffmatrix)<mztolerance,arr.ind=T)
  
  similarity<-sum(w1[sameindex[,2]]*w2[sameindex[,1]])/(sqrt(sum(w2^2))*sqrt(sum(w1^2)))
  
  return(similarity)
}

MSMScosine2_df <- function(df){
  scan1 <- scantable(df["scan1"])
  scan2 <- scantable(df["scan2"])
  mass1 <- as.numeric(df["mass1"])
  mass2 <- as.numeric(df["mass2"])
  
  mztolerance<-0.02
  
  loss1<-scan1
  loss1[,1] <- mass1 - loss1[, 1]
  loss2<-scan2
  loss2[,1]<-mass2-loss2[,1]
  w1<-(loss1[,1]^2)*sqrt(loss1[,2])
  w2<-(loss2[,1]^2)*sqrt(loss2[,2])
  
  diffmatrix<-sapply(loss1[,1], function(x) loss2[,1]-x)
  sameindex<-which(abs(diffmatrix)<mztolerance,arr.ind=T)
  
  similarity2<-sum(w1[sameindex[,2]]*w2[sameindex[,1]])/(sqrt(sum(w2^2))*sqrt(sum(w1^2)))
  
  return(similarity2)
}

# Get a scan table from a concatenated scanlist---- (Scan is output from getms2spectra function)
scantable <- function(Scan) {
  Try <- read.table(text=as.character(Scan),
                    col.names=c('mz','intensity')) %>% 
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", ""))) %>%
    mutate(intensity = round(intensity/max(intensity)*100, digits = 1))%>%
    filter(intensity > 0.5) %>%
    arrange(desc(intensity)) 
  
  return(Try)
}    

