library(tidyverse)
source("Functions.R")
options(digits = 6)

Area_Angie <- read.delim("data_from_lab_members/Area_0_202010201051.txt")

# Set header, filter unknowns ---------------------------------------
columns.to.drop <- c('Formula', 'Ontology', 'INCHIKEY', 'SMILES', 
                     'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
                     'MS.MS.spectrum', 'Post.curation.result',
                     'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
                     'm.z.matched', 'Dot.product', 'Reverse.dot.product')

headers.set <- SetHeader(Area_Angie) %>%
  select(1:143) %>%
  #select(-one_of(columns.to.drop)) %>%
  rename(Metabolite.Name = Metabolite.name)

classes.changed <- ChangeXClasses(headers.set)

# Rearrange data and combine to one dataframe -------------------------------------------------'
Area <- RearrangeDatasets(classes.changed, parameter = "Area.Value")

# Standardize dataset --------------------------------------------------
Area.Standardized <- StandardizeMetabolites(Area)
###
Ingalls.Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                              stringsAsFactors = FALSE, header = TRUE)


# cosine similarity -------------------------------------------------------
# Remember, you can get the score w/out MS2! The eventual full function should be able to do both.

# "Scan" is ms2, a tiny dataframe of ms2 where the first column is mass and the second is intensity.
# The below function only works on one comparison at a time, there is another function in KRH file
# for making a dataframe of ms2. Keep digging in there.

## Currently only looking at Uracil
Uracil.standard.msp <- read.delim("data_extra/Ingalls_HILICNeg_Standards.msp", header = FALSE, sep = "") %>%
  slice(56:66) %>%
  rename(mz = 1, intensity = 2)

#Pull an MS2 spectra from a DDA file according to the a targeted mass and time------
#Sums intensities and filters out fragments like in Tabb2003
getms2spectra <-
  function(xs=xs, mass=mass, time=time, timetol=15, masstolLR=0.4, masstolHR=0.02){
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
#Close Funtion

#Get a scan table from a concatenated scanlist---- (Scan is output from getms2spectra function)

scantable <- function(Scan) {
  Try <- read.table(text=as.character(Scan),
                    col.names=c('mz','intensity')) %>% 
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", "")))
  return(Try)
}     

MSMScosine_1 <- function(scan1, scan2, mass1, mass2) {
  
  mztolerance <- 0.02
  
  w1 <- (scan1[, 1] ^ 2) * sqrt(scan1[, 2])
  w2 <- (scan2[, 1] ^ 2) * sqrt(scan2[, 2])
  
  diffmatrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x)
  sameindex <- which(abs(diffmatrix) < mztolerance, arr.ind = T)
  
  similarity <- sum(w1[sameindex[, 2]] * w2[sameindex[, 1]]) / (sqrt(sum(w2^2)) * sqrt(sum(w1^2)))
  
  return(similarity)
}

MSMSconsine1_df <- function(df){
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




## NEED TO FIX THE COMPOUND NAME REPLACEMENTS

# update.names <- Area.Standardized %>%
#   select(Metabolite.Name) %>%
#   rename(Compound.Name_old = Metabolite.Name) %>%
#   left_join(Ingalls.Standards %>% select(Compound.Name, Compound.Name_old)) %>%
#   rename(Compound.Name_new = Compound.Name,
#          Metabolite.Name = Compound.Name_old)
# 
# combined.final <- Area.Standardized %>%
#   left_join(update.names) %>%
#   select(-Metabolite.Name) %>%
#   rename(Metabolite.Name = Compound.Name_new) %>%
#   unique()
