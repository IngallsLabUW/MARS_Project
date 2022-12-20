library(RJSONIO)
# library(data.table)
# library(plyr)
library(tidyverse)
# library(devtools)
library(tidyjson)

#No scraping necessary - the .json files are found here: https://mona.fiehnlab.ucdavis.edu/downloads

#MoNAraw <- fromJSON("data_extra/MoNA-export-LC-MS-MS_Negative_Mode.json")
# MoNAraw <- fromJSON("data_extra/MoNA-export-LC-MS-MS_Positive_Mode.json")
# MoNAraw <- fromJSON("MoNADat/MoNA-export-HMDB.json") #getsPosandNeg

MoNAraw <- jsonlite::fromJSON("data_extra/MoNA-export-LC-MS_Spectra.json", simplifyVector = FALSE) %>%
  spread_all()

FullMoNA <- MoNAraw %>%
  select(id, spectrum, score = score.score, library = library.tag.text, ..JSON)

testMoNA <- FullMoNA %>%
  select(id, ..JSON) %>%
  filter(id == "AC000001") %>%
  unnest(cols = c(..JSON))
MetaDat <- sapply(testMoNA, function(x) x["metaData"])
MetaDat[["category"]] <- sapply(MetaDat[[1]], function(x) x["category"])
MetaDatTemp[["computed"]] <- sapply(MetaDat[[i]], function(x) x[["computed"]])
MetaDatTemp[["hidden"]] <- sapply(MetaDat[[i]], function(x) x[["hidden"]])
MetaDatTemp[["name"]] <- sapply(MetaDat[[i]], function(x) x[["name"]])
MetaDatTemp[["value"]] <- sapply(MetaDat[[i]], function(x) x[["value"]])


# Once run recall here -------
#Spectra <- read_csv("POS_Spectra.csv")
#MetaDatCMPDDFALL <- read_csv("POS_CmpInfo.csv")
#MetaDatDFALL <- read_csv("POS_MetaData.csv")
#MetaDatNamesDFALL <- read_csv("POS_Names.csv")

#Extract MoNA IDs-----
SpectrumID <- sapply(MoNAraw, function(x) x[["id"]])

#Extract MoNA MS2s------
Spectra <- sapply(MoNAraw, function(x) x[["spectrum"]])

#Extract MoNA score-------
SpectraScore <- sapply(MoNAraw, function(x) x[["score"]][["score"]])

#Get the database its from-------
SpectraDatabase <- sapply(MoNAraw, function(x) x[["library"]][["tag"]][["text"]])

#Get all the MetaData for each spectra------
MetaDat <- sapply(MoNAraw, function(x) x[["metaData"]])
MetaDatDFs <- list()
for (i in 1:length(MetaDat)) {
  MetaDatTemp <- list()
  MetaDatTemp[["category"]] <- sapply(MetaDat[[i]], function(x) x[["category"]])
  MetaDatTemp[["computed"]] <- sapply(MetaDat[[i]], function(x) x[["computed"]])
  MetaDatTemp[["hidden"]] <- sapply(MetaDat[[i]], function(x) x[["hidden"]])
  MetaDatTemp[["name"]] <- sapply(MetaDat[[i]], function(x) x[["name"]])
  MetaDatTemp[["value"]] <- sapply(MetaDat[[i]], function(x) x[["value"]])
  MetaDatTempDF <- do.call(cbind, MetaDatTemp) %>% as.data.frame()
  colnames(MetaDatTempDF) <- names(MetaDatTemp)
  MetaDatTempDF$SpectraID <- SpectrumID[i]
  MetaDatDFs[[i]] <-  MetaDatTempDF
}
MetaDatDFALL <- ldply(MetaDatDFs, rbind) %>% as.data.frame()

#Get all the MetaData for each compound------
MetaDat <- sapply(MoNAraw, function(x) x[["compound"]][[1]][["metaData"]])
MetaDatDFs <- list()
for (i in 1:length(MetaDat)) {
  MetaDatTemp <- list()
  MetaDatTemp[["category"]] <- sapply(MetaDat[[i]], function(x) x[["category"]])
  MetaDatTemp[["computed"]] <- sapply(MetaDat[[i]], function(x) x[["computed"]])
  MetaDatTemp[["hidden"]] <- sapply(MetaDat[[i]], function(x) x[["hidden"]])
  MetaDatTemp[["name"]] <- sapply(MetaDat[[i]], function(x) x[["name"]])
  MetaDatTemp[["value"]] <- sapply(MetaDat[[i]], function(x) x[["value"]])
  MetaDatTempDF <- do.call(cbind, MetaDatTemp) %>% as.data.frame()
  colnames(MetaDatTempDF) <- names(MetaDatTemp)
  MetaDatTempDF$SpectraID <- SpectrumID[i] ## error here
  MetaDatDFs[[i]] <-  MetaDatTempDF
}
MetaDatCMPDDFALL <- ldply(MetaDatDFs, rbind) %>% as.data.frame()
  
#Get all Names for each compound-----
MetaDat <- sapply(MoNAraw, function(x) x[["compound"]][[1]][["names"]])
MetaDatDFs <- list()
for (i in 1:length(MetaDat)) {
  MetaDatTemp <- list()
  MetaDatTemp[["computed"]] <- sapply(MetaDat[[i]], function(x) x[["computed"]])
  MetaDatTemp[["name"]] <- sapply(MetaDat[[i]], function(x) x[["name"]])
  MetaDatTemp[["score"]] <- sapply(MetaDat[[i]], function(x) x[["score"]])
  MetaDatTempDF <- do.call(cbind, MetaDatTemp) %>% as.data.frame()
  colnames(MetaDatTempDF) <- names(MetaDatTemp)
  MetaDatTempDF$SpectraID <- SpectrumID[i]
  MetaDatDFs[[i]] <-  MetaDatTempDF
}
MetaDatNamesDFALL <- ldply(MetaDatDFs, rbind) %>% as.data.frame()


#Make a searchable database with pertinent information-----
Spectra <- cbind(SpectrumID, Spectra, SpectraScore) %>% as.data.frame()
colnames(Spectra) <- c("ID", "spectrum", "score")
Spectra <- Spectra %>% mutate(spectrum_KRHform = spectrum %>% 
                            str_replace_all(" ", "; ") %>%
                            str_replace_all(":", ", "))


#This part gets names - EVERYONE GETS A NAME!
Spectra <- MetaDatNamesDFALL %>%
  group_by(SpectraID) %>%
  summarize(Names = as.character(paste(name,  collapse="; "))) %>%
  rename(ID = SpectraID) %>%
  right_join(., Spectra)

#This part gets exact masses - MOST HAVE EXACT MASSES
Spectra <- MetaDatCMPDDFALL %>%
  filter(name == "total exact mass") %>%
  select(SpectraID, value) %>%
  rename(ID = SpectraID) %>%
  right_join(., Spectra) %>%
  mutate(M_mass = as.numeric(value)) %>%
  select(-value)


#Now we want collision energy - MOST HAVE CEs
Spectra <- MetaDatDFALL %>%
  filter(name == "collision energy") %>%
  select(SpectraID, value) %>%
  rename(ID = SpectraID) %>%
  right_join(., Spectra) %>%
  rename(CE = value) %>%
  select(ID, M_mass, Names, CE, score, spectrum, spectrum_KRHform)


#Rework the spectra so they are in the same format as my spectra-------
library(RCurl)
source_url('https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/XCMS_funs.R')
#library("parallel")
#detectCores()
#mclapply(1:5, function(x) Sys.getpid(), mc.cores = 5)

#Initialize function
reArrangeScans <- function(df){
  RawSpectra <- df["spectrum_KRHform"]
  bestscan <- scantable(RawSpectra) %>%
    arrange(desc(intensity)) %>%
    filter(intensity > 0.5) %>%
    mutate(intensity = round(intensity/max(intensity)*100, digits = 1),
           mz = round (mz, digits = 5),
           mash = paste(mz, intensity, sep = ", " ))
  return(paste(bestscan$mash, collapse = "; "))
}

Spectra$spectrum_KRHform_filtered <- NA
time <- Sys.time()
Spectra$spectrum_KRHform_filtered <- apply(Spectra, 1, FUN=function(x) reArrangeScans(x)) 
print(Sys.time() - time)
setwd("~/Google_Drive/00_XCMS_Working")
write_csv(Spectra, "MoNADat/POS_Spectra.csv")


Spectra <- read_csv("MoNADat/NEG_Spectra.csv")
Spectra$spectrum_KRHform_filtered <- NA
time <- Sys.time()
Spectra$spectrum_KRHform_filtered <- apply(Spectra, 1, FUN=function(x) reArrangeScans(x)) 
print(Sys.time() - time)
setwd("~/Google_Drive/00_XCMS_Working")
write_csv(Spectra, "MoNADat/NEG_Spectra.csv")


##Write out the files ------
#setwd("~/Google_Drive/00_XCMS_Working")

##FOR NEGATIVE MS-MS .json database
#write_csv(Spectra, "MoNADat/NEG_Spectra.csv")
#write_csv(MetaDatCMPDDFALL,  "MoNADat/NEG_CmpInfo.csv")
#write_csv(MetaDatDFALL, "MoNADat/NEG_MetaData.csv")
#write_csv(MetaDatNamesDFALL, "MoNADat/NEG_Names.csv")


##FOR POSITIVE MS-MS .json database
#write_csv(Spectra, "MoNADat/POS_Spectra.csv")

#write_csv(MetaDatCMPDDFALL,  "MoNADat/POS_CmpInfo.csv")
#write_csv(MetaDatDFALL, "MoNADat/POS_MetaData.csv")
#write_csv(MetaDatNamesDFALL, "MoNADat/POS_Names.csv")

