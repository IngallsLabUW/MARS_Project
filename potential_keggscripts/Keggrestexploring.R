## setup -----
library(tidyverse)
library(stats)
library(cowplot)
library(graphics)
#require(rain)
library(KEGGREST)
require(readr)
#library(lattice)
#library("ChemmineR")
library("CHNOSZ")
library(OrgMassSpecR)

#Exact masses from http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some-------
Hmass=1.00782503223
Cmass=12.0000000 
Omass=15.99491461957 
Nmass=14.00307400443
Pmass=30.97376199842
Smass=31.9720711744
Clmass=34.968852682
Brmass= 78.9183376
Fmass=18.99840316273
Imass=126.904457
Semass=79.916519
Bmass=11.007660
Momass=97.905411
Femass=55.933296
Comass=58.933201
Mgmass=23.985043
Nimass=57.935349
Asmass=74.921600

#Make a list of all compounds in all pathways in KEGG database! This takes a while :) -----
Paths <- names(keggList("pathway"))
badPaths <- c("path:map04080", "path:map04614") #These pathway is full of peptides, pretty problematic
Paths <- Paths[!grepl(paste(badPaths, collapse = "|"), Paths)]

Cmpds <- keggLink("compound", as.character(Paths[1])) 
for (i in 2:length(Paths)){
Cmpds2 <- keggLink("compound", as.character(Paths[i])) 
Cmpds <- c(Cmpds, Cmpds2)
}
Cmpds <- unique(Cmpds) #remove duplicates, leaves us with 5883 compounds
badCmpds <- c("cpd:C00464", "cpd:C01355", "cpd:C00340", "cpd:C00435", "cpd:C00406", "cpd:C01629", "cpd:C19630",
              "cpd:C05524", "cpd:C05996", "cpd:C17936", "cpd:C17937", "cpd:C05677", "cpd:C01667", "cpd:C16663",
              "cpd:C16664", "cpd:C00996", "cpd:C00999", "cpd:C05896", "cpd:C03688", "cpd:C01708", "cpd:C03029",
              "cpd:C03179", "cpd:C05779", "cpd:C05780", "cpd:C05781", "cpd:C05782", "cpd:C05783", "cpd:C05784",
              "cpd:C05785", "cpd:C15804", "cpd:C15805", "cpd:C15806", "cpd:C15807", "cpd:C00338", "cpd:C15922",
              "cpd:C15936", "cpd:C16054", "cpd:C16059", "cpd:C18313", "cpd:C19596", "cpd:C20930", "cpd:C01834",
              "cpd:C06704", "cpd:C06705", "cpd:C06706", "cpd:C06707", "cpd:C21220", "cpd:C21221", "cpd:C21222",
              "cpd:C21241", "cpd:C21242", "cpd:C21243", "cpd:C21245", "cpd:C21223", "cpd:C21224", "cpd:C21225",
              "cpd:C21226", "cpd:C21227", "cpd:C21228", "cpd:C21229", "cpd:C21230", "cpd:C21231", "cpd:C21202",
              "cpd:C21232", "cpd:C21233", "cpd:C21234", "cpd:C21235", "cpd:C21236", "cpd:C21237", "cpd:C21238",
              "cpd:C21239", "cpd:C21244", "cpd:C21246", "cpd:C01647", "cpd:C16003", "cpd:C16004", "cpd:C06865", 
              "cpd:C16005", "cpd:C02017", "cpd:C15906", "cpd:C15949", "cpd:C16010", "cpd:C16012", "cpd:C16013",
              "cpd:C16021", "cpd:C16022", "cpd:C16048", "cpd:C16088", "cpd:C16089", "cpd:C16119", "cpd:C18182",
              "cpd:C00752", "cpd:C01574", "cpd:C01836", "cpd:C02002", "cpd:C02210", "cpd:C02758", "cpd:C15855",
              "cpd:C15856", "cpd:C15869", "cpd:C15890", "cpd:C15891", "cpd:C15901", "cpd:C15995", "cpd:C16016",
              "cpd:C16025", "cpd:C16030", "cpd:C16032", "cpd:C16033", "cpd:C16034", "cpd:C16035", "cpd:C16041",
              "cpd:C16042", "cpd:C16125", "cpd:C16127", "cpd:C00290", "cpd:C16055", "cpd:C16056", "cpd:C16057",
              "cpd:C16058", "cpd:C16060", "cpd:C16061", "cpd:C16062", "cpd:C16114", "cpd:C05606", "cpd:C16051", 
              "cpd:C02049", "cpd:C08822", "cpd:C05005", "cpd:C01074", "cpd:C01498", "cpd:C18287", "cpd:C18288",
              "cpd:C18289", "cpd:C18290", "cpd:C18291", "cpd:C18292", "cpd:C18991", "cpd:C18997", "cpd:C18998",
              "cpd:C18994", "cpd:C18995", "cpd:C18996", "cpd:C19313") #this compounds don't have formulas, so I'm cutting them before getting info.
Cmpds <- Cmpds[!grepl(paste(badCmpds, collapse = "|"), Cmpds)]


KEGGCompounds <- data.frame(Cmpds, Cmpds) %>%
  select(Cmpds) %>%
  mutate(Name = NA) %>%
  mutate(Formula = NA) %>%
  mutate(AllNames = NA) %>%
  mutate(AllReactions = NA) %>%
  mutate(AllPath = NA) %>%
  mutate(AllEnzyme = NA) %>%
  mutate(AllModule = NA) %>%
  mutate(AllEnzyme = NA) %>%
  mutate(AllBrite = NA) %>%
  mutate(AllDB = NA) %>%
  mutate(OtherCmpds = NA)


for (i in 2:length(Cmpds)){
  CmpInfo <- keggGet(Cmpds[i])
  KEGGCompounds[i, "Name"] <- CmpInfo[[1]]$NAME[1]
  KEGGCompounds[i, "Formula"] <- CmpInfo[[1]]$FORMULA
  KEGGCompounds[i, "AllNames"] <- paste(CmpInfo[[1]]$NAME, sep="", collapse=" ")
  KEGGCompounds[i, "AllReactions"] <- paste(CmpInfo[[1]]$REACTION, sep="", collapse="; ")
  KEGGCompounds[i, "AllPath"] <- paste(CmpInfo[[1]]$PATHWAY, sep="", collapse="; ")
  KEGGCompounds[i, "AllEnzyme"] <- paste(CmpInfo[[1]]$ENZYME, sep="", collapse="; ")
  KEGGCompounds[i, "AllModule"] <- paste(CmpInfo[[1]]$MODULE, sep="", collapse="; ")
  KEGGCompounds[i, "AllDB"] <- paste(CmpInfo[[1]]$DBLINKS, sep="", collapse="; ")
  KEGGCompounds[i, "AllBrite"] <- paste(CmpInfo[[1]]$BRITE, sep="", collapse="; ")
  KEGGCompounds[i, "OtherCmpds"] <- paste(names(keggFind("compound",
                                                         as.character(KEGGCompounds[i, "Formula"]),
                                                         "formula")),  sep="", collapse="; ")
}

#setwd("/Users/katherine2/Google_Drive/00_XCMS_Working/KEGG_Explore")
#write.csv(KEGGCompounds, "allKEGGCompounds.csv") #SAVE INTERMEDIATE WITH ALL THIS JUNK IN IT!!

#Calculating exact masses for these compounds in Pos, Neg, Neutral-------

#KEGGCompounds <- read_csv("~/Google_Drive/00_XCMS_Working/KEGG_Explore/allKEGGCompounds.csv")

#Get rid of compounds that we will not be able to predict a mass for
badforms <- c("R", "n", ")", ". ")
KEGGCompounds2 <- KEGGCompounds[!grepl(paste(badforms, collapse = "|"), KEGGCompounds$Formula) ,]
KEGGCompounds2 <- KEGGCompounds2[grepl("C", KEGGCompounds2$Formula) ,]
KEGGCompounds2 <- KEGGCompounds2[grepl("H", KEGGCompounds2$Formula) ,]

#Get formulas split to calculate +H or -H masses
Formulas <- as.data.frame(t(count.elements( KEGGCompounds2[1, "Formula"])))
Formulas$Cmpds <- KEGGCompounds2[1, "Cmpds"]
Cmpds <- KEGGCompounds2$Cmpds
for (i in 2:length(Cmpds)){
  thisFormulaDF <- as.data.frame(t(count.elements(KEGGCompounds2[i, "Formula"])))
  thisFormulaDF$Cmpds <- KEGGCompounds2[i, "Cmpds"]
  Formulas <-bind_rows(Formulas, thisFormulaDF)
}
Formulas$Cmpds <- as.character(Formulas$Cmpds) 
Formulas2 <- Formulas %>%
  mutate(PosForm = paste("C", C, "H", (H+1), "O", O, "N", N,
                         "P", P, "S", S, "Cl", Cl, "I", I, "Br", Br,   
                         "F", F, "Se", Se, "B", B, "Mo", Mo, "Fe", Fe,
                         "Co", Co, "Mg", Mg, "Ni", Ni, "As", As,
                         sep = "", collapse = NULL)) %>%
  mutate(NegForm = paste("C", C, "H", (H-1), "O", O, "N", N,
                         "P", P, "S", S, "Cl", Cl, "I", I, "Br", Br,   
                         "F", F, "Se", Se, "B", B, "Mo", Mo, "Fe", Fe,
                         "Co", Co, "Mg", Mg, "Ni", Ni, "As", As,
                         sep = "", collapse = NULL)) %>%
  mutate(PosForm = PosForm %>% 
           str_replace("CNA","") %>%
           str_replace("HNA","") %>%
           str_replace("ONA","") %>%
           str_replace("NNA","") %>%
           str_replace("PNA","") %>%
           str_replace("SNA","") %>%
           str_replace("ClNA","") %>%
           str_replace("INA","") %>%
           str_replace("BrNA","") %>%
           str_replace("FNA","") %>%
           str_replace("SeNA","") %>%
           str_replace("BNA","") %>%
           str_replace("MoNA","") %>%
           str_replace("FeNA","") %>%
           str_replace("CoNA","") %>%
           str_replace("MgNA","") %>%
           str_replace("NiNA","") %>%
           str_replace("AsNA","")) %>%
  mutate(NegForm = NegForm %>% 
           str_replace("CNA","") %>%
           str_replace("HNA","") %>%
           str_replace("ONA","") %>%
           str_replace("NNA","") %>%
           str_replace("PNA","") %>%
           str_replace("SNA","") %>%
           str_replace("ClNA","") %>%
           str_replace("INA","") %>%
           str_replace("BrNA","") %>%
           str_replace("FNA","") %>%
           str_replace("SeNA","") %>%
           str_replace("BNA","") %>%
           str_replace("MoNA","") %>%
           str_replace("FeNA","") %>%
           str_replace("CoNA","") %>%
           str_replace("MgNA","") %>%
           str_replace("NiNA","") %>%
           str_replace("AsNA","")) %>%
  mutate(NeutalMass = NA )%>%
  mutate(PosMZ = NA )%>%
  mutate(NegMZ = NA )
Formulas2[is.na(Formulas2)] <- 0

for (i in 2:length(Formulas2$Cmpds)){
  Formulas2$NeutalMass[i] <- (Hmass*Formulas2$H[i] + Cmass*Formulas2$C[i]+ Omass*Formulas2$O[i]+ 
                                Nmass*Formulas2$N[i]+ Pmass*Formulas2$P[i]+  Smass*Formulas2$S[i]+ 
                                Clmass*Formulas2$Cl[i]+ Brmass* Formulas2$Br[i]+ Fmass*Formulas2$F[i]+ 
                                Imass*Formulas2$I[i]+  Semass*Formulas2$Se[i]+ Bmass*Formulas2$B[i]+ 
                                Momass*Formulas2$Mo[i]+ Femass*Formulas2$Fe[i]+ Comass*Formulas2$Co[i]+ 
                                Mgmass*Formulas2$Mg[i]+ Nimass*Formulas2$Ni[i]+ Asmass*Formulas2$As[i])
  Formulas2$PosMZ[i] <- (Hmass*(Formulas2$H[i]+1) + Cmass*Formulas2$C[i]+ Omass*Formulas2$O[i]+ 
                                Nmass*Formulas2$N[i]+ Pmass*Formulas2$P[i]+  Smass*Formulas2$S[i]+ 
                                Clmass*Formulas2$Cl[i]+ Brmass* Formulas2$Br[i]+ Fmass*Formulas2$F[i]+ 
                                Imass*Formulas2$I[i]+  Semass*Formulas2$Se[i]+ Bmass*Formulas2$B[i]+ 
                                Momass*Formulas2$Mo[i]+ Femass*Formulas2$Fe[i]+ Comass*Formulas2$Co[i]+ 
                                Mgmass*Formulas2$Mg[i]+ Nimass*Formulas2$Ni[i]+ Asmass*Formulas2$As[i])
  Formulas2$NegMZ[i] <- (Hmass*(Formulas2$H[i]-1) + Cmass*Formulas2$C[i]+ Omass*Formulas2$O[i]+ 
                                Nmass*Formulas2$N[i]+ Pmass*Formulas2$P[i]+  Smass*Formulas2$S[i]+ 
                                Clmass*Formulas2$Cl[i]+ Brmass* Formulas2$Br[i]+ Fmass*Formulas2$F[i]+ 
                                Imass*Formulas2$I[i]+  Semass*Formulas2$Se[i]+ Bmass*Formulas2$B[i]+ 
                                Momass*Formulas2$Mo[i]+ Femass*Formulas2$Fe[i]+ Comass*Formulas2$Co[i]+ 
                                Mgmass*Formulas2$Mg[i]+ Nimass*Formulas2$Ni[i]+ Asmass*Formulas2$As[i])
  }
    
KEGGCompounds2 <- full_join(KEGGCompounds2, Formulas2, by = "Cmpds") 
#write.csv(KEGGCompounds2, "KEGGCompounds_withMasses.csv")

#Example:

#Search for "Cysteine and methionine metabolism" in Allpaths of big dataset.

