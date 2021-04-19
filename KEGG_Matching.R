library(fuzzyjoin)
library(janitor)
library(tidyverse)

# Require ppm tolerance, if user doesn't enter it then use 15
# Figure out z/hilic stuff. What about reverse phase/cyano?
# fuzzyjoin to KEGG list
# combine kegg names to kegg ids


## Need to figure out how to download KEGGwithMasses
KEGG.Data <- read.csv("data_extra/KEGGCompounds_withMasses.csv", header = TRUE)

Experimental.Data <- read.csv("data_processed/Example_ConfidenceLevel2.csv") %>%
  rename(mz = mz_experimental) %>%
  select(MassFeature, compound_experimental, mz, column_experimental, z_experimental) %>%
  unique()


#matchedKEGGs <- list() # Why is this here
keggCompounds <- KEGG.Data %>%
  rename(Compound = OtherCmpds)

keggPos <- keggCompounds %>% 
  select(Compound, PosMZ) %>% 
  mutate(column = "HILIC",
         z = 1) %>%
  # mutate(Fraction1 = "HILICPos", Fraction2="CyanoAq", Fraction3= "CyanoDCM") %>% unique() %>% Why is this here too
  rename(mz = PosMZ)

keggNeg <- keggCompounds %>% 
  select(Compound, NegMZ) %>% 
  mutate(column = "HILIC",
         z = -1) %>%
  # mutate(Fraction1 = "HILICPos", Fraction2="CyanoAq", Fraction3= "CyanoDCM") %>% unique() %>% Why is this here too
  rename(mz = NegMZ)

keggCompoundsShort <- keggNeg %>%
  rbind(keggPos)

#MFstry <- Experimental.Data %>% mutate(MF_Frac2 = MF_Frac) %>% separate(MF_Frac2, c("MF", "Frac"), sep =  "_X_") %>% select(-MF) 
matched <- difference_inner_join(x = keggCompoundsShort, y = Experimental.Data, # originally mfstry
                                 by = "mz", max_dist = .02, distance_col = NULL) %>%
  rename(mz_KEGG = mz.x,
         mz_experimental = mz.y) %>%
  filter(z == z_experimental) %>%
  unique() %>% # can maybe drop 
  mutate(ppm = (abs(mz_KEGG - mz_experimental) / mz_KEGG *10^6)) %>%
  filter(ppm < 15) %>% #ppmtol, #Frac == Fraction1 | Frac == Fraction2 | Frac == Fraction3) %>%
  mutate(mz_similarity_scoreKEGG = exp(-0.5 * (((mz_experimental - mz_KEGG) / 0.02) ^ 2))) %>%
  #mutate(Frac = as.factor(Frac), mz = mz.x) %>% 
  #select(MF_Frac, Compound, ppm ) %>% 
  rename(Keggppm = ppm)


matchedNames <- left_join(matched, keggCompounds) %>% 
  select(Compound, Name) %>%
  group_by(Compound) %>%
  summarise(KEGGMatchesNames = as.character(paste(Name,  collapse=" ")))

# matchedKEGGs[[1]] <- left_join(matched, matchedNames) %>% 
#   rename(KeggMatches = Compound)
final <- matched %>%
  left_join(matchedNames) %>%
  arrange(Keggppm) %>%
  group_by(MassFeature) 




## Will's code
#########################################################3
library(tidyverse)
library(httr)
library(xml2)
library(pbapply)
# Biosynthesis of amino acids
base_url <- "https://www.genome.jp/kegg-bin/show_pathway?map01230"
base_data <- GET(base_url) %>% content()
mol_ids <- base_data %>%
  xml_find_all("//map[@id='mapdata']/area[@shape='circle']") %>%
  xml_attr("data-entry")
pathway_cmpds <- pbsapply(mol_ids, function(mol_id){
  cmpd_data <- "https://www.genome.jp/dbget-bin/www_bget?cpd:" %>%
    paste0(mol_id) %>%
    GET() %>%
    stop_for_status() %>%
    content()
  cmpd_mz <- xml_find_first(cmpd_data, "//nobr[text() = 'Exact mass']/../../td/div") %>%
    xml_text() %>% as.numeric()
  if(is.na(cmpd_mz)){return(NULL)}
  cmpd_name <- xml_find_first(cmpd_data, "//nobr[text() = 'Name']/../../td/div") %>%
    xml_text() %>% strsplit(";") %>% unlist() %>% `[`(1)
  data.frame(mol_id, cmpd_name, cmpd_mz)
}) %>%
  do.call(what = "rbind")
##############################################################
