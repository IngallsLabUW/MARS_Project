library(fuzzyjoin)
library(tidyverse)
options(scipen = 999)

# Require ppm tolerance, if user doesn't enter it then use 15
# Figure out z/hilic stuff. What about reverse phase/cyano?
# Need to figure out how to download KEGGwithMasses

# Define experimental and theoretical data
Experimental.Data <- read.csv("data_processed/Example_ConfidenceLevel2.csv") 

KEGG.Data <- read.csv("data_extra/KEGGCompounds_withMasses.csv", header = TRUE) %>%
  rename(Compound_KEGG = OtherCmpds)

ToJoin <- Experimental.Data %>%
  rename(mz = mz_experimental) %>%
  select(MassFeature, compound_experimental, mz, column_experimental, z_experimental) %>%
  unique()

# Isolate pos and neg theoretical data
keggPos <- KEGG.Data %>% 
  select(Compound_KEGG, PosMZ) %>% 
  mutate(column_KEGG = "HILIC",
         z = 1) %>%
  rename(mz = PosMZ)

keggNeg <- KEGG.Data %>% 
  select(Compound_KEGG, NegMZ) %>% 
  mutate(column_KEGG = "HILIC",
         z = -1) %>%
  rename(mz = NegMZ)

# Combine to full theoretical dataset
keggCompounds <- keggNeg %>%
  rbind(keggPos)

# Fuzzyjoin datasets
My.Fuzzy.Join <- difference_inner_join(x = keggCompounds, y = ToJoin, 
                                 by = "mz", max_dist = .02, distance_col = NULL) %>%
  rename(mz_KEGG = mz.x,
         mz_experimental = mz.y,
         z_KEGG = z) %>%
  mutate(ppm = (abs(mz_KEGG - mz_experimental) / mz_KEGG *10^6), 
         mz_similarity_scoreKEGG = exp(-0.5 * (((mz_experimental - mz_KEGG) / 0.02) ^ 2))) %>%
  filter(ppm < 15, 
         z_KEGG == z_experimental) %>%
  unique() %>% 
  rename(KEGGppm = ppm) %>%
  arrange(compound_experimental)

# Match KEGG names to IDs
matchedNames <- left_join(My.Fuzzy.Join, KEGG.Data) %>% 
  select(Compound, Name) %>%
  group_by(Compound) %>%
  summarise(KEGGMatchesNames = as.character(paste(Name,  collapse=" ")))

# Combine to full matched df
final <- My.Fuzzy.Join %>%
  left_join(matchedNames) %>%
  select(MassFeature, compound_experimental, mz_KEGG, z_KEGG, Compound, KEGGMatchesNames, KEGGppm, mz_similarity_scoreKEGG)


t <- final %>%
  difference_left_join(Experimental.Data, by = c("compound_experimental"))




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
