## Will's code
#########################################################
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

path <- "http://rest.kegg.jp/list/compound"

r <- GET(url = path) %>%
  content()

t <- read.csv(text = gsub("\t\n", "", r), sep = "\t", header = FALSE, col.names = c("cmpd", "name","hmm"))

masspath <- "http://rest.kegg.jp/find/compound/exact_mass"	
emass <- GET(url = masspath) %>%
  content()
tmass <- read.csv(text = gsub("\t\n", "", emass), sep = "\t", header = FALSE, col.names = c("cmpd", "exact_mass"))
##############################################################
