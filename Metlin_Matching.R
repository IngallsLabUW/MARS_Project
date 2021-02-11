library(tidyverse)
source("Functions.R")
source("getMetlinFunctions.R")

## Some notes:
# Start with scraping the KRH MFCluster compounds. Metlin should be lower priority than MoNA, maybe after 
# scraping MoNA return some stats and ask, do you want to try metlin? Then follow up with "this compound
# wasn't found as of the mm/dd/yy scrape", etc.

# Make a baseline dataframe of results from each successive scrape for now.
# Later we will try a larger scale scrape to see how it goes.

getMetlinName <- function(name) {
  metlin_data <- paste0("https://metlin.scripps.edu/advanced_search_result.php?",
                        "name=", name, "&AminoAcid=add&drug=add&toxinEPA=add&",
                        "keggIDFilter=add") %>%
    GET(add_headers(referer = "https://metlin.scripps.edu/advanced_search.php")) %>%
    #stop_for_status() %>% ## ISSUE IS HERE
    content(encoding = "UTF-8") %>%
    xml_find_all(xpath = "//tbody")

  search_ids <- metlin_data %>%
    xml_find_all(xpath = "//th[@scope]/a") %>%
    xml_text()

  search_data <- metlin_data %>%
    xml_find_all(xpath="//td") %>%
    xml_text() %>%
    matrix(ncol=7, byrow=T) %>%
    cbind(search_ids, .) %>%
    as.data.frame() %>%
    `names<-`(c("cmpd_id", "exact_mass", "cmpd_name", "formula",
                "CAS", "KEGG", "MSMS", "Structure")) %>%
    select(-Structure)

  MSMS_subset <- subset(search_data, MSMS=="experimental")

  print(paste("Metlin returned", nrow(search_data), "compound(s) with name",
              name, "with",
              length(unique(search_data$formula)), "unique formula(s):",
              paste(search_data$formula, collapse = ", ")))
  if(nrow(MSMS_subset)){
    print(paste("Of those,", nrow(MSMS_subset),
                "have experimental MS/MS data:",
                paste(as.character(MSMS_subset$cmpd_name), collapse = ", ")))
  } else {
    print("None of these compounds have experimental MS/MS data")
  }

  return(search_data)
}

#######

file.pattern <- "bet_"

filenames <- RemoveCsv(list.files(path = "data_extra/Metlin_RelationalSpreadsheets/", pattern = file.pattern))
filepath <- file.path("data_extra/Metlin_RelationalSpreadsheets/", paste(filenames, ".csv", sep = ""))

for (i in filenames) {
  filepath <- file.path("data_extra/Metlin_RelationalSpreadsheets/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}

Experimental.Values <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") %>%
  select(-contains("cluster"))

Standards.to.Scrape <- read.csv("data_extra/Ingalls_Lab_Standards_Jan11.csv", stringsAsFactors = FALSE) %>%
  select("m.z", contains("Name"))

####

argo_data_by_name <- getMetlinName("Argininosuccinate")
argo_data_by_mz <- getMetlinMz(291.13046 - 1.007276)
argo_MSMS_data <- argo_data_by_name %>%
  filter(cmpd_name == "Argininosuccinic acid") %>%
  filter(MSMS == "experimental") %>%
  pull(cmpd_id) %>%
  as.numeric() %>%
  getMetlinMS2()

aspartic_data_by_name <- getMetlinName("Aspartic acid")

####

## Mini test run

Test.Standards <- Standards.to.Scrape %>%
  slice(1:10) 

t <- getMetlinName(Test.Standards[1, 1])
t2 <- lapply(unique(Test.Standards[, 3]), getMetlinName)


for (i in Test.Standards) {
  print(i)
  getMetlinName(i)
}