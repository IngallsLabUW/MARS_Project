library(tidyverse)
source("Functions.R")

## Some notes:
# Start with scraping the KRH MFCluster compounds. Metlin should be lower priority than MoNA, maybe after 
# scraping MoNA return some stats and ask, do you want to try metlin? Then follow up with "this compound
# wasn't found as of the mm/dd/yy scrape", etc.

# Make a baseline dataframe of results from each successive scrape for now.
# Later we will try a larger scale scrape to see how it goes.

auth_url <- "https://metlin.scripps.edu/lib/json/user.php"
resp <- POST(auth_url, body = list(
  user="wkumler@uw.edu",
  password="password",
  action="login"
), add_headers(referer="https://metlin.scripps.edu/landing_page.php?pgcontent=mainPage"))
content(resp)


run_number <- 0

getMetlinName_rmledit <- function(name) {
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
              length(unique(search_data$formula)), "unique formula(s)."))
  
  if (nrow(MSMS_subset)) {
    print(paste("Of those,", nrow(MSMS_subset),
                "have experimental MS/MS data", collapse = ", "))
  } else {
    print("None of these compounds have experimental MS/MS data")
  }
  run_number <<- run_number + 1
  print(paste("Run number: ", run_number))

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
  select("m.z", contains("Name")) %>%
  slice(1:5)


## Mini test run

t2 <- lapply(unique(Standards.to.Scrape[, 3]), getMetlinName_rmledit)
df <- bind_rows(t2, .id = "column_label")
#write.csv(df, "data_underway/first_metlin_scrape.csv")
####

# argo_data_by_name <- getMetlinName("Argininosuccinate")
# argo_data_by_mz <- getMetlinMz(291.13046 - 1.007276)
# argo_MSMS_data <- argo_data_by_name %>%
#   filter(cmpd_name == "Argininosuccinic acid") %>%
#   filter(MSMS == "experimental") %>%
#   pull(cmpd_id) %>%
#   as.numeric() %>%
#   getMetlinMS2()
# 
# aspartic_data_by_name <- getMetlinName("Aspartic acid")

####
