library(httr)
library(tidyverse)
library(xml2)
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
    `names<-`(c("cmpd_id", "exact_mass", "cmpd_name", "formula", #### exact mass already accounts for proton?
                "CAS", "KEGG", "MSMS", "Structure")) %>%
    select(-Structure) %>%
    mutate(Ingalls_cmpd_name = name)

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

getMetlinMz_rmledit <- function(cmpd_mz, ppm = 2.5) {
  if(!is.numeric(cmpd_mz) | cmpd_mz <= 0) {
    stop("Mass must be positive and numeric")
  }
  if(!is.numeric(ppm) | ppm <= 0) {
    stop("Mass must be positive and numeric")
  }
  
  
  mz_range <- c(cmpd_mz-cmpd_mz*ppm/1000000, cmpd_mz+cmpd_mz*ppm/1000000)
  metlin_data <- paste0("https://metlin.scripps.edu/advanced_search_result.php?",
                        "molid=&mass_min=", min(mz_range),
                        "&mass_max=", max(mz_range), "&Amino",
                        "Acid=add&drug=add&toxinEPA=add&keggIDFilter=add") %>%
    GET() %>%
    stop_for_status() %>%
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
    select(-Structure) %>%
    mutate(Ingalls_cmpd_mz = cmpd_mz)
  
  MSMS_subset <- subset(search_data, MSMS=="experimental")
  
  print(paste("Metlin returned", nrow(search_data), "compound(s) between",
              min(mz_range), "and", max(mz_range), "m/z with",
              length(levels(search_data$formula)), "unique formula(s):",
              paste(levels(search_data$formula), collapse = ", ")))
  if(nrow(MSMS_subset)){
    print(paste("Of those,", nrow(MSMS_subset),
                "have experimental MS/MS data:",
                paste(as.character(MSMS_subset$cmpd_name), collapse = ", ")))
  } else {
    print("None of these compounds have experimental MS/MS data")
  }
  
  return(search_data)
}

getMetlinMS2 <- function(cmpd_id){ # punched in as the key
  if(!is.numeric(cmpd_id)|cmpd_id<=0){
    stop("Mass must be positive and numeric")
  }
  metlin_data <- paste0("https://metlin.scripps.edu/showChart.php?molid=",
                        cmpd_id,"&etype=experimental") %>%
    GET() %>% stop_for_status()
  ms2_raw <- metlin_data %>%
    content(encoding = "UTF-8") %>%
    xml_find_all(xpath = '/html/head/script[5]') %>%
    xml_text() %>%
    gsub(pattern = ".*series: \\[", replacement = "") %>%
    gsub(pattern = " ]}]\n                }); \n\n            });\n        ",
         replacement = "") %>%
    gsub(pattern = "&nbsp;", replacement = "") %>%
    strsplit(split = "\\{name: ") %>%
    unlist() %>%
    `[`(-1) %>%
    lapply(strsplit, split=",data:\\[")
  
  if(!length(ms2_raw)){stop("No MS2 data found. Are you sure Metlin has this data?")}
  
  ms2_titles <- sapply(ms2_raw, `[[`, 1)[1,] #Not sure why this works but OK
  ms2_polarities <- ifelse(grepl(ms2_titles, pattern = "(\\+)"), "+", "-")
  ms2_voltages <- gregexpr("\\d+\\.?\\d*", ms2_titles) %>%
    regmatches(x=ms2_titles) %>%
    unlist() %>% as.numeric()
  ms2_adducts <- gregexpr("\\[M.*\\]", ms2_titles) %>%
    regmatches(x = ms2_titles) %>%
    unlist()
  
  ms2_list <- sapply(ms2_raw, `[[`, 1)[2,] %>%
    gsub(pattern = " ]},", replacement = "") %>%
    lapply(function(x){
      wo1 <- substring(x, 2)
      wo2 <- substring(wo1, 1, nchar(wo1))
      ms2_list <- unlist(strsplit(wo2, "},\\{"))
      ms2_list <- lapply(ms2_list, function(x)
        as.numeric(unlist(regmatches(x, m = gregexpr("\\d+\\.?\\d*", x)))))
      ms2_df <- as.data.frame(do.call(rbind, ms2_list))
      names(ms2_df) <- c("mass", "rel_intensity")
      return(ms2_df)
    })
  for(i in seq_along(ms2_list)){
    ms2_list[[i]] <- cbind(ms2_polarities[i], ms2_adducts[i],
                           ms2_voltages[i], ms2_list[[i]])
  }
  ms2_df <- as.data.frame(do.call(rbind, ms2_list))
  names(ms2_df) <- c("polarity", "adduct", "voltage", "frag_mass", "frag_int")
  
  print(paste("Metlin had", length(ms2_titles), "MS2 records for this compound,",
              "with collision energies of",
              paste(paste0(ms2_polarities, ms2_voltages), collapse=", ")))
  
  return(ms2_df)
}
#######

Experimental.Values <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") %>%
  separate(MassFeature_Column, into = c("MassFeature", "column"), sep = "_") %>%
  select(-contains("cluster"), -column) %>%
  slice(1:5)

Standards.to.Scrape <- read.csv("data_extra/Ingalls_Lab_Standards_Jan11.csv", stringsAsFactors = FALSE) %>%
  select("m.z", contains("Name")) %>%
  transform(m.z = as.numeric(m.z)) %>%
  slice(1:3, 24)


## Mini test run, some compounds are missing.

test_databyname <- lapply(unique(Standards.to.Scrape[, 3]), getMetlinName_rmledit)
test_databyname_df <- bind_rows(test_databyname, .id = "column_label")

test_databymz <- lapply(unique(Standards.to.Scrape[, 1]), getMetlinMz_rmledit)
test_databymz_df <- bind_rows(test_databymz)


experimental_databyname <- lapply(unique(Experimental.Values[, 1]), getMetlinName_rmledit)
experimental_databyname_df <- bind_rows(experimental_databyname)

experimental_databymz <- lapply(unique(Experimental.Values[, 4]), getMetlinMz_rmledit)
experimental_databymz_df <- bind_rows(experimental_databymz)

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
