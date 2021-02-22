library(httr)
library(tidyverse)
library(xml2)
source("Functions.R")

## Some notes:k
# Would love to use a closure for that run_number variable
# Mz function needs to drop nas, for now I'm doing it manually
# Once the output dataframes are settled, these scrapes will form the base of the
# google drive metlin dataframe.
# Need to make the scraped database for this: decide on a format and make the first one.

auth_url <- "https://metlin.scripps.edu/lib/json/user.php"
resp <- POST(auth_url, body = list(
  user="wkumler@uw.edu",
  password="password",
  action="login"
), add_headers(referer="https://metlin.scripps.edu/landing_page.php?pgcontent=mainPage"))
content(resp)


# Functions ---------------------------------------------------------------
run_number <- 0

getMetlinName_rmledit <- function(mf.name) {
  metlin.data <- paste0("https://metlin.scripps.edu/advanced_search_result.php?",
                        "name=", mf.name, "&AminoAcid=add&drug=add&toxinEPA=add&",
                        "keggIDFilter=add") %>%
    GET(add_headers(referer = "https://metlin.scripps.edu/advanced_search.php")) %>%
    content(encoding = "UTF-8") %>%
    xml_find_all(xpath = "//tbody")

  search.ids <- metlin.data %>%
    xml_find_all(xpath = "//th[@scope]/a") %>%
    xml_text()

  search.data <- metlin.data %>%
    xml_find_all(xpath="//td") %>%
    xml_text() %>%
    matrix(ncol=7, byrow=T) %>%
    cbind(search.ids, .) %>%
    as.data.frame() %>%
    `names<-`(c("cmpd_id", "exact_mass", "cmpd_name", "formula", #### exact mass already accounts for proton?
                "CAS", "KEGG", "MSMS", "Structure")) %>%
    select(-Structure) %>%
    mutate(Ingalls_cmpd_name = mf.name) 

  MSMS.subset <- subset(search.data, MSMS == "experimental")

  print(paste("Metlin returned", nrow(search.data), "compound(s) with name",
              mf.name, "with",
              length(unique(search.data$formula)), "unique formula(s)."))
  
  if (nrow(MSMS.subset)) {
    print(paste("Of those,", nrow(MSMS.subset),
                "have experimental MS/MS data", collapse = ", "))
  } else {
    print("None of these compounds have experimental MS/MS data")
  }
  run_number <<- run_number + 1
  print(paste("Run number: ", run_number))

  return(search.data)
}

getMetlinrt_rmledit <- function(rt, ppm = 2.5) {
  if(!is.numeric(cmpd.mz) | cmpd.mz <= 0) {
    stop("Mass must be positive and numeric")
  }
  if(!is.numeric(ppm) | ppm <= 0) {
    stop("Mass must be positive and numeric")
  }
  
  
  mz.range <- c(cmpd.mz - cmpd.mz * ppm / 1000000, cmpd.mz + cmpd.mz * ppm / 1000000)
  metlin.data <- paste0("https://metlin.scripps.edu/advanced_search_result.php?",
                        "molid=&mass_min=", min(mz.range),
                        "&mass_max=", max(mz.range), "&Amino",
                        "Acid=add&drug=add&toxinEPA=add&keggIDFilter=add") %>%
    GET() %>%
    stop_for_status() %>%
    content(encoding = "UTF-8") %>%
    xml_find_all(xpath = "//tbody")
  
  search.ids <- metlin.data %>%
    xml_find_all(xpath = "//th[@scope]/a") %>%
    xml_text()
  
  search.data <- metlin.data %>%
    xml_find_all(xpath="//td") %>%
    xml_text() %>%
    matrix(ncol=7, byrow=T) %>%
    cbind(search.ids, .) %>%
    as.data.frame() %>%
    `names<-`(c("cmpd_id", "exact_mass", "cmpd_name", "formula",
                "CAS", "KEGG", "MSMS", "Structure")) %>%
    select(-Structure) %>%
    mutate(Ingalls_cmpd.mz = cmpd.mz)
  
  MSMS.subset <- subset(search.data, MSMS=="experimental")
  
  print(paste("Metlin returned", nrow(search.data), "compound(s) between",
              min(mz.range), "and", max(mz.range), "m/z with",
              length(levels(search.data$formula)), "unique formula(s):",
              paste(levels(search.data$formula), collapse = ", ")))
  if(nrow(MSMS.subset)){
    print(paste("Of those,", nrow(MSMS.subset),
                "have experimental MS/MS data:",
                paste(as.character(MSMS.subset$cmpd_name), collapse = ", ")))
  } else {
    print("None of these compounds have experimental MS/MS data")
  }
  
  return(search.data)
}

getMetlinMS2 <- function(cmpd_id){ # punched in as the key
  if(!is.numeric(cmpd_id)|cmpd_id<=0){
    stop("Mass must be positive and numeric")
  }
  metlin.data <- paste0("https://metlin.scripps.edu/showChart.php?molid=",
                        cmpd_id,"&etype=experimental") %>%
    GET() %>% stop_for_status()
  ms2_raw <- metlin.data %>%
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


# Data --------------------------------------------------------------------
Experimental.Values <- read.csv("data_processed/confidence_level1.csv") %>%
  select(-X) %>%
  filter(!is.na(compound_unknown)) %>%
  slice(1:66)

# Scrape Metlin -----------------------------------------------------------
# Mini experimental testrun: TAKES A LONG TIME AND CURRENTLY CAN'T PARALLELLIZE IT FOR SOME REASON
system.time(
  experimental_databyname <- lapply(unique(Experimental.Values[, 3]), getMetlinName_rmledit)  # basically a standards scrape
  ) 
experimental_databyname_df <- bind_rows(experimental_databyname)


experimental_databymz <- lapply(unique(Experimental.Values[, 4]), getMetlinMz_rmledit)
experimental_databymz_df <- bind_rows(experimental_databymz) %>% unique() %>%
  rename(mz_unknown = Ingalls_cmpd.mz)

t <- experimental_databymz_df %>% left_join(Experimental.Values, by = "mz_unknown")
#write.csv(df, "data_underway/first_metlin_scrape.csv")

