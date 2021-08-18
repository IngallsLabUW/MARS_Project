library(httr)
library(tidyverse)
library(xml2)
source("Functions.R")


### DO NOT USE THIS SCRIPT ###

# Functions ---------------------------------------------------------------
run_number <- 0

getMetlinMz_rmledit <- function(cmpd_mz, ppm = 2.5) {
  if(!is.numeric(cmpd_mz) | cmpd_mz <= 0) {
    na.omit(cmpd_mz)
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
    #select(-Structure) %>%
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

getMetlinMS2 <- function(cmpd_id) { # punched in as the key
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
Experimental.Values <- read.csv("~/work/phobos/example_data/Example_ConfidenceLevel1.csv") %>%
  mutate(mz_experimental = ifelse(is.na(mz_experimental) & !is.na(mz_theoretical), mz_theoretical, mz_experimental)) %>%
  select(MassFeature, compound_theoretical, mz_experimental)

# Scrape Metlin -----------------------------------------------------------
# Mini experimental test runs: TAKES A LONG TIME AND CURRENTLY CAN'T PARALLELLIZE IT FOR SOME REASON

# Scrape for mz
system.time(
  experimental_databymz <- lapply(unique(Experimental.Values[, 3]), getMetlinMz_rmledit)
)
experimental_databymz_df <- bind_rows(experimental_databymz) %>% 
  rename(mz_experimental = Ingalls_cmpd_mz) %>%
  left_join(Experimental.Values, by = "mz_experimental") %>%
  select(MassFeature, everything()) %>%
  unique()

ethyltests <- experimental_databymz_df %>%
  filter(cmpd_name %in% c("Ethyl sulfate", "Ethyl methyl phosphate")) %>%
  select(-CAS, -KEGG, -Structure)


# This is more a scrape for names we already have
Experimental.Values.forname <- Experimental.Values %>%
  filter(!is.na(compound_known))
system.time(
  experimental_databyname <- lapply(unique(Experimental.Values.forname[, 3]), getMetlinName_rmledit)  # basically a standards scrape
) 
experimental_databyname_df <- bind_rows(experimental_databyname)
t <- experimental_databymz_df %>% 
  select(compound_unknown, cmpd_name, formula, KRH_identification, compound_known, mz_unknown, exact_mass) %>% 
  rename(compound_stds = compound_known, mz_experimental = mz_unknown, exact_mass_metlin = exact_mass)
