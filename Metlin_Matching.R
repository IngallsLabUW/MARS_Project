library(tidyverse)
source("Functions.R")
source("KRH_Functions.R")

## Some notes:
# Start with scraping the KRH MFCluster compounds. Metlin should be lower priority than MoNA, maybe after 
# scraping MoNA return some stats and ask, do you want to try metlin? Then follow up with "this compound
# wasn't found as of the mm/dd/yy scrape", etc.

# Make a baseline dataframe of results from each successive scrape for now.
# Later we will try a larger scale scrape to see how it goes.

file.pattern <- "bet_"

filenames <- RemoveCsv(list.files(path = "data_from_lab_members", pattern = file.pattern))
filepath <- file.path("data_from_lab_members", paste(filenames, ".csv", sep = ""))

for (i in filenames) {
  filepath <- file.path("data_from_lab_members", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}