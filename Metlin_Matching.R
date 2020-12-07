library(tidyverse)
source("Functions.R")
source("KRH_Functions.R")


file.pattern <- "bet_"

filenames <- RemoveCsv(list.files(path = "data_from_lab_members", pattern = file.pattern))
filepath <- file.path("data_from_lab_members", paste(filenames, ".csv", sep = ""))

for (i in filenames) {
  filepath <- file.path("data_from_lab_members", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}