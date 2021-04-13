library(splitstackshape)
library(tidyverse)

Cyano.MS2 <- read.csv("data_extra/Standards_MS2/Cyano_stds_withMS2s.csv") %>%
  mutate(column = "RP",
         z = NA)

NonScaled.MS2 <- read.csv("~/work/Ingalls_Standards/MSMS/data_processed/Ingalls_Lab_Standards_MSMS.csv")

# Functions ---------------------------------------------------------------
ConcatToScan <- function(concat.format) {
  scantable <- cSplit(concat.format, "MS2", sep = ";") %>%
    pivot_longer(cols = starts_with("MS2")) %>%
    separate(value, sep = ",", into = c("mz", "intensity")) %>%
    select(compound, mz, intensity, file_name) #%>%
    #drop_na()
  
  return(scantable)
}

ScanToConcat <- function(scantable.format) {
  concatenated <- scantable.format %>%
    unique() %>%
    group_by(compound, file_name) %>%
    unite("MS2", mz:intensity, sep = ", ") %>%
    ungroup() %>%
    group_by(compound, file_name) %>%
    mutate(MS2 = paste(MS2, collapse = "; ")) %>%
    unique() %>%
    as.data.frame() %>%
    mutate(MS2 = ifelse(MS2 == ("NA, NA"), NA, MS2))
  
  return(concatenated)
}

ScaleMS2 <- function(scan) {
  scantable <- read.table(text = as.character(scan),
                          col.names = c("mz", "intensity"), fill = TRUE) %>%
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", "")),
           intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
    filter(intensity > 0.5) %>%
    arrange(desc(intensity)) %>%
    unite("MS2", mz:intensity, sep = ", ") %>%
    mutate(MS2 = paste(MS2, collapse = "; ")) %>%
    unique()
  
  return(scantable)
}

  
# Both Layouts ------------------------------------------------------------
Concatenated.Format <- Cyano.MS2 %>%
  rename(MS2 = MS2s) %>%
  select(compound, MS2, file_name)

Scantable.Format <- cSplit(Concatenated.Format, "MS2", sep = ";") %>%
  group_by(compound, file_name) %>%
  pivot_longer(cols = starts_with("MS2")) %>%
  separate(value, sep = ",", into = c("mz", "intensity")) %>%
  select(compound, mz, intensity, file_name)


# Switch Between Layouts --------------------------------------------------
# Scantable -> Concatenated
Scantable.to.Concatenated <- Scantable.Format %>%
  ungroup() %>%
  group_modify(~ScanToConcat(.x))

# Concatenated -> Scantable
Concatenated.to.Scantable <- Concatenated.Format %>%
  group_modify(~ConcatToScan(.x))

# Scale raw MS2 data --------------------------------------------------
ScaledforMS2 <- NonScaled.MS2 %>%
  rowwise() %>%
  mutate(MS2_2 = ScaleMS2(MS2))
