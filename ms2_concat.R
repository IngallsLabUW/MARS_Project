library(splitstackshape)
library(tidyverse)

Cyano.MS2 <- read.csv("data_extra/Standards_MS2/Cyano_stds_withMS2s.csv") %>%
  mutate(column = "RP",
         z = NA)


# Functions ---------------------------------------------------------------
ConcatToScan <- function(concat.format) {
  final <- cSplit(concat.format, "MS2", sep = ";") %>%
    pivot_longer(cols = starts_with("MS2")) %>%
    separate(value, sep = ",", into = c("mz", "intensity")) %>%
    select(compound, mz, intensity) %>%
    drop_na()
  
  return(final)
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
    as.data.frame()
  
  return(concatenated)
}

# mutate(mz = as.numeric(mz %>% str_replace(",", "")),
#        intensity = as.numeric(intensity %>% str_replace(";", "")),
#        intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
# filter(intensity > 0.5) %>%


# Both Layouts ------------------------------------------------------------
Concatenated.Format <- Cyano.MS2 %>%
  rename(MS2 = MS2s) %>%
  select(compound, MS2, file_name)

Scantable.Format <- cSplit(Concatenated.Format, "MS2", sep = ";") %>%
  group_by(compound, file_name) %>%
  pivot_longer(cols = starts_with("MS2")) %>%
  separate(value, sep = ",", into = c("mz", "intensity")) %>%
  select(compound, mz, intensity, file_name) #%>%
#drop_na()

# Switch Between Layouts --------------------------------------------------
# Scantable -> Concatenated
Scantable.to.Concatenated <- Scantable.Format %>%
  ungroup() %>%
  group_modify(~ScanToConcat(.x))

# Concatenated -> Scantable
Concatenated.to.Scantable <- Concatenated.Format %>%
  group_modify(~ConcatToScan(.x))
