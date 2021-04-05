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

ScanToConcat <- function(scantable) {
  concatenated <- scantable %>%
    group_by(compound) %>%
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", "")),
           intensity = round(intensity / max(intensity) * 100, digits = 1)) %>%
    filter(intensity > 0.5)
    # select("mz", "intensity") %>%
    # unite("MS2", 1:2, sep = ", ") %>%
    # rownames_to_column() %>%
    # pivot_wider(names_from = 1, values_from = 2) %>%
    # unite("MS2", 1:ncol(.), sep = "; ") %>%
    # as.data.frame() #%>%

  return(concatenated)
}


# Both Layouts ------------------------------------------------------------
Concatenated.Format <- Cyano.MS2 %>%
  rename(MS2 = MS2s) %>%
  select(compound, MS2)

Scantable.Format <- cSplit(Concatenated.Format, "MS2", sep = ";") %>%
  pivot_longer(cols = starts_with("MS2")) %>%
  separate(value, sep = ",", into = c("mz", "intensity")) %>%
  select(compound, mz, intensity) %>%
  drop_na()


# Switch Between Layouts --------------------------------------------------
# Scantable -> Concatenated
Scantable.to.Concatenated <- Scantable.Format %>%
  group_modify(~ScanToConcat(.x))

# Concatenated -> Scantable
Concatenated.to.Scantable <- Concatenated.Format %>%
  group_modify(~ConcatToScan(.x))


