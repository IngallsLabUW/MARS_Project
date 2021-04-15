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
    separate(value, sep = ",", into = c("mz", "intensity")) #%>%
    #select(compound, mz, intensity, file_name)
  
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
  mutate(MS2_scaled = ScaleMS2(MS2))

################################################################################################################

PlotScaled <- function(compound.name) {
  scaled <- ScaledforMS2 %>%
    filter(grepl("Mix1", filename)) %>%
    filter(compound_name == compound.name & voltage == 20) %>%
    select("compound" = compound_name, "file_name" = filename, voltage, "MS2" = MS2_scaled) %>%
    group_modify(~ConcatToScan(.x)) %>%
    mutate(B = substr(file_name, nchar(file_name)-9+1, nchar(file_name))) %>%
    mutate(run = substr(B,1,nchar(B)-5)) %>%
    select(-name, -B, -file_name) %>%
    unique() %>%
    drop_na() %>%
    mutate(intensity = as.numeric(intensity)) %>%
    mutate(mz = as.numeric(mz))
  
  plot <- ggplot(scaled, aes(x = mz, y = intensity, group = run)) +
    #geom_line(aes(color = run), size=1.2) +
    geom_point(aes(color = run)) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(paste(compound.name, "20 volts, PosMix1, Scaled"))
  
  print(plot)
}

plot1 <- PlotScaled("4-Aminobutyric acid")
plot2 <- PlotScaled("5-Hydroxyectoine")
plot3 <- PlotScaled("Adenine")
plot4 <- PlotScaled("Isocitric acid")

require(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4, ncol=2)





## garbage

adenine.raw <- NonScaled.MS2 %>%
  filter(compound_name == "Adenine" & voltage == 20) %>%
  select("compound" = compound_name, "file_name" = filename, voltage, MS2) %>%
  group_modify(~ConcatToScan(.x)) %>%
  filter(grepl("Mix1", file_name)) %>%
  mutate(B = substr(file_name, nchar(file_name)-9+1, nchar(file_name))) %>%
  mutate(run = substr(B,1,nchar(B)-5)) %>%
  select(-name, -B, -file_name) %>%
  unique() %>%
  drop_na()

t.raw <- ggplot(adenine.raw, aes(x = mz, y = intensity, group = run)) +
  geom_line(aes(color = run), size=1.2) +
  geom_point() +
  ggtitle("Adenine, 20 volts, Mix1, Raw Data")
t.raw

