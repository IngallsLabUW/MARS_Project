library(splitstackshape)
library(tidyverse)
options(scipen = 999)

Cyano.MS2 <- read.csv("data_extra/Standards_MS2/Cyano_stds_withMS2s.csv") %>%
  mutate(column = "RP",
         z = NA)

#NonScaled.MS2 <- read.csv("~/work/Ingalls_Standards/MSMS/data_processed/Ingalls_Lab_Standards_MSMS.csv")
NonScaled.MS2 <- read.csv("~/Downloads/Ingalls_Lab_Standards_MSMS.csv") #5x

#Fivetimes <- read.csv("data_extra/Ingalls_Lab_Standards_MSMS5x.csv")
Fivetimes <- read.csv("~/Downloads/Ingalls_Lab_Standards_MSMS.csv")

# Functions ---------------------------------------------------------------
ConcatToScan <- function(concat.format) {
  scantable <- cSplit(concat.format, "MS2", sep = ";") %>%
    pivot_longer(cols = starts_with("MS2")) %>%
    separate(value, sep = ",", into = c("mz", "intensity")) #%>%
  #select(compound, mz, intensity, file_name)
  
  return(scantable)
}

PlotScaled <- function(compound.name, df) {
  # df must be in concatenated, scaled form
  scaled <- df %>%
    filter(grepl("Mix1", filename)) %>%
    filter(compound_name == compound.name & voltage == 20) %>%
    select("compound" = compound_name, "file_name" = filename, voltage, "MS2" = contains("MS2")) %>%
    group_modify(~ConcatToScan(.x)) %>%
    mutate(B = substr(file_name, nchar(file_name)-9+1, nchar(file_name))) %>%
    mutate(run = substr(B,1,nchar(B)-5)) %>%
    select(-name, -B, -file_name) %>%
    unique() %>%
    drop_na() %>%
    mutate(intensity = as.numeric(intensity)) %>%
    mutate(mz = as.numeric(mz)) %>%
    filter(mz > 84 & mz < 86)
  
  plot <- ggplot(scaled) +
    geom_segment(aes(x = mz, y = intensity, group = run, xend = mz, yend = 0)) +
    facet_wrap(~run) +
    #geom_line(aes(color = run), size=1.2) +
    #geom_point(aes(color = run)) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(paste(compound.name, "20 volts, PosMix1, Scaled"))
  
  print(plot)
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
  drop_na() %>%
  rowwise() %>%
  mutate(MS2_scaled = ScaleMS2(MS2))

################################################################################

# MS2 five time run variation --------------------------------------------------
Fivetimes.Scantable <- Fivetimes %>%
  rowwise() %>%
  mutate(MS2 = ScaleMS2(MS2))

t <- Fivetimes.Scantable %>%
  mutate(B = substr(filename, nchar(filename)-9+1, nchar(filename))) %>%
  mutate(run = substr(B,1,nchar(B)-5)) %>%
  select(-B, -filename) %>%
  cSplit(2, sep = ";") %>%
  cSplit(4:ncol(.), sep = ", ")
names(t) <- gsub(x = names(t), pattern = "_1", replacement = "_mz")  
names(t) <- gsub(x = names(t), pattern = "_2", replacement = "_int")  

t2 <- t %>%
  arrange(voltage) %>%
  group_by(voltage, compound_name) %>%
  mutate_at(vars(contains("MS2")), funs(var(.))) %>%
  select(-run) %>%
  unique()

t3 <- t2 %>%
  ungroup() %>%
  select(voltage, compound_name, contains("mz")) %>%
  pivot_longer(
    cols = starts_with("MS2"),
    names_to = "mz",
    names_prefix = "wk",
    values_to = "mz_variation",
    values_drop_na = TRUE
  ) %>%
  mutate(voltage = factor(voltage)) %>%
  mutate(sqrt_mz_variation = sqrt(mz_variation)) %>%
  filter(!grepl("int", mz))

p <- ggplot(t3, aes(x = compound_name, y = voltage, fill = mz_variation)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))
p

p.sqrt <- ggplot(t3, aes(x = compound_name, y = voltage, fill = sqrt_mz_variation)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))
p.sqrt

ggplot(data = t3,
       mapping = aes(x = reorder(compound_name, mz_variation), mz_variation)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))

# Plot attempts -----------------------------------------------------------
# Base R
compound.name <- "Adenine"
voltage <- 50
My.Compound <- Fivetimes.Scaled %>%
  filter(grepl("Mix1", filename)) %>%
  filter(compound_name == compound.name) %>%
  select("compound" = compound_name, "file_name" = filename, voltage, "MS2" = contains("MS2")) %>%
  group_modify(~ConcatToScan(.x)) %>%
  mutate(B = substr(file_name, nchar(file_name)-9+1, nchar(file_name))) %>%
  mutate(run = substr(B,1,nchar(B)-5)) %>%
  select(-name, -B, -file_name) %>%
  unique() %>%
  drop_na() %>%
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(mz = as.numeric(mz)) %>%
  filter(voltage == 50)

# Single MS2 plot
plot(intensity~mz, type="h", data=My.Compound, ylab="Intensity", xlab="Fragment m/z", col = "blue", lwd = 4)
title(paste("Scaled,", compound.name, voltage, "volts, Mix 1"))
legend("topleft", legend = unique(My.Compound$compound))

plot1 <- PlotScaled("Isocitric acid", Fivetimes.Scaled)
plot2 <- PlotScaled("Succinylglycine", Fivetimes.Scaled)
plot3 <- PlotScaled("Adenine", Fivetimes.Scaled)
plot4 <- PlotScaled("Isocitric acid", Fivetimes.Scaled)

require(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4, ncol=2)

################################################################################################################

######## Example code
library(tidyverse)
library(data.table)
pmppm <- function (mass, ppm = 4) {
  return(c(mass * (1 - ppm/1e+06), mass * (1 + ppm/1e+06)))
}
frag1_msms_data <- data.frame(
  mz=100+runif(5, -0.00001, 0.00001),
  int=100+runif(5, -2, 0)
)
frag2_msms_data <- data.frame(
  mz=50+runif(5, -0.00001, 0.00001),
  int=50+runif(5, -2, 0)
)
frag3_msms_data <- data.frame(
  mz=75,
  int=1
)

msms_data <- rbind(frag1_msms_data, frag2_msms_data, frag3_msms_data) %>%
  arrange(desc(int))

TakeConsensus <- function(msms_df) {
  
  consensus <- list()
  while(nrow(msms_df)>0){
    first_frag_mz <- msms_df$mz[1]
    frag_data <- msms_df %>%
      filter(mz%between%pmppm(first_frag_mz, 5))
    if(nrow(frag_data)>=3){
      consensus[[length(consensus)+1]] <- c(mean(frag_data$mz), mean(frag_data$int))
    }
    msms_df <- anti_join(msms_df, frag_data)
  }
  return(consensus)
}

test_msms <- t %>%
  filter(compound_name == "(3-Carboxypropyl)trimethylammonium") %>%
  group_by(voltage)

