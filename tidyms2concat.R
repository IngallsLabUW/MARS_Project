library(splitstackshape)
library(tidyverse)

## Tidy MS2 consensus
ConcatToScan <- function(concat.format) {
  scantable <- cSplit(concat.format, "MS2", sep = ";") %>%
    pivot_longer(cols = starts_with("MS2")) %>%
    separate(value, sep = ",", into = c("mz", "intensity")) %>%
    rename(fragment = name)
  
  return(scantable)
}

Filterppm <- function(mass, ppm) {
  
  filtered.mass <- c(mass * (1 - ppm/1e+06), mass * (1 + ppm/1e+06))
  
  return(filtered.mass)
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

ScanToConcat <- function(scantable.format) {
  concatenated <- scantable.format %>%
    unique() %>%
    group_by(compound_name, filename) %>%
    unite("MS2", mz:intensity, sep = ",") %>%
    ungroup() %>%
    group_by(compound_name, voltage) %>%
    mutate(MS2 = paste(MS2, collapse = "; ")) %>%
    select(-fragment) %>%
    unique() %>%
    as.data.frame()
  
  return(concatenated)
}

TakeConsensus <- function(msms_df) { ## Different length issues- check
  
  msms_df <- msms_df %>%
    mutate(mz = as.numeric(mz),
           intensity = as.numeric(intensity))
  
  consensus <- list()
  while(nrow(msms_df) > 0) {
    first_frag_mz <- msms_df$mz[1]
    frag_data <- msms_df %>%
      filter(mz%between%Filterppm(first_frag_mz, 5)) ## Two fragments within the same window
    if(nrow(frag_data) >= 3){
      consensus[[length(consensus)+1]] <- c(
        unique(frag_data$compound_name),
        unique(frag_data$voltage),
        unique(frag_data$fragment),
        mean(frag_data$mz), mean(frag_data$intensity))
    }
    msms_df <- anti_join(msms_df, frag_data)
  }
  consensus.df <- data.frame(matrix(unlist(consensus), nrow=length(consensus), byrow=TRUE)) #%>%
  # rename(compound_name = 1,
  #        voltage = 2,
  #        mz = 3,
  #        intensity = 4)
  
  return(consensus.df)
}


# Import MSMS data that has been run 5 times
# Scale and adjust format
Fivetimes <- read.csv("~/Downloads/Ingalls_Lab_Standards_MSMS.csv") %>%
  rowwise() %>%
  mutate(MS2 = ScaleMS2(MS2))

Fivetimes.Scantable <- Fivetimes %>%
  ConcatToScan() %>%
  drop_na()

Fivetimes.Concat <- Fivetimes.Scantable %>%
  ScanToConcat()
  

Fivetimes.Wide <- Fivetimes %>%
  mutate(B = substr(filename, nchar(filename)-9+1, nchar(filename))) %>%
  mutate(run = substr(B,1,nchar(B)-5)) %>%
  select(-B, -filename) %>%
  cSplit(2, sep = ";") %>%
  cSplit(4:ncol(.), sep = ", ")
names(Fivetimes.Wide) <- gsub(x = names(Fivetimes.Wide), pattern = "_1", replacement = "_mz")  
names(Fivetimes.Wide) <- gsub(x = names(Fivetimes.Wide), pattern = "_2", replacement = "_int")  

Fivetimes.Var <- Fivetimes.Wide %>%
  group_by(voltage, compound_name) %>%
  mutate_at(vars(contains("MS2")), funs(var(.))) %>%
  select(-run) %>%
  unique() %>%
  select(voltage, compound_name, contains("mz")) %>%
  pivot_longer(
    cols = starts_with("MS2"), names_to = "mz", names_prefix = "wk",
    values_to = "mz_variation", values_drop_na = TRUE
  ) %>%
  mutate(voltage = factor(voltage)) %>%
  mutate(sqrt_mz_variation = sqrt(mz_variation)) %>%
  filter(!grepl("int", mz))

## Visualize variation
var.heatplot <- ggplot(Fivetimes.Var, aes(x = compound_name, y = voltage, fill = mz_variation)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))
var.heatplot

var.heatplot.sqrt <- ggplot(Fivetimes.Var, aes(x = compound_name, y = voltage, fill = sqrt_mz_variation)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))
var.heatplot.sqrt

var.line <- ggplot(data = Fivetimes.Var,
                   mapping = aes(x = reorder(compound_name, mz_variation), mz_variation)) +
  geom_line() +
  theme(axis.text.x = element_blank()) +
  xlab("compound")
var.line

## Specific troublesome compounds
# All compounds with a variation over 10,000
issues <- Fivetimes.Var %>%
  arrange(desc(mz_variation)) %>%
  filter(mz_variation < 1)
print(unique(issues$compound_name))
## Check this for non-manual adjustment of incorrect fragments

# Taking consensus --------------------------------------------------------

# Success on mini df: But for some reason 35 creates unusual intensity average??
my.compound <- "(3-Carboxypropyl)trimethylammonium"
my.voltage <- 20
my.compound <- "Creatine"
my.voltage <- 100
compound.subset <- Fivetimes.Scantable %>%
  filter(voltage == my.voltage,
         compound_name == my.compound) %>%
  arrange(fragment)

compound.consensus <- TakeConsensus(compound.subset)


# Full attempt
# full <- Fivetimes.Scantable %>%
#   filter(compound_name %in% c("(3-Carboxypropyl)trimethylammonium", "Isocitric acid")) %>%
#   group_by(voltage, compound_name) %>%
#   arrange(voltage, fragment) %>%
#   group_split()
#
# full.split <- lapply(full, TakeConsensus)
#
# df <- bind_rows(t)

