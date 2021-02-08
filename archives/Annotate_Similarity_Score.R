library(tidyverse)
source("KRH_Functions.R")
options(digits = 6)


# cosine1: dot product, compares, works well when you have fragments that match between two compounds. 
# sometimes when things fragments, you see difference between original/molecular mass, 
# and that difference is indicitive of class.

# cosine2: instead of comparing spectra on fragments, you compare on "shadow fragments", 
# which helps classify the compound. It's not useful when asking "are they the same compound", 
# but mostly when comparing two compounds in the same class; unknown spectra with similar neutral losses;
# trying to figure out what the loss is. We should keep it around, but can comment out for now.
# base mass - fragments for that. 

# Candidates$Cosine2 <- apply(Candidates, 1, FUN=function(x) MSMScosine2_df(x)). 
# Originally from the massbank mass positive match function in the MoNA_Matching script.

## Katherine csv: fully curated
MF_ClusterAssignments_Katherine <- read.csv("data_from_lab_members/MFCluster_Assignments_Katherine.csv") 
Cyano.MS2 <- read.csv("data_extra/Cyano_stds_withMS2s.csv") %>%
  mutate(Column = "Cyano")
HILICPos.MS2 <- read.csv("data_extra/HILICPos_stds_withMS2s.csv") %>%
  mutate(Column = "HILICPos")
HILICNeg.MS2 <- read.csv("data_extra/HILICNeg_stds_withMS2s.csv") %>%
  mutate(Column = "HILICNeg")

Standards.MS2 <- Cyano.MS2 %>% 
  rbind(HILICPos.MS2, HILICNeg.MS2) 

# Betaine Example ---------------------------------------------------
betaine_exp <- MF_ClusterAssignments_Katherine %>%
  filter(Identification == "Betaine") %>%
  select(-contains("cluster"))


# Uracil Example (from old MSP format) ---------------------------------------------------
Uracil.standard.msp <- read.delim("data_extra/Ingalls_HILICNeg_Standards(OLD_FORMAT).msp", header = FALSE, sep = "") %>%
  slice(56:66) %>%
  rename(mz = 1, intensity = 2) %>%
  mutate(mz = as.numeric(mz),
         intensity = as.numeric(intensity))

filtered.krh <- MF_ClusterAssignments_Katherine %>%
  filter(Identification == "Uracil") %>%
  select(MS2)
Uracil.KRH.msp <- scantable(filtered.krh)

test.uracil.msp.standards <- Uracil.standard.msp %>%
  mutate(intensity = round(intensity/max(intensity)*100, digits = 1)) %>%
  arrange(desc(intensity)) 

uracil.MS2cosine.sim <- MSMScosine_1(scan1 = Uracil.KRH.msp, scan2 = Uracil.standard.msp)

# MS1 Similarity
uracil.krh <- MF_ClusterAssignments_Katherine %>%
  filter(Identification == "Uracil") %>%
  select(Identification, mz, rt)

uracil.standard <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                            stringsAsFactors = FALSE, header = TRUE) %>%
  mutate(RT.seconds = RT..min. * 60) %>%
  select(-RT..min.) %>%
  filter(Compound.Name == "Uracil") %>%
  select(Compound.Name, m.z, RT.seconds)

exp.value.mz = uracil.krh$mz
theor.value.mz = uracil.standard$m.z

exp.value.rt = uracil.krh$rt
theor.value.rt = uracil.standard$RT.seconds

MS1.mz.similarity <- exp(-0.5 * (((exp.value.mz - theor.value.mz) / 0.02) ^ 2))
MS1.rt.similarity <- exp(-0.5 * (((exp.value.rt - theor.value.rt) / 0.4) ^ 2))

# Total Similarity Score
# TS = ((MS2 Similarity + MS1 Similarity) / 2) * 100
Uracil.Total.Similarity <- ((uracil.MS2cosine.sim + MS1.mz.similarity + MS1.rt.similarity) / 3) * 100

# Citrulline Example --------------------------------------------------------
Citrulline.standard.msp <- read.csv("data_extra/HILICPosMS2.csv") %>%
  filter(str_detect(Compound.Name, "Citrulline")) %>%
  select(MS2)
Citrulline.standard.msp <- scantable(Citrulline.standard.msp)

citrulline.filtered.krh <- MF_ClusterAssignments_Katherine %>%
  filter(Identification == "Citrulline") %>%
  select(MS2)
Citrulline.KRH.msp <- scantable(citrulline.filtered.krh)

citrulline.MS2cosine.sim <- MSMScosine_1(scan1 = Citrulline.KRH.msp, scan2 = Citrulline.standard.msp)

# MS1 Similarity
citrulline.krh <- MF_ClusterAssignments_Katherine %>%
  filter(Identification == "Citrulline") %>%
  select(Identification, mz, rt)

citrulline.standard <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                            stringsAsFactors = FALSE, header = TRUE) %>%
  mutate(RT.seconds = RT..min. * 60) %>%
  select(-RT..min.) %>%
  filter(Compound.Name == "Citrulline") %>%
  select(Compound.Name, m.z, RT.seconds)

exp.value.mz = citrulline.krh$mz
theor.value.mz = citrulline.standard$m.z

exp.value.rt = citrulline.krh$rt
theor.value.rt = citrulline.standard$RT.seconds

MS1.mz.similarity <- exp(-0.5 * (((exp.value.mz - theor.value.mz) / 0.02) ^ 2))
MS1.rt.similarity <- exp(-0.5 * (((exp.value.rt - theor.value.rt) / 0.4) ^ 2))

# Total Similarity Score
# TS = ((MS2 Similarity + MS1 Similarity) / 2) * 100
Citrulline.Total.Similarity <- ((uracil.MS2cosine.sim + MS1.mz.similarity + MS1.rt.similarity) / 3) * 100
