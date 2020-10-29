library(tidyverse)


# Various untargeted outputs from lab members -----------------------------

# Will csv: output from xcms
Addiso_Will <- read.csv("data_from_lab_members/addiso_features_Will.csv")
Complete_Features_Will <- read.csv("data_from_lab_members/complete_features_Will.csv")
Final_Peaks_Will <- read.csv("data_from_lab_members/final_peaks_Will.csv")

# Josh text: MSDial
ET_HILICPos_Josh <- read.delim("data_from_lab_members/ET_HILIC_Pos_Area_Josh.txt")
ET_HILICNeg_Josh <- read.delim("data_from_lab_members/ET_HILIC_Neg_Area_Josh.txt")

# Angie csv + text: MSDial and xcms
XCMS_features_info_Angie <- read.csv("data_from_lab_members/XCMS-filled-features-info-RTSubset_.csv")
XCMS_features_values_Angie <- read.csv("data_from_lab_members/XCMS-filled-features-values-RTSubset_.csv")
Area_Angie <- read.delim("data_from_lab_members/Area_0_202010201051.txt")
