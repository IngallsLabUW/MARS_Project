library(tidyverse)

## Issues with standards sheet: isolating rows 


Standards <- read.csv("data_extra/Ingalls_Lab_Standards_Jan11.csv", stringsAsFactors = FALSE) %>%
  select(Compound.Name, Compound.Name_old, Column, Fraction1, Fraction2, Liquid.Fraction, ionization_form, z)

# Issues ------------------------------------------------------------------

# Glutamylphenylalanine: HILICPos, negative z, HILICPos, postive z. Not a real issue due to ionizatin forms, but it is the only compound that does that.
Fraction1Z <- Standards %>%
  mutate(test = ifelse(Fraction1 == "HILICNeg" & z == 1, TRUE,
                       ifelse(Fraction1 == "HILICPos" & z == -1, TRUE, FALSE))) %>%
  filter(test == TRUE)

# Any way to check that the [z] is correct on these?
OnlyM <- Standards %>%
  filter(ionization_form == "[M]")

# What about when ionization_form is NA?
NoM <- Standards %>%
  filter(is.na(ionization_form))

# Mismatching ionization_forms and z
IonForms  <- Standards %>%
  filter(!Compound.Name %in% OnlyM$Compound.Name) %>%
  mutate(PosNeg = ifelse(str_detect(ionization_form, "-"), TRUE, FALSE)) %>%
  filter(PosNeg == TRUE & z == 1 | PosNeg == FALSE & z == -1) %>%
  select(-PosNeg)

# No problems in the Column, Fraction 1/2 setup in HILIC or Cyano, fraction 1 and 2 do not overlap at any point.
# In the RP Column:
# CyanoDCM --> CyanoAq (fraction 1 --> fraction 2) has only Org liquid fractions, but the other way around,
# CyanoAq --> CyanoDCM has a fair amount of both Aq and Org liquid fractions
# Is this correct? If we ditch fracitons1/2 will that be a problem?

LiquidFractions <- Standards %>%
  filter(Column == "RP")
