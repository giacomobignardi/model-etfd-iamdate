# Author: Giacomo Bignardi
# Adapted from Bignardi et al., https://doi.org/10.1038/s41598-022-07161-z
# prepare NTR data for etfd-based modeling
rm(list = ls())
library(tidyverse)
library(haven) # required to read .sav

# set open access working directories
wd_oa <- getwd()
wd_oa_scripts <- "R"
wd_oa_figures <- "Figures"

# set other working directories
wd_noa <- substr(
  getwd(),
  0,
  nchar(getwd()) - nchar("11_github") - 1
)
wd_noa_data <- "02_data"
wd_noa_output <- "03_output"

# load twin data obtained from the NTR
data_aes_chills <- read_sav(sprintf("%s/%s/Giacomo\ Bignardi_Heritability\ of\ Aesthetic\ chills.sav", wd_noa, wd_noa_data))

# n indiviudals (per sex)
table(data_aes_chills$sex)

# range age
range(data_aes_chills$age7, na.rm = T)
range(data_aes_chills$age8, na.rm = T)
range(data_aes_chills$age10, na.rm = T)

# n families
data_aes_chills %>%
  distinct(FID) %>%
  nrow()

# select family members
aes_chills <- data_aes_chills %>%
  dplyr::filter(Extension == 1 | Extension == 2 | Extension == 10 | Extension == 61 | Extension == 62) %>% # select twins
  # 1 first multiple in family
  # 2 2nd multiple in family
  # 10 sibling
  # 61 first registered spouse of person 1 in family
  # 62 first registered spouse of person 2 in family
  dplyr::filter(multiple_type <= 2) %>% # filter for < triplet
  dplyr::select(c(FID, Extension, twzyg, sex, age7, age8, age10, neo43_7, neo43_8, neo43_10)) %>%
  dplyr::mutate(neo43_7 = as.numeric(neo43_7), neo43_8 = as.numeric(neo43_8), neo43_10 = as.numeric(neo43_10)) %>%
  #   value label
  # 1   MZM
  # 2   DZM
  # 3   MZF
  # 4   DZF
  # 5 DOSmf
  # 6 DOSfm
  dplyr::filter(twzyg == 1 | twzyg == 2 | twzyg == 3 | twzyg == 4 | twzyg == 5 | twzyg == 6 | is.na(twzyg)) %>%
  # 1 male 2 female
  dplyr::filter(sex == 1 | sex == 2)

# SAVE --------------------------------------------------------------
## family twin data ####
write.csv(as.data.frame(aes_chills),sprintf("%s/%s/01_aes_chills_estd.csv",wd_noa,wd_noa_output))