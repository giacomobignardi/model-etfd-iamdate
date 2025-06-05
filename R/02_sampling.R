# Author: Giacomo Bignardi
# Pseudo-randomized (optimized) sampling
# This scripts generate Supplementary Table 1
rm(list = ls())
library(tidyverse)

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

# load selected fam of twin data
dat <- read.csv(sprintf("%s/%s/01_aes_chills_estd.csv", wd_noa, wd_noa_output))
set.seed(42)
# whole sample descriptive
nrow(dat %>% distinct(FID))
# 10311 families
nrow(dat)
# 19404 total individuals
nrow(dat %>% filter(sex == 2))
# 12351 women
dat %>% reframe(range(age7, na.rm = T))
dat %>% reframe(range(age8, na.rm = T))
dat %>% reframe(range(age10, na.rm = T))
# > dat %>% reframe(range(age8, na.rm = T))
# range(age8, na.rm = T)
# 1                     11
# 2                     97

# sample descriptives for individuals with at least one measure across surveys
dat_complete <- dat %>% filter(!(is.na(neo43_7) & is.na(neo43_8) & is.na(neo43_10)))
nrow(dat_complete)
# 18474
nrow(dat_complete %>% filter(sex == 2))
# 11868

table(dat_complete$Extension)
# 1    2   10    61   62
# 7526 7377 2210 691  670

# 1 first multiple in family
# 2 2nd multiple in family
# 10 sibling
# 61 first registered spouse of person 1 in family
# 62 first registered spouse of person 2 in family

# SAMPLING --------------------------------------------------------------
# sampling procedure
# select twins and make a wide dataframe with one row per twin pair across waves of data collection
dat_tw <- dat %>%
  filter(!is.na(twzyg)) %>% # keep only twins with known zygosity
  rename(age_7 = age7, age_8 = age8, age_10 = age10) %>%
  dplyr::select(c(FID, Extension, twzyg, age_7, age_8, age_10, neo43_7, neo43_8, neo43_10)) %>%
  pivot_longer(age_7:neo43_10, names_to = c("var", "wave"), names_sep = "_", values_to = "score") %>%
  mutate(score = as.numeric(score)) %>%
  pivot_wider(names_from = c("var"), values_from = "score") %>%
  pivot_wider(names_from = "Extension", values_from = "age":"neo43") %>%
  rename(
    age_tw1 = age_1, age_tw2 = age_2,
    tw1 = neo43_1, tw2 = neo43_2
  ) %>%
  mutate(complete_data = is.na(tw1) + is.na(tw2)) # index for complete pair: 2 = both missing, 1 = one complete, 0 = both complete

# select pairs with complete data
comp_tw <- dat_tw %>%
  dplyr::filter(complete_data == 0) %>%
  group_by(FID) %>%
  slice_sample(n = 1) # sample only one per family across waves

# select pairs with data only for only one twin
sing_tw <- dat_tw %>%
  dplyr::filter(complete_data == 1) %>%
  group_by(FID) %>%
  slice_sample(n = 1) %>% # sample only one per family across waves
  filter(!(FID %in% comp_tw$FID))

# final data frame for twins
dat_twf <- rbind(comp_tw, sing_tw) %>%
  arrange(FID) %>%
  mutate( # add sex (to use as a covariate later MALE 1, FEMALE 2)
    sex_tw1 = ifelse(twzyg == 1 | twzyg == 2 | twzyg == 5, 1, 2),
    sex_tw2 = ifelse(twzyg == 1 | twzyg == 2 | twzyg == 6, 1, 2)
  ) %>%
  dplyr::select(-complete_data)

# final sample
table(dat_twf$twzyg)
# 1    2    3    4    5    6
# 1331  851 2883 1704 1270 1262

# impute missing age
# how many first-born twins have missing age? 18
sum(!is.na(dat_twf$tw1) & is.na(dat_twf$age_tw1))
# for these for how many I have age infomration on the other twin? 7
sum(!is.na(dat_twf$tw1) & is.na(dat_twf$age_tw1) & !is.na(dat_twf$age_tw2))

# how many second-born twins have missing age? 10
sum(!is.na(dat_twf$tw2) & is.na(dat_twf$age_tw2))
# for these for how many I have age infomration on the other twin? 7
sum(!is.na(dat_twf$tw2) & is.na(dat_twf$age_tw2) & !is.na(dat_twf$age_tw1))

# impute missing ages for the 8 first-born and 7 second born twins with missing ages
dat_twf[which(!is.na(dat_twf$tw1) & is.na(dat_twf$age_tw1) & !is.na(dat_twf$age_tw2)), ]$age_tw1 <- dat_twf[which(!is.na(dat_twf$tw1) & is.na(dat_twf$age_tw1) & !is.na(dat_twf$age_tw2)), ]$age_tw2
dat_twf[which(!is.na(dat_twf$tw2) & is.na(dat_twf$age_tw2) & !is.na(dat_twf$age_tw1)), ]$age_tw2 <- dat_twf[which(!is.na(dat_twf$tw2) & is.na(dat_twf$age_tw2) & !is.na(dat_twf$age_tw1)), ]$age_tw1

# remove individual with no age information
nrow(dat_twf)
dat_twf <- dat_twf %>%
  filter(!(!is.na(tw1) & is.na(age_tw1))) %>%
  filter(!(!is.na(tw2) & is.na(age_tw2)))

nrow(dat_twf)
table(dat_twf$twzyg)
# 1    2    3    4    5    6
# 1330  850 2882 1702 1267 1258

# make age equal across twins
dat_twf <- dat_twf %>%
  mutate(
    age_tw1 = ifelse(is.na(age_tw1), age_tw2, age_tw1),
    age_tw2 = ifelse(is.na(age_tw2), age_tw1, age_tw2)
  )

# remove individuals of age less then 18
nrow(dat_twf %>% filter(age_tw1 < 18))
table(dat_twf[dat_twf$age_tw1 < 18, ]$twzyg)
dat_twf <- dat_twf %>% filter(age_tw1 > 17)
table(dat_twf$twzyg)
# 1     2   3    4    5    6
# 1311  839 2857 1688 1264 1255
# simplify descriptives as we do not need information of opposite sex twins
dat_twf$twzyg <- ifelse(dat_twf$twzyg == 1 | dat_twf$twzyg == 3, "MZ", "DZ")

# descriptive twins
tw_des <- dat_twf %>%
  group_by(twzyg) %>%
  reframe(
    n = sum(!is.na(tw1)),
    mean = mean(age_tw1, na.rm = T),
    sd = sd(age_tw1, na.rm = T)
  )

# extend to family (sibs,spouses)
dat_fam <-
  dat %>%
  dplyr::filter(is.na(twzyg) & !(Extension == 1 | Extension == 2)) %>%
  rename(age_7 = age7, age_8 = age8, age_10 = age10) %>%
  dplyr::select(c(FID, sex, Extension, age_7, age_8, age_10, neo43_7, neo43_8, neo43_10)) %>%
  pivot_longer(age_7:neo43_10, names_to = c("var", "wave"), names_sep = "_", values_to = "score") %>%
  mutate(score = as.numeric(score)) %>%
  pivot_wider(names_from = c("var"), values_from = "score") %>%
  pivot_wider(names_from = "Extension", values_from = c("sex", "age":"neo43")) %>%
  rename(
    age_s1 = age_10, age_stw1 = age_61, age_stw2 = age_62,
    sex_s1 = sex_10, sex_stw1 = sex_61, sex_stw2 = sex_62,
    s1 = neo43_10, stw1 = neo43_61, stw2 = neo43_62
  ) %>%
  semi_join(dat_twf, by = c("FID", "wave")) # match to wave sampled for twins (twins are sampling units)
sum(!is.na(dat_fam$s1))
#[1] 1348
sum(!is.na(dat_fam$stw1))
#[1] 511
sum(!is.na(dat_fam$stw2))
#[1] 478

# create final extended twin family including partners
dat_twfam <- merge(dat_twf, dat_fam, c("FID", "wave"), all.x = T) %>%
  arrange(FID) %>%
  # årrange columns to be easier to read
  dplyr::select(
    FID, wave,
    twzyg,
    age_tw1, age_tw2, age_s1, age_stw1, age_stw2,
    sex_tw1, sex_tw2, sex_s1, sex_stw1, sex_stw2,
    tw1, tw2, s1, stw1, stw2
  )

# check for missing ages (here imputation is not possible)
sum(!is.na(dat_twfam$s1) & is.na(dat_twfam$age_s1))
# 1
dat_twfam[which(!is.na(dat_twfam$s1) & is.na(dat_twfam$age_s1)), ]$sex_s1
# male
# remove
dat_twfam[which(!is.na(dat_twfam$s1) & is.na(dat_twfam$age_s1)), ]$s1 <- NA
sum(!is.na(dat_twfam$stw1) & is.na(dat_twfam$age_stw1))
# 2
sum(!is.na(dat_twfam$stw2) & is.na(dat_twfam$age_stw2))
# 3
dat_twfam[which(!is.na(dat_twfam$stw1) & is.na(dat_twfam$age_stw1)), ]$sex_stw1
# 1 1
dat_twfam[which(!is.na(dat_twfam$stw2) & is.na(dat_twfam$age_stw2)), ]$sex_stw2
# 1 1 2

# remove
dat_twfam[which(!is.na(dat_twfam$stw1) & is.na(dat_twfam$age_stw1)), ]$stw1 <- NA
dat_twfam[which(!is.na(dat_twfam$stw2) & is.na(dat_twfam$age_stw2)), ]$stw2 <- NA

# remove individuals of age less then 18
nrow(dat_twfam %>% filter(age_s1 < 18))
table(dat_twfam[dat_twfam$age_s1 < 18, ]$sex_s1)
dat_twfam[which(!is.na(dat_twfam$s1) & dat_twfam$age_s1 < 18), ]$s1 <- NA
dat_twfam[which(dat_twfam$age_s1 < 18), ]$age_s1 <- NA

# SAMPLE --------------------------------------------------------------
# total sample
sum(!is.na(dat_twfam$tw1)) +
  sum(!is.na(dat_twfam$tw2)) +
  sum(!is.na(dat_twfam$s1)) +
  sum(!is.na(dat_twfam$stw1)) +
  sum(!is.na(dat_twfam$stw2))
#[1] 16583

# n of women
table(dat_twfam[which(!is.na(dat_twfam$tw1)), ]$sex_tw1) +
  table(dat_twfam[which(!is.na(dat_twfam$tw2)), ]$sex_tw2) +
  table(dat_twfam[which(!is.na(dat_twfam$s1)), ]$sex_s1) +
  table(dat_twfam[which(!is.na(dat_twfam$stw1)), ]$sex_stw1) +
  table(dat_twfam[which(!is.na(dat_twfam$stw2)), ]$sex_stw2)
# 1     2 
# 5819 10764 

# n of families
nrow(dat_twfam %>% distinct(FID))
#[1] 9214

# check that there is no individual younger then 18
range(dat_twfam$age_tw1, na.rm = T)
range(dat_twfam$age_tw2, na.rm = T)
range(dat_twfam$age_s1, na.rm = T)
range(dat_twfam$age_stw1, na.rm = T)
range(dat_twfam$age_stw2, na.rm = T)

# MZ
n_MZ_m <- table(dat_twfam[which(!is.na(dat_twfam$tw1)), ][dat_twfam$twzyg == "MZ", ]$sex_tw1)[1] + table(dat_twfam[which(!is.na(dat_twfam$tw2)), ][dat_twfam$twzyg == "MZ", ]$sex_tw2)[1]
n_MZ_f <- table(dat_twfam[which(!is.na(dat_twfam$tw1)), ][dat_twfam$twzyg == "MZ", ]$sex_tw1)[2] + table(dat_twfam[which(!is.na(dat_twfam$tw2)), ][dat_twfam$twzyg == "MZ", ]$sex_tw2)[2]
# DZ
n_DZ_m <- table(dat_twfam[which(!is.na(dat_twfam$tw1)), ][dat_twfam$twzyg == "DZ", ]$sex_tw1)[1] + table(dat_twfam[which(!is.na(dat_twfam$tw2)), ][dat_twfam$twzyg == "DZ", ]$sex_tw2)[1]
n_DZ_f <- table(dat_twfam[which(!is.na(dat_twfam$tw1)), ][dat_twfam$twzyg == "DZ", ]$sex_tw1)[2] + table(dat_twfam[which(!is.na(dat_twfam$tw2)), ][dat_twfam$twzyg == "DZ", ]$sex_tw2)[2]
# Full-sibs
n_s_m <- table(dat_twfam[which(!is.na(dat_twfam$s1)), ]$sex_s1)[1]
n_s_f <- table(dat_twfam[which(!is.na(dat_twfam$s1)), ]$sex_s1)[2]
# husbands
n_h <- table(dat_twfam[which(!is.na(dat_twfam$stw1)), ]$sex_stw1)[1] + table(dat_twfam[which(!is.na(dat_twfam$stw2)), ]$sex_stw2)[1]
# wifes
n_w <- table(dat_twfam[which(!is.na(dat_twfam$stw1)), ]$sex_stw1)[2] + table(dat_twfam[which(!is.na(dat_twfam$stw2)), ]$sex_stw2)[2]

# final sample size
dem <- data.frame(
  N = c(n_MZ_m, n_MZ_f, n_DZ_m, n_DZ_f, n_s_m, n_s_f, n_h, n_w),
  Type = c("MZ", "MZ", "DZ", "DZ", "Full sib", "Full sib", "partner", "partner"),
  Sex = c("men", "women", "men", "women", "men", "women", "man", "women")
)

# SAVE --------------------------------------------------------------
# demographic
write_csv(dem, sprintf("%s/SI/02_table1_dem.csv", wd_oa))

## sampled family twin data ####
write_csv(dat_twfam, sprintf("%s/%s/02_dat_twfam_estd.csv", wd_noa, wd_noa_output))