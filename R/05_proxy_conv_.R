# Author: Giacomo Bignardi
# Proxy test for convergence
# Use to calculate correlations between younger and older group of spouses of twins
library(tidyverse)
require(OpenMx)
set.seed(42)

# optimizer
mxOption(NULL, "Default optimizer", "NPSOL")

# clear wd
rm(list = ls())

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

# load extended twin fadmily data
dat_fam2g <- read.csv(sprintf("%s/%s/02_dat_twfam_estd.csv", wd_noa, wd_noa_output))
dat_fam2g <- dat_fam2g %>%
  rename(aes_tw1 = tw1, aes_tw2 = tw2, aes_s1 = s1, aes_stw1 = stw1, aes_stw2 = stw2) %>%
  select(-ends_with("s1"))

# pivot_longer to work at the partner level
dat_fam2g_long <- dat_fam2g %>%
  pivot_longer(
    cols = starts_with(c("age_", "sex_", "aes_")),
    names_to = c(".value", "type", "n"),
    names_pattern = "(age|sex|aes)_(tw|stw)([12])"
  ) %>%
  pivot_wider(
    names_from = "type",
    values_from = c("age", "sex", "aes")
  ) %>%
  # add median age of the couple
  rowwise() %>%
  mutate(age = median(age_tw, age_stw, na.rm = T)) %>%
  # split younger and older cohort
  ungroup() %>%
  mutate(cohort = ifelse(age < 50, 1, 2))

# number of dyads
table(dat_fam2g_long[!is.na(dat_fam2g_long$aes_tw) & !is.na(dat_fam2g_long$aes_stw), ]$cohort)

# select Variables for Analysis
vars <- "aes" # list of variables names
sel_vars <- paste0("aes", "_", c("tw", "stw"))
# select Data for Analysis
coh1 <- subset(dat_fam2g_long, cohort == 1, sel_vars)
coh2 <- subset(dat_fam2g_long, cohort == 2, sel_vars)

# number of complete observation
coh1_sample <- psych::pairwiseCount(coh1)
coh2_sample <- psych::pairwiseCount(coh2)
colnames(coh1_sample) <- sel_vars
rownames(coh1_sample) <- sel_vars
colnames(coh2_sample) <- sel_vars
rownames(coh2_sample) <- sel_vars

# MODEL--------------------------------------------------------------
# set starting values
m.va <- 2.5 # start value for means
v.va <- 1 # start value for variance
## model specification ####
# Create lgebra for expected mean matrices
meancoh1 <- mxMatrix(type = "Full", nrow = 1, ncol = 2, free = TRUE, values = m.va, labels = c("mtw1", "mstw1"), name = "meancoh1")
meancoh2 <- mxMatrix(type = "Full", nrow = 1, ncol = 2, free = TRUE, values = m.va, labels = c("mtw2", "mstw2"), name = "meancoh2")

# create algebra for expected variance/covariance matrices
# note that order on numbers index, number of pair (first and second) and number of cohort
covcoh1 <- mxMatrix(
  type = "Full", nrow = 2, ncol = 2, free = TRUE, values = diag(2),
  labels = c(
    "vtw1", "cstw1",
    "cstw1", "vstw1"
  ), name = "covcoh1"
)

# create data objects for multiple groups
covcoh2 <- mxMatrix(
  type = "Full", nrow = 2, ncol = 2, free = TRUE, values = diag(2),
  labels = c(
    "vtw2", "cstw2",
    "cstw2", "vstw2"
  ), name = "covcoh2"
)

# correlations
corcoh1 <- mxAlgebra(cov2cor(covcoh1), name = "corcoh1")
corcoh2 <- mxAlgebra(cov2cor(covcoh2), name = "corcoh2")

# create data objects for multiple groups
datacoh1 <- mxData(observed = coh1, type = "raw")
datacoh2 <- mxData(observed = coh2, type = "raw")

# create expectation objects for multiple groups
expcoh1 <- mxExpectationNormal(covariance = "covcoh1", means = "meancoh1", dimnames = sel_vars)
expcoh2 <- mxExpectationNormal(covariance = "covcoh2", means = "meancoh2", dimnames = sel_vars)
funML <- mxFitFunctionML()

# create model objects for multiple groups
modelcoh1 <- mxModel(meancoh1, covcoh1, corcoh1, datacoh1, expcoh1, funML, name = "coh1")
modelcoh2 <- mxModel(meancoh2, covcoh2, corcoh2, datacoh2, expcoh2, funML, name = "coh2")

multi <- mxFitFunctionMultigroup(c("coh1", "coh2"))

# create confidence interval objects
ciCov <- mxCI(c("coh1.covcoh1", "coh2.covcoh2"))
ciCor <- mxCI(c("coh1.corcoh1", "coh2.corcoh2"))
ciMean <- mxCI(c("coh1.meancoh1", "coh2.meancoh2"))

# build the extended saturated Model with confidence intervals
esat_mod <- mxModel("saturated", modelcoh1, modelcoh2, multi, ciCov, ciCor, ciMean)

# FIT--------------------------------------------------------------
# run saturated model
esat_fit <- mxRun(esat_mod, intervals = T)
esat_sumy <- summary(esat_fit)
esat_sumy$parameters
esat_sumy$CI
esat_sumy$parameters %>% arrange(name)

# CONSTRAINS--------------------------------------------------------------
esat_cage_mod <- mxModel(esat_fit, name = "cov age")
esat_cage_mod <- omxSetParameters(esat_cage_mod, label = c("cstw1", "cstw2"), free = TRUE, values = 1, newlabels = "cstw")

# fit and test model fit
esat_cage_fit <- mxRun(esat_cage_mod, intervals = T)
mxCompare(esat_fit, esat_cage_fit)
summary(esat_fit)$CI
