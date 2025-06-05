# Author: Giacomo Bignardi
# Saturated model
# Use to calculate correlations between family members
# Adapted from https://hermine-maes.squarespace.com/#/one/ saturated models
rm(list = ls())
library(tidyverse)
require(OpenMx)

# Optimizer
mxOption(NULL, "Default optimizer", "NPSOL")
set.seed(42)

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
dat_fam2g <- dat_fam2g %>% rename(aes_tw1 = tw1, aes_tw2 = tw2, aes_s1 = s1, aes_stw1 = stw1, aes_stw2 = stw2)

# make data long
dat_fam2g_long <- dat_fam2g %>% pivot_longer(cols = age_tw1:aes_stw2, names_to = c(".value", "rel"), names_sep = "_")

# residualise for sex and age
dat_fam2g_long$aes_res <- residuals(lm(aes ~ sex + age, data = dat_fam2g_long, na.action = na.exclude))

# pivot wider
dat_fam2g <- dat_fam2g_long %>%
  dplyr::select(FID, twzyg, rel, aes_res) %>%
  pivot_wider(names_from = rel, values_from = aes_res)
dat_fam2g <- dat_fam2g %>% rename(aes_tw1 = tw1, aes_tw2 = tw2, aes_s1 = s1, aes_stw1 = stw1, aes_stw2 = stw2)

# select Variables for Analysis
vars <- "aes" # list of variables names
sel_vars <- paste0("aes", "_", c("stw1", "tw1", "s1", "tw2", "stw2"))
# select Data for Analysis
mz <- subset(dat_fam2g, twzyg == "MZ", sel_vars)
dz <- subset(dat_fam2g, twzyg == "DZ", sel_vars)

# DYADS--------------------------------------------------------------
# number of complete observation
mz_sample <- psych::pairwiseCount(mz)
dz_sample <- psych::pairwiseCount(dz)

colnames(mz_sample) <- sel_vars
rownames(mz_sample) <- sel_vars
colnames(dz_sample) <- sel_vars
rownames(dz_sample) <- sel_vars

# cor MZ twin
nmz_comp_pair <- mz_sample["aes_tw1", "aes_tw2"]

# cor DZ twin
ndz_comp_pair <- dz_sample["aes_tw1", "aes_tw2"]

# cor FS
nfs_comp_pair <- mz_sample["aes_s1", "aes_tw1"] + mz_sample["aes_s1", "aes_tw2"] + dz_sample["aes_s1", "aes_tw1"] + dz_sample["aes_s1", "aes_tw2"]

# cor spouses (am)
nstw_comp_pair <- mz_sample["aes_tw1", "aes_stw1"] + dz_sample["aes_tw2", "aes_stw2"] + mz_sample["aes_tw2", "aes_stw2"] + dz_sample["aes_tw1", "aes_stw1"]

# cor MZ twin in law
nmzil_comp_pair <- mz_sample["aes_tw1", "aes_stw2"] + mz_sample["aes_tw2", "aes_stw1"]

# cor DZ twin in law
ndzil_comp_pair <- dz_sample["aes_tw1", "aes_stw2"] + dz_sample["aes_tw2", "aes_stw1"]
nfsil_comp_pair <- mz_sample["aes_s1", "aes_stw1"] + mz_sample["aes_s1", "aes_stw2"] + dz_sample["aes_s1", "aes_stw1"] + dz_sample["aes_s1", "aes_stw2"]

# cor twin in co-law
nmzcil_comp_pair <- mz_sample["aes_stw1", "aes_stw2"]

# cor twin in co-law
ndzcil_comp_pair <- dz_sample["aes_stw1", "aes_stw2"]

# save for later
n_compair <- data.frame(
  # cor MZ 
  nmz_comp_pair,
  # cor DZ 
  ndz_comp_pair,
  # cor FS 
  nfs_comp_pair,
  # cor spouses (am)
  nstw_comp_pair,
  # cor MZ twin in law
  nmzil_comp_pair,
  # cor DZ twin in law
  ndzil_comp_pair,
  # cor FS twin in law
  nfsil_comp_pair,
  # cor mz in co-law
  nmzcil_comp_pair,
  # cor dz in co-law
  ndzcil_comp_pair
)

# save complete pair to add to later table
write_csv(n_compair, sprintf("%s/%s/03_pair_number.csv", wd_noa, wd_noa_output))

# check that we are capturing all relationships (total number of pairwise complete observation)
sum(mz_sample[lower.tri(mz_sample)]) + sum(dz_sample[lower.tri(dz_sample)]) ==

  # cor MZ 
  nmz_comp_pair +
    # cor DZ 
    ndz_comp_pair +
    # cor FS 
    nfs_comp_pair +
    # cor spouses (am)
    nstw_comp_pair +
    # cor MZ twin in law
    nmzil_comp_pair +
    # cor DZ twin in law
    ndzil_comp_pair +
    # cor FS twin in law
    nfsil_comp_pair +
    # cor twin in co-law
    nmzcil_comp_pair +
    ndzcil_comp_pair

# total number of individuals
sum(diag(mz_sample)) + sum(diag(dz_sample))
# 16583

# MODEL--------------------------------------------------------------
# set starting Values
m.va <- 0 # start value for means
v.va <- 1 # start value for variance
lb.Va <- .0001 # lower bound for variance
## Model specification####
# create algebra for expected mean matrices
meanmz <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = TRUE, values = m.va, labels = c("mmzstw1", "mmz1", "mmzfs", "mmz2", "mmzstw2"), name = "meanmz")
meandz <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = TRUE, values = m.va, labels = c("mdzstw1", "mdz1", "mdzfs", "mdz2", "mdzstw2"), name = "meandz")

# create algebra for expected variance/covariance matrices
covmz <- mxMatrix( type="Full", nrow=5, ncol=5, free=TRUE, values=diag(5),
                   labels=c("vmzstw1",   "cmzstw1", "cmzfsil1",  "cmzil2",   "cmzstwcil",
                            "cmzstw1",   "vmz1",    "cmzfs1",    "cmz",      "cmzil1",
                            "cmzfsil1",  "cmzfs1",  "vmzfs",     "cmzfs2",  "cmzfsil2",
                            "cmzil2",    "cmz",     "cmzfs2",   "vmz2",     "cmzstw2",
                            "cmzstwcil", "cmzil1",  "cmzfsil2",  "cmzstw2",  "vmzstw2"), name="covmz" )

# create data objects for multiple groups
covdz <- mxMatrix( type="Full", nrow=5, ncol=5, free=TRUE, values=diag(5),
                   labels=c("vdzstw1",   "cdzstw1", "cdzfsil1",   "cdzil2",  "cdzstwcil",
                            "cdzstw1",   "vdz1",    "cdzfs1",     "cdz",     "cdzil1",
                            "cdzfsil1",  "cdzfs1",  "vdzfs",      "cdzfs2", "cdzfsil2",
                            "cdzil2",    "cdz",     "cdzfs2",     "vdz2",    "cdzstw2",
                            "cdzstwcil", "cdzil1",  "cdzfsil2",   "cdzstw2", "vdzstw2"), name="covdz" )

# correlations
cormz <- mxAlgebra(cov2cor(covmz), name = "cormz")
cordz <- mxAlgebra(cov2cor(covdz), name = "cordz")

# create data objects for multiple groups
datamz <- mxData(observed = mz, type = "raw")
datadz <- mxData(observed = dz, type = "raw")

# create expectation objects for multiple groups
expmz <- mxExpectationNormal(covariance = "covmz", means = "meanmz", dimnames = sel_vars)
expdz <- mxExpectationNormal(covariance = "covdz", means = "meandz", dimnames = sel_vars)
funML <- mxFitFunctionML()

# create model objects for multiple groups
modelmz <- mxModel(meanmz, covmz, cormz, datamz, expmz, funML, name = "mz")
modeldz <- mxModel(meandz, covdz, cordz, datadz, expdz, funML, name = "dz")

multi <- mxFitFunctionMultigroup(c("mz", "dz"))

# create confidence interval objects
ciCov <- mxCI(c("mz.covmz", "dz.covdz"))
ciCor <- mxCI(c("mz.cormz", "dz.cordz"))
ciMean <- mxCI(c("mz.meanmz", "dz.meandz"))

# build the extended saturated model with confidence intervals
esat_mod <- mxModel("saturated", modelmz, modeldz, multi, ciCov, ciCor, ciMean)

# FIT--------------------------------------------------------------
# run saturated model
esat_fit <- mxRun(esat_mod, intervals = F)
esat_sumy <- summary(esat_fit)
esat_sumy$parameters
esat_sumy$CI
esat_sumy$parameters %>% arrange(name)

save(esat_fit, file = sprintf("%s/%s/03_saturated.Rdata", wd_noa, wd_noa_output))

# CONSTRAINS--------------------------------------------------------------
# constrain mean and variances and covariances to have a reduced number of correlations
# constrain means across twin order
esat_mvar_mod <- mxModel(esat_fit, name = "means twin order")

# equate means and variances to obtain summary statistics on the following:
# twin order means
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("mmz1", "mmz2"), free = TRUE, values = 1, newlabels = "mmz1")
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("mdz1", "mdz2"), free = TRUE, values = 1, newlabels = "mdz1")
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("mdzstw1", "mdzstw2"), free = TRUE, values = 1, newlabels = "mdzstw1")
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("mmzstw1", "mmzstw2"), free = TRUE, values = 1, newlabels = "mmzstw1")

# twin order var
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("vmz1", "vmz2"), free = TRUE, values = 1, newlabels = "vmz1")
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("vdz1", "vdz2"), free = TRUE, values = 1, newlabels = "vdz1")
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("vdzstw1", "vdzstw2"), free = TRUE, values = 1, newlabels = "vdzstw1")
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("vmzstw1", "vmzstw2"), free = TRUE, values = 1, newlabels = "vmzstw1")

# sibling  means
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("mmzfs", "mdzfs"), free = TRUE, values = 1, newlabels = "mfs1")
# sibling  variances
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("vmzfs", "vdzfs"), free = TRUE, values = 1, newlabels = "vfs1")

# fit and test model fit
esat_mvar_fit <- mxRun(esat_mvar_mod, intervals = F)
mxCompare(esat_fit, esat_mvar_fit)

# extract parameters
esat_mvar_sumy <- summary(esat_mvar_fit)
esat_sumy$parameters

# element to later augment table 1.0
var_mz <- esat_mvar_sumy$parameters %>%
  filter(name == "vmz1") %>%
  pull(Estimate)
var_dz <- esat_mvar_sumy$parameters %>%
  filter(name == "vdz1") %>%
  pull(Estimate)
var_sib <- esat_mvar_sumy$parameters %>%
  filter(name == "vfs1") %>%
  pull(Estimate)
mean_mz <- esat_mvar_sumy$parameters %>%
  filter(name == "mmz1") %>%
  pull(Estimate)
mean_dz <- esat_mvar_sumy$parameters %>%
  filter(name == "mdz1") %>%
  pull(Estimate)
mean_sib <- esat_mvar_sumy$parameters %>%
  filter(name == "mfs1") %>%
  pull(Estimate)
var_mzpartners <- esat_mvar_sumy$parameters %>%
  filter(name == "vmzstw1") %>%
  pull(Estimate)
mean_mzpartners <- esat_mvar_sumy$parameters %>%
  filter(name == "mmzstw1") %>%
  pull(Estimate)
var_dzpartners <- esat_mvar_sumy$parameters %>%
  filter(name == "vdzstw1") %>%
  pull(Estimate)
mean_dzpartners <- esat_mvar_sumy$parameters %>%
  filter(name == "mdzstw1") %>%
  pull(Estimate)

# equate means and variances to test for various differences
# equate mean and variances across zygosities
esat_mvartw1_mod <- mxModel(esat_mvar_mod, name = "mean and var across zyg and sib")
esat_mvartw1_mod <- omxSetParameters(esat_mvartw1_mod, label = c("mmz1", "mdz1", "mfs1"), free = TRUE, values = 1, newlabels = "msib1")
esat_mvartw1_mod <- omxSetParameters(esat_mvartw1_mod, label = c("vmz1", "vdz1", "vfs1"), free = TRUE, values = 1, newlabels = "vsib1")
esat_mvartw1_fit <- mxRun(esat_mvartw1_mod, intervals = F)
mxCompare(esat_fit, esat_mvartw1_fit)

# equate mean and variances across zygosities for spouses of twins
esat_mvastw_mod <- mxModel(esat_mvartw1_mod, name = "mean and var across zyg for partners")
esat_mvastw_mod <- omxSetParameters(esat_mvastw_mod, label = c("mmzstw1", "mdzstw1"), free = TRUE, values = 1, newlabels = "mstw")
esat_mvastw_mod <- omxSetParameters(esat_mvastw_mod, label = c("vmzstw1", "vdzstw1"), free = TRUE, values = 1, newlabels = "vstw")
esat_mvastw_fit <- mxRun(esat_mvastw_mod, intervals = F)
mxCompare(esat_fit, esat_mvastw_fit)

# equate variances
esat_var1_mod <- mxModel(esat_mvastw_mod, name = "var partners and sib")
esat_var1_mod <- omxSetParameters(esat_var1_mod, label = c("vstw", "vsib1"), free = TRUE, values = 1, newlabels = "vp1")
esat_var1_fit <- mxRun(esat_var1_mod, intervals = F)
mxCompare(esat_fit, esat_var1_fit)

# equate means
esat_mean1_mod <- mxModel(esat_var1_mod, name = "mean partners and sib")
esat_mean1_mod <- omxSetParameters(esat_mean1_mod, label = c("mstw", "msib1"), free = TRUE, values = 1, newlabels = "mp1")
esat_mean1_fit <- mxRun(esat_mean1_mod, intervals = F)
mxCompare(esat_fit, esat_mean1_fit)

# extract parameters
esat_mean1_sumy <- summary(esat_mean1_fit)
esat_mean1_sumy$parameters

# sibling correlations
esat_sib_mod <- mxModel(esat_mean1_fit, name = "covar full siblings across zyg")
esat_sib_mod <- omxSetParameters(esat_sib_mod, label = c("cdzfs1", "cdzfs2", "cmzfs1", "cmzfs2"), free = TRUE, values = .1, newlabels = "cfs1")
esat_sib_fit <- mxRun(esat_sib_mod, intervals = F)
mxCompare(esat_fit, esat_sib_fit)

# assortative mating across zygosities
esat_asmtw_mod <- mxModel(esat_sib_mod, name = "cov partners across zyg")
esat_asmtw_mod <- omxSetParameters(esat_asmtw_mod, label = c("cmzstw1", "cmzstw2", "cdzstw1", "cdzstw2"), free = TRUE, values = .1, newlabels = "cprtw")
esat_asmtw_fit <- mxRun(esat_asmtw_mod, intervals = F)
esat_asmtw_fit <- mxTryHard(esat_asmtw_mod, intervals = F)
mxCompare(esat_fit, esat_asmtw_fit)

# in law
esat_law_mod <- mxModel(esat_asmtw_mod, name = "cov in law across zyg")
esat_law_mod <- omxSetParameters(esat_law_mod, label = c("cmzil1", "cmzil2"), free = TRUE, values = .1, newlabels = "cmzinlaw")
esat_law_mod <- omxSetParameters(esat_law_mod, label = c("cdzil1", "cdzil2", "cdzfsil1", "cdzfsil2", "cmzfsil1", "cmzfsil2"), free = TRUE, values = .1, newlabels = "cdzinlaw")
esat_law_fit <- mxRun(esat_law_mod, intervals = F)
esat_law_fit <- mxTryHard(esat_law_mod, intervals = T)
mxCompare(esat_fit, esat_law_fit)
esat_law_sumy <- summary(esat_law_fit)
esat_law_sumy$parameters %>% arrange(name)
esat_law_sumy$CI


# TEST --------------------------------------------------------------
# Sibs
esat_csib_mod <- mxModel(esat_law_mod, name = "test1: full-sibling vs dz")
esat_csib_mod <- omxSetParameters(esat_csib_mod, label = c("cdz", "cfs1"), free = TRUE, values = .1, newlabels = "cfs1")
esat_csib_fit <- mxRun(esat_csib_mod, intervals = F)
esat_csib_fit <- mxTryHard(esat_csib_mod, intervals = F)
mxCompare(esat_law_fit, esat_csib_fit)

# AM (in law)
esat_am_mod <- mxModel(esat_law_mod, name = "test2: mz vs dz in law")
esat_am_mod <- omxSetParameters(esat_am_mod, label = c("cmzinlaw", "cdzinlaw"), free = TRUE, values = .1, newlabels = "cinlaw")
esat_am_fit <- mxRun(esat_am_mod, intervals = F)
esat_am_fit <- mxTryHard(esat_am_mod, intervals = F)
mxCompare(esat_law_fit, esat_am_fit)

# SAVE --------------------------------------------------------------
save(esat_law_fit, file = sprintf("%s/%s/03_saturated_reduced.Rdata", wd_noa, wd_noa_output))

## mx compare list ####
esat_model_comparsion <- mxCompare(esat_fit, c(
  esat_mvar_fit,
  esat_mvartw1_fit,
  esat_mvastw_fit,
  esat_var1_fit,
  esat_mean1_fit,
  esat_sib_fit,
  esat_asmtw_fit,
  esat_law_fit
))

esat_model_comparsion <- esat_model_comparsion %>%
  as.data.frame() %>%
  dplyr::select(-c(fitUnits, fit, diffFit, chisq, SBchisq)) %>%
  mutate(minus2LL = round(minus2LL, 0), AIC = round(AIC, 0), diffLL = round(diffLL, 2), p = round(p, 5))

## supplementary table ####
write_csv(as.data.frame(esat_model_comparsion), sprintf("%s/SI/03_table2_esat_mxcompare.csv", wd_oa))
