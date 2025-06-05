# Author: Giacomo Bignardi
# Correlations between family members
# library(ggplot2)
# library(patchwork)
library(tidyverse)
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

# load estimates from the constrained sat models
load(file = sprintf("%s/%s/03_saturated_reduced.Rdata", wd_noa, wd_noa_output))

# load number of complete pairs
n_compair <- read.csv(sprintf("%s/%s/03_pair_number.csv", wd_noa, wd_noa_output))
esat_law_sumy <- summary(esat_law_fit)
fam_cov <- esat_law_sumy$parameters %>%
  filter(grepl("^[c]", name)) %>%
  dplyr::select(name, Estimate, Std.Error) %>%
  mutate(name = recode(name,
    "cprtw" = "am",
    "cmz" = "mz",
    "cdz" = "dz",
    "cfs1" = "sib",
    "cmzinlaw" = "mz-in-law",
    "cdzinlaw" = "sib-in-law",
    "cmzstwcil" = "mz-co-in-law",
    "cdzstwcil" = "sib-co-in-law",
  ))

# rename for later plotting and ease of interpretability
fam_cor <- esat_law_sumy$CI %>%
  rownames_to_column() %>%
  mutate(type = recode(rowname,

    # Cor partners
    "mz.cormz[2,1]" = "ram",
    # Cor MZ
    "mz.cormz[4,2]" = "rmz",
    # Cor DZ
    "dz.cordz[4,2]" = "rdz",
    # Cor Sib
    "mz.cormz[3,2]" = "rsib",
    # Cor MZ in law
    "mz.cormz[4,1]" = "rmz-in-law",
    # Cor Sibs in-law
    "dz.cordz[4,1]" = "rsib-in-law",
    # Cor Sibs in-colaw
    "mz.cormz[5,1]" = "rmz-co-in-law",
    # Cor Sibs in-colaw
    "dz.cordz[5,1]" = "rsib-co-in-law"
  )) %>%
  # dplyr::select only correlations
  filter(type %in% c(
    "ram",
    "rmz",
    "rdz",
    "rsib",
    "rmz-in-law",
    "rsib-in-law",
    "rmz-co-in-law",
    "rsib-co-in-law"
  ))

# total number of sib il law pairs (including dz)
n_compair$nsil_comp_pair <- n_compair$ndzil_comp_pair + n_compair$nfsil_comp_pair

# rename sample size to match correlations
n_compair_long <- n_compair %>%
  rename(
    "ram" = nstw_comp_pair,
    "rmz" = nmz_comp_pair,
    "rdz" = ndz_comp_pair,
    "rsib" = nfs_comp_pair,
    "rmz-in-law" = nmzil_comp_pair,
    "rsib-in-law" = nsil_comp_pair,
    "rmz-co-in-law" = nmzcil_comp_pair,
    "rsib-co-in-law" = ndzcil_comp_pair
  ) %>%
  dplyr::select(-c(
    "ndzil_comp_pair",
    "nfsil_comp_pair"
  )) %>%
  pivot_longer(names_to = "type", values_to = "n", rmz:"rsib-in-law")

# merge correlations and sample
fam_cor_fin <- rbind(merge(fam_cor, n_compair_long, by = c("type"), all = T))
# create a column to merge with covariances
fam_cor_fin$name <- substr(fam_cor_fin$type, 2, nchar(fam_cor_fin$type))
fam_cor_fin$type <- factor(fam_cor_fin$type, levels = rev(c(
  "rmz", "rdz", "rsib", "ram", "rmz-in-law", "rsib-in-law", "rmz-co-in-law", "rsib-co-in-law"
)))

# prepare table with cov cor estimates
fam_est <- merge(fam_cov, fam_cor_fin)
fam_est$name <- sapply(strsplit(fam_cor_fin$name, "\\."), `[`, 1)
fam_est <- fam_est %>%
  as.data.frame() %>%
  arrange(desc(type)) %>%
  dplyr::select(name, type, Estimate, Std.Error, estimate, lbound, ubound) %>%
  rename("family relationship" = name, cov = Estimate, SE = Std.Error, r = estimate) %>%
  mutate(cov = round(cov, 5), SE = round(SE, 5), r = round(r, 5), lbound = round(lbound, 5), ubound = round(ubound, 5))

# Figure #####
# This is an additonal figure (not present in the article)
# it plots the correlations in fam_est (table 2)
# fam_corplot_ntr <- fam_cor_fin %>%
#   filter(!type %in% c("rmz-co-in-law","rsib-co-in-law") ) %>%
#   ggplot(aes(type, y = estimate)) +
#   geom_linerange(aes(ymin = lbound, ymax = ubound),color = "#06175F", size = .75,position = position_dodge(.9)) +
#   geom_point(fill = "#06175F", size = 2, shape=21) +
#   ylim(c(-.1,.4))+
#   geom_text(aes(label = paste0("n=",n), y = -.075),
#             nudge_x = -0.2,
#             hjust = 0,
#             size = 3.5)+
#   # scale_fill_manual(values = c("#fea24f"))+
#   scale_x_discrete(labels = c(
#     "ram"= "partners",
#     "rmz"= "MZ twins",
#     "rdz"= "DZ twins",
#     "rsib"= "siblings",
#     "rmz-in-law" = "MZ in-law",
#     "rsib-in-law" = "siblings in-law" ,
#     "rtw-co-in-law" = "siblings co-in-law")) +
#   labs(x = expression(italic(r)))+
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_classic(base_size = 12)+
#   labs(y =  expression(paste(italic(r), " aesthetic chills")),
#        x =  "Familial relationship") +
#   theme(strip.text.x = element_text(angle = 0),
#         legend.title.align=0.5)
