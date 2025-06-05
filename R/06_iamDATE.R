# Author: Giacomo Bignardi
# iam-DATE model
# Test for dAM and estimate paths to variance components accounting for AM
library(tidyverse)
require(OpenMx)
library(patchwork)
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

# load extended twin family data
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

# load saturated model for later model comparison
load(sprintf("%s/%s/03_saturated.Rdata", wd_noa, wd_noa_output))

# source ctd and estd inspired model
source(sprintf("%s/%s/function/iam.date.fun.R", wd_oa, wd_oa_scripts))
source(sprintf("%s/%s/function/acde.fun.R", wd_oa, wd_oa_scripts))

# select family members for analysis and adapt to specifics of iAM model
# remember that the order needs to match exactly the order in the function
vars <- "aes"
sel_vars <- paste0("aes", "_", c("stw1", "tw1", "s1", "tw2", "stw2"))
# select data for analysis
mz <- subset(dat_fam2g, twzyg == "MZ", sel_vars)
dz <- subset(dat_fam2g, twzyg == "DZ", sel_vars)
mz_sample <- psych::pairwiseCount(mz)
dz_sample <- psych::pairwiseCount(dz)

# name for readability
colnames(mz_sample) <- sel_vars
rownames(mz_sample) <- sel_vars
colnames(dz_sample) <- sel_vars
rownames(dz_sample) <- sel_vars

# ADE --------------------------------------------------------------
# first fit classical twin design
# with dominance
ade_mod <- acde.fun(mz, dz, c("aes_tw1", "aes_tw2"), dp.free = T)
ade_fit <- mxRun(ade_mod, intervals = T)
ade_fit$US$result

# without dominance
ae_mod <- acde.fun(mz, dz, c("aes_tw1", "aes_tw2"))
ae_fit <- mxRun(ae_mod, intervals = T)
ae_fit$US$result

# compare significance of the D parameter
mxCompare(ade_fit, ae_fit)

# iam_date --------------------------------------------------------------
# then fit extended twin family design
# fit the model with indirect assortment and dominance
iam_date_mod <- iam.date.fun(mz,
  dz,
  sel_vars,
  indirect.assortment = T,
  # focal paths
  d1.free = T,
  n1.free = F,
  # sorting path
  d1s.free = T,
  n1s.free = F,
  # starting values
  n1.val = 0,
  # name of the model
  name = "iam_date"
)
iam_date_fit <- mxRun(iam_date_mod, intervals = T)
iam_date_fit$est_var
iam_date_fit$est_path

# fit the model with direct assortment and dominance
dam_date_mod <- iam.date.fun(mz,
  dz,
  sel_vars,
  indirect.assortment = F,
  # focal paths
  d1.free = T,
  n1.free = F,
  # sorting path
  d1s.free = T,
  n1s.free = F,
  # starting values
  n1.val = 0,
  # name of the model
  name = "dam_date"
)
dam_date_fit <- mxRun(dam_date_mod, intervals = T)
dam_date_fit$est_var
dam_date_fit$est_path

# test if indirect assortment can be constrained
mxCompare(iam_date_fit, dam_date_fit) # indirect assortment can be constrained

# fit the model with direct assortment and dominance but not twin-specific effects
dam_ade_mod <- iam.date.fun(mz,
  dz,
  sel_vars,
  indirect.assortment = F,
  # focal paths
  d1.free = T,
  n1.free = F,
  t1.free = F,
  # sorting path
  d1s.free = T,
  n1s.free = F,
  t1s.free = F,
  # starting values
  n1.val = .0,
  t1.val = .0,
  # name of the model
  name = "dam_ade"
)
dam_ade_fit <- mxRun(dam_ade_mod, intervals = T)
dam_ade_fit <- mxTryHard(dam_ade_mod, intervals = T)
dam_ade_fit$est_var

# compare models
mxCompare(dam_date_fit, dam_ade_fit) # twin effects can be dropped

# test if dominance effects is significant
dam_ae_mod <- iam.date.fun(mz,
  dz,
  sel_vars,
  indirect.assortment = F,
  # focal paths
  d1.free = F,
  n1.free = F,
  t1.free = F,
  # sorting path
  d1s.free = F,
  n1s.free = F,
  t1s.free = F,
  # starting values
  n1.val = .0,
  t1.val = .0,
  # name of the model
  name = "dam_ae"
)
dam_ae_fit <- mxRun(dam_ae_mod, intervals = T)
dam_ae_fit <- mxTryHard(dam_ae_mod, intervals = T)

# Dominance cannot be dropped
mxCompare(dam_ade_fit, dam_ae_fit) # dominance cannot be dropped

# final model
iam_date_mxcompare <- rbind(
  as.data.frame(mxCompare(esat_fit, iam_date_fit)),
  as.data.frame(mxCompare(iam_date_fit, dam_date_fit))[2, ],
  as.data.frame(mxCompare(dam_date_fit, dam_ade_fit))[2, ],
  as.data.frame(mxCompare(dam_ade_fit, dam_ae_fit))[2, ]
)

# comparison with final model
mxCompare(esat_fit, c(dam_ade_fit))
iam_date_mxcompare <- iam_date_mxcompare %>%
  dplyr::select(-c(fitUnits, fit, diffFit, chisq, SBchisq)) %>%
  mutate(minus2LL = round(minus2LL, 0), AIC = round(AIC, 0), diffLL = round(diffLL, 2), p = round(p, 5))

# save for supplementary table
write_csv(iam_date_mxcompare, sprintf("%s/SI/06_table3_mod_comp.csv", wd_oa))

# iam_nate --------------------------------------------------------------
# repeat using  boundary scenarios of rNA = 0(conservative estimates for significance testing)
# fit the model with indirect assortment and non-additive
iam_nate_mod <- iam.date.fun(mz,
  dz,
  sel_vars,
  indirect.assortment = T,
  # focal paths
  d1.free = F,
  n1.free = T,
  # sorting path
  d1s.free = F,
  n1s.free = T,
  # starting values
  n1.val = .2,
  # name of the model
  name = "iam_nate"
)
iam_nate_fit <- mxRun(iam_nate_mod, intervals = T)
iam_nate_fit <- mxTryHard(iam_nate_mod, intervals = T)

# fit the model with direct assortment and non-additive
dam_nate_mod <- iam.date.fun(mz,
  dz,
  sel_vars,
  indirect.assortment = F,
  # focal paths
  d1.free = F,
  n1.free = T,
  # sorting path
  d1s.free = F,
  n1s.free = T,
  # starting values
  n1.val = 0,
  # name of the model
  name = "amNATE"
)
dam_nate_fit <- mxRun(dam_nate_mod, intervals = T)

# test if indirect assortment can be constrained
mxCompare(iam_nate_fit, dam_nate_fit) # indirect assortment can be constrained

# fit the model with direct assortment and non-additive but not twin-specific effects
dam_ane_mod <- iam.date.fun(mz,
  dz,
  sel_vars,
  indirect.assortment = F,
  # focal paths
  d1.free = F,
  n1.free = T,
  t1.free = F,
  # sorting path
  d1s.free = F,
  n1s.free = T,
  t1s.free = F,
  # starting values
  n1.val = .0,
  t1.val = .0,
  # name of the model
  name = "dam_ane"
)
dam_ane_fit <- mxRun(dam_ane_mod, intervals = T)
dam_ane_fit <- mxTryHard(dam_ane_mod, intervals = T)

dam_ane_fit$est_var
# compare models
mxCompare(dam_nate_fit, dam_ane_fit) # twin effects can be dropped
# non-additive cannot be dropped
mxCompare(dam_ane_fit, dam_ae_fit)

# final model
mxCompare(esat_fit, c(iam_nate_fit, dam_ane_fit))

# iam_nate with 0<=rNA<=.25 --------------------------------------------------------------
# estimate variety of scenarios
# fit the model with indirect assortment and non-additive but not twin-specific effects
iam_nate_est <- c()
for (i in seq(0, .24, .01)) {
  iam_nate_mod <- iam.date.fun(mz,
    dz,
    sel_vars,
    indirect.assortment = T,
    # focal paths
    d1.free = F,
    n1.free = T,
    t1.free = T,
    # sorting path
    d1s.free = F,
    n1s.free = T,
    t1s.free = T,
    # starting values
    n1.val = .0,
    # name of the model
    name = "iam_nate",
    # change degree of r epistatic
    rI.f = i
  )
  iam_nate_fit <- mxTryHard(iam_nate_mod, intervals = T)
  iam_nate_est_i <- summary(iam_nate_fit)$CI %>%
    as.data.frame() %>%
    mutate(parameter = c(colnames(iam_nate_fit$est_var$result), colnames(iam_nate_fit$est_path$result)), model = "iam_nate", ri = i)
  iam_nate_est <- rbind(iam_nate_est, iam_nate_est_i)
}

# fit the model with direct assortment and dominance but not twin-specific effects
dam_ane_est <- c()
for (i in seq(0, .24, .01)) {
  dam_ane_mod <- iam.date.fun(mz,
    dz,
    sel_vars,
    indirect.assortment = F,
    # focal paths
    d1.free = F,
    n1.free = T,
    t1.free = F,
    # sorting path
    d1s.free = F,
    n1s.free = T,
    t1s.free = F,
    # starting values
    n1.val = .0,
    t1.val = .0,
    # name of the model
    name = "dam_ane",
    # change degree of r epistatic
    rI.f = i
  )
  dam_ane_fit <- mxTryHard(dam_ane_mod, intervals = T)
  dam_ane_est_i <- summary(dam_ane_fit)$CI %>%
    as.data.frame() %>%
    mutate(parameter = c(colnames(dam_ane_fit$est_var$result), colnames(dam_ane_fit$est_path$result)), model = "dam_ane", ri = i)
  dam_ane_est <- rbind(dam_ane_est, dam_ane_est_i)
}

# SUMMARY####
# create summary of estimates for models with direct and indirect assortment
ade_sumy <- summary(ade_fit)
ade_est <- ade_sumy$CI %>%
  as.data.frame() %>%
  mutate(parameter = colnames(ade_fit$US$result), model = "ADE")
ae_sumy <- summary(ae_fit)
ae_est <- ae_sumy$CI %>%
  as.data.frame() %>%
  mutate(parameter = colnames(ae_fit$US$result), model = "AE")
iam_date_sumy <- summary(iam_date_fit)
iam_date_est <- iam_date_sumy$CI %>%
  as.data.frame() %>%
  mutate(parameter = c(colnames(iam_date_fit$est_var$result), colnames(iam_date_fit$est_path$result)), model = "iam_date")
dam_ade_sumy <- summary(dam_ade_fit)
dam_ade_est <- dam_ade_sumy$CI %>%
  as.data.frame() %>%
  mutate(parameter = c(colnames(dam_ade_fit$est_var$result), colnames(dam_ade_fit$est_path$result)), model = "dam_ade")
est <- rbind(ade_est, ae_est, iam_date_est, dam_ade_est)

# Save table 4
write_csv(est %>% filter(!parameter %in% c("varP_N", "varS_N", "n1s", "n1", "varP_C", "varS_C", "c1s", "c1")) %>% select(model, parameter, estimate, lbound, ubound) %>% mutate(lbound = round(lbound, 5), estimate = round(estimate, 5), ubound = round(ubound, 5)), sprintf("%s/SI/06_table4_est.csv", wd_oa))

# summary for non-additive effects
nate_est <- rbind(iam_nate_est, dam_ane_est)

# comparison
mxCompare(esat_fit, iam_date_fit)
mxCompare(esat_fit, iam_nate_fit)

# save for later
nate_est_tb <- nate_est %>%
  filter(parameter %in% c("varP_A", "varP_N", "varP_E")) %>%
  filter(model == "dam_ane") %>%
  mutate(
    parameter = case_match(
      parameter,
      "varP_E" ~ "E",
      "varP_N" ~ "N",
      "varP_A" ~ "A"
    ),
    parameter = factor(parameter, levels = c("A", "N", "E"))
  ) %>%
  rename(rNA = ri) %>%
  select(model, rNA, parameter, estimate, lbound, ubound) %>%
  arrange(rNA) %>%
  mutate(lbound = round(lbound, 5), estimate = round(estimate, 5), ubound = round(ubound, 5))

# Save table 5
write_csv(nate_est_tb, sprintf("%s/SI/06_table5_est.csv", wd_oa))

# FIGURES
# note this figures are not included in the article but
# since I have made them for a previous draft I left the code
# here in case anyone (including myself) would like to plot them
# in the future
# AE
# plot_a =  nate_est %>%
#   filter(parameter %in% c("varP_A","varP_N")) %>%
#   filter(ri == .0) %>%
#   filter(model %in% c("dam_ane")) %>%
#   mutate(model = factor(model, levels = c("iam_nate","dam_ane"))) %>%
#   mutate(parameter_est = case_match(  parameter,
#                                       "varP_A" ~ "Additive genetic",
#                                       "varP_N" ~ "Interaction deviation"),
#          parameter = case_match(  parameter,
#                                   "varP_A" ~ "A",
#                                   "varP_N" ~ "N"),
#          parameter_est = factor(parameter_est, levels = rev(c("Interaction deviation","Residual environmental", "Twin specific","Dominance deviation","Additive genetic"))),
#          parameter = factor(parameter, levels =c("A","N","T")))%>%
#   ggplot(aes(parameter_est,estimate, fill = parameter)) +
#   geom_bar(stat = "identity", position= position_dodge(), color = "black", width = .5) +
#   geom_errorbar(aes(ymin = lbound, ymax = ubound), position = position_dodge(), width = .1)+
#   scale_x_discrete(labels = c("Additive genetic" = expression(hat(italic(h))[twin-ped]^{"2"}),
#                               "Interaction deviation" = expression(hat(delta)[twin-ped]^{"2"}),
#                               "Twin specific" = expression(hat(tau)[twin-ped]^{"2"})))+
#   # geom_text(aes( y = .025, label = round(estimate,2)), vjust = 0, position = position_dodge(.8),color = "white") +
#   scale_fill_manual(values = c("#B41F01","#FAA53C","#5c53a5"))+ #FAA53C
#   labs(
#     fill = "parameters (epistasis)",
#     x = "Estimates",
#     y= "")+
#   scale_y_continuous(limits = c(0,.35))+
#   theme_classic(base_size = 12)+
#   theme(strip.background = element_blank(),
#         legend.position = "none")
#
# plot_c =  est %>%
#   filter(parameter %in% c("varP_A","varP_D")) %>%
#   filter(model %in% c("dam_ade")) %>%
#   mutate(parameter = factor(parameter, levels = c("varP_D","varP_A","varP_N","varP_T","varP_E"))) %>%
#   ggplot(aes(parameter,estimate, fill = parameter)) +
#   geom_bar(stat = "identity", position= position_dodge(),color = "black", width = .5) +
#   geom_errorbar(aes(ymin = lbound, ymax = ubound), position = position_dodge(), width = .3)+
#   scale_x_discrete(labels = c("varP_A" = expression(hat(italic(h))[twin-ped]^{"2"}),
#                               "varP_D" = expression(hat(delta)[twin-ped]^{"2"}),
#                               "varP_N" = expression(hat(iota)[twin]^{"2"}),
#                               "varP_T" = expression(hat(tau)[twin]^{"2"}),
#                               "varP_E" = expression(hat(epsilon)[twin]^{"2"})))+
#   # geom_text(aes( y = .025, label = round(estimate,2)), vjust = 0, position = position_dodge(.8),color = "white") +
#   scale_fill_manual(values = c("#FAA53C","#B41F01","#3E59DC","#041474"),
#                     labels = c("varP_A" = "A",
#                                "varP_D" = "D",
#                                "varP_N" = "N",
#                                "varP_T" = "T",
#                                "varP_E" = "E"))+ #FAA53C
#   labs(    fill = "Parameters",
#            x = "Estimates",
#            y= "")+
#   theme_classic(base_size = 12)+
#   ylim(c(0,.35))+
#   # scale_y_continuous(limits = c(0,.4))+
#   theme(strip.background = element_blank(),
#         legend.position = "none")
#
# plot_b <- nate_est %>%
#   filter(parameter %in% c("varP_A","varP_N")) %>%
#   filter(model ==  "dam_ane") %>%
#   mutate(parameter_est = case_match(  parameter,
#                                       "SE1" ~ "Residual",
#                                       "ST1" ~ "Twin specific",
#                                       "SD1" ~ "Dominance deviation",
#                                       "SD"  ~ "Dominance deviation",
#                                       "varP_N" ~ "Interaction deviation",
#                                       "SI"  ~ "Interaction deviation",
#                                       "varP_A" ~ "Additive genetic",
#                                       "SA"  ~ "Additive genetic"),
#          parameter = case_match(  parameter,
#                                   "SE1" ~ "E",
#                                   "ST1" ~ "T",
#                                   "SD1" ~ "D",
#                                   "SD"  ~ "D",
#                                   "varP_N" ~ "N",
#                                   "SI"  ~ "N",
#                                   "varP_A" ~ "A",
#                                   "SA"  ~ "A"),
#          model = case_match(  model,
#                               "amAIE" ~ "interaction deviation"),
#          parameter_est = factor(parameter_est, levels = rev(c("Interaction deviation","Residual environmental", "Twin specific","Dominance deviation","Additive genetic"))),
#          parameter = factor(parameter, levels =c("A","N","T"))) %>%
#   ggplot(aes(ri,estimate, fill = parameter, color = parameter)) +
#   geom_ribbon(aes(ymin = lbound, ymax = ubound, fill = parameter), size = 1, shape=21, alpha = .75)+
#   geom_line(size = 1.5)+
#   geom_vline(xintercept = .25, linetype = "dashed")+
#   geom_vline(xintercept = .00, linetype = "dashed")+
#   # geom_text(aes( y = .025, label = round(estimate,2)), vjust = 0, position = position_dodge(.8),color = "white") +
#   scale_fill_manual(values = c("#B41F01","#FAA53C"),##e34f6f
#                     labels = c("A" = expression(hat(italic(h))[twin-ped]^{"2"}),
#                                "D" = expression(hat(delta(h))[twin-ped]^{"2"}),
#                                "N" = expression(hat(delta)[twin]^{"2"})))+ #FAA53C
#   scale_color_manual(values = c("#B41F01","#FAA53C"),guide = "none")+
#   labs(
#     x = expression(paste(italic(r)["NA"])),
#     y = "Values",
#     fill = "Estimates"
#     )+
#   theme_classic(base_size = 12)+
#   scale_y_continuous(limits = c(0,.35))+
#   theme(strip.background = element_blank())
# pdf(file = ?
#     width = 8,
#     height = 2.5)
# (plot_a|plot_b|plot_c) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")
# dev.off()