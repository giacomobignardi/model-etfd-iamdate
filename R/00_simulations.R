# SCRIPT: simulations
# Author: G.B. (Giaco)
# part of the multivariate model to simulate twin data is based on a script
# shared by S.B. (originally written by C.V.D)
rm(list = ls())
library(tidyverse)
library(OpenMx) # to fit SEM model
library(MASS)   # to simulate data

# Optimizer
# use NPSOL (non linear constrains)
# ?omxGetNPSOL()
mxOption(NULL, "Default optimizer", "NPSOL")
set.seed(42)

# set open access working directories
wd_oa = getwd()
wd_oa_scripts = "R"
wd_oa_figures = "Figures"

# set other working directories
wd_noa = substr(
  getwd(),
  0,
  nchar(getwd())-nchar("11_github")-1
)
wd_noa_Data = "03_data"
wd_noa_output = "03_outputs"

# FUNCTIONS --------------------------------------------------------------
# source ctd-based model (classical twin design)
source(sprintf("%s/%s/function/acde.fun.R", wd_oa,wd_oa_scripts))

# source etfd-based model (extended twin family design)
source(sprintf("%s/%s/function/iam.date.fun.R", wd_oa,wd_oa_scripts))

# SIMULATION --------------------------------------------------------------
## starting values ####
# select n of families (here 5k with MZ and 5k with DZ twins)
Nmz <- 5000
Ndz <- 5000

# as raw ML implies division by N not (N-1)
cor_Nmz <- (Nmz - 1) / Nmz
cor_Ndz <- (Ndz - 1) / Ndz
N <- Nmz + Ndz
Nmem <- 5 # 2 partners + 2 twins + 1 full-sib
Nmem1 <- Nmem + 1
exact <- TRUE # if TRUE the parameter recovery should be exact (sample parameter, not estimate from population sample)

# value of assortment
mu.vals <- seq(0, .9, .1) # mu co-path for assortative mating
est_am <- c()

# set squared paths coefficients for non-additive effects
d2 <- .15 # dominance genetics
i2 <- .00 # other non-additive genetic effects beyond dominance - here interaction deviation
n2 <- d2 + i2 # non-additive genetics
t2 <- .00 # social hom.

# set incremental for additive effects
a.vals <- seq(2 * (n2 + t2), 1 - (.1 + n2 + t2), .05)

# type of assortment
am.type <- c(
  "dAM",
  "narrow-sense hom.",
  "broad-sense hom.",
  "social hom.",
  "idio. hom.",
  "mixed hom."
)

## starting values for AM ####
# loop across all possible scenarios
# loop 1: across type of homogamy
# loop 2: across values of h2
# loop 3: across values of mu
for (am in am.type) {
  est_h2 <- c()
  for (a2 in a.vals) {
    est_mu <- c()
    for (mu in mu.vals) {
      mu.val <- mu # mu co-path for assortative mating
      # focal phenotype
      a1 <- sqrt(a2) # additive genetics
      c1 <- log(1.0) # shared environment (always assumed to be equal to 0)
      t1 <- sqrt(t2) # twin environment
      d1 <- sqrt(d2) # dominance genetics
      i1 <- sqrt(i2) # interaction deviation
      e1 <- sqrt(1 - a2 - d2 - i2 - t2) # unique environment

      if (am == "dAM") {
        # sorting phenotype
        a1s <- sqrt(a2) # narrow-sense hom.
        c1s <- log(1.0) # social hom. (shared; (always assumed to be equal to 0))
        t1s <- sqrt(t2) # social hom. (twin-shared)
        d1s <- sqrt(d2)
        is1 <- sqrt(i2)
        e1s <- sqrt(1 - a2 - d2 - i2 - t2) # idyosincratic hom.
      } else if (am == "narrow-sense hom.") {
        # sorting phenotype (narrow-sense hom.)
        a1s <- sqrt(1) # purely correlate via A
        c1s <- log(1.0)
        t1s <- log(1.0)
        d1s <- log(1.0)
        is1 <- log(1.0)
        e1s <- log(1.0)
      } else if (am == "broad-sense hom.") {
        # sorting phenotype (broad-sense hom.)
        a1s <- sqrt(2 / 3) # correlate via A
        c1s <- log(1.0)
        t1s <- log(1.0)
        d1s <- sqrt(1 / 3) # correlate via D
        is1 <- log(1.0)
        e1s <- log(1.0)
      } else if (am == "social hom.") {
        # s orting phenotype (twin-shared hom.)
        a1s <- log(1.0)
        c1s <- log(1.0)
        t1s <- sqrt(1) # purely correlate via T
        d1s <- log(1.0)
        is1 <- log(1.0)
        e1s <- log(1.0)
      } else if (am == "mixed hom.") {
        # sorting phenotype (mixed hom.)
        a1s <- sqrt(2 / 6) # correlate via A
        c1s <- log(1.0)
        t1s <- sqrt(.25) # correlate via T
        d1s <- sqrt(1 / 6) # correlate via D
        is1 <- log(1.0)
        e1s <- sqrt(.25) # correlate via E
      } else if (am == "idio. hom.") {
        # sorting phenotype (idio. hom. homogamy)
        a1s <- log(1.0)
        c1s <- log(1.0)
        t1s <- log(1.0)
        d1s <- log(1.0)
        is1 <- log(1.0)
        e1s <- sqrt(1) # purely correlate via E
      }

      ## generative model ####
      # bulid variances and copahts
      Vp <- a1^2 + c1^2 + t1^2 + d1^2 + i1^2 + e1^2 # phenotypic variance focal phenotype (trait)
      Vs <- a1s^2 + c1s^2 + t1s^2 + d1s^2 + is1^2 + e1s^2 # phenotypic variance sorting phenotype (sorting factor)
      mu <- mu.val / Vs^2

      # partner genetic correlations
      partner_rG <- mu * a1s^2

      # DZ and sib correlations
      fs_rA <- (1 + partner_rG) / 2 # A correlation
      fs_rD <- .25 # D correlation
      fs_rN <- .0 # N high-order epistasis boundary

      # build the MZ and DZ covariance matrices based on the parameter values given above
      Vp <- a1^2 + c1^2 + t1^2 + d1^2 + i1^2 + e1^2 # phenotypic variance focal phenotype
      Vs <- a1s^2 + c1s^2 + t1s^2 + d1s^2 + is1^2 + e1s^2 # phenotypic variance sorting phenotype
      mu <- mu.val / Vs^2
      Cps <- a1 * a1s + c1 * c1s + t1 * t1s + d1 * d1s + i1 * is1 + e1 * e1s # covariance between focal phenotype and sorting factor
      Cpar <- Cps * mu * Cps # covariance between partners
      Cmz <- 1 * a1^2 + c1^2 + 1 * d1^2 + 1 * i1^2 + t1^2 # mz covariance
      Cdz <- fs_rA * a1^2 + c1^2 + fs_rD * d1^2 + fs_rN * i1^2 + t1^2 # dz covariance
      Cfs <- fs_rA * a1^2 + c1^2 + fs_rD * d1^2 + fs_rN * i1^2 # fs covariance
      Cmzil <- mu * Cps * (1 * a1 * a1s + c1 * c1s + 1 * d1 * d1s + 1 * i1 * is1 + t1 * t1s) # mz in-law
      Cdzil <- mu * Cps * (fs_rA * a1 * a1s + c1 * c1s + fs_rD * d1 * d1s + fs_rN * i1 * is1 + t1 * t1s) # dz in-law
      Cfsil <- mu * Cps * (fs_rA * a1 * a1s + c1 * c1s + fs_rD * d1 * d1s + fs_rN * i1 * is1) # fs in-law
      Cmzcil <- Cps^2 * (mu^2 * (1 * a1s^2 + c1s^2 + 1 * d1s^2 + 1 * is1^2 + t1s^2)) # mz co-in-law
      Cdzcil <- Cps^2 * (mu^2 * (fs_rA * a1s^2 + c1s^2 + fs_rD * d1s^2 + fs_rN * is1^2 + t1s^2)) # dz co-in-law

      # variance components
      # focal phenotype
      SA1 <- a1^2 / Vp
      SC1 <- c1^2 / Vp
      ST1 <- t1^2 / Vp
      SD1 <- d1^2 / Vp
      SI1 <- i1^2 / Vp
      SE1 <- e1^2 / Vp

      # MZ covariance matrix
      # p1, tw1, s1, tw2, p2
      Smz <- matrix(
        c(
          Vp,     Cpar,  Cfsil, Cmzil, Cmzcil,
          Cpar,   Vp,    Cfs,   Cmz,   Cmzil,
          Cfsil,  Cfs,   Vp,    Cfs,   Cfsil,
          Cmzil,  Cmz,   Cfs,   Vp,    Cpar,
          Cmzcil, Cmzil, Cfsil, Cpar,  Vp
        ),
        5, 5,
        byrow = T
      )

      # DZ covariance matrix
      Sdz <- matrix(
        c(
          Vp,     Cpar,  Cfsil, Cdzil, Cdzcil,
          Cpar,   Vp,    Cfs,   Cdz,   Cdzil,
          Cfsil,  Cfs,   Vp,    Cfs,   Cfsil,
          Cdzil,  Cdz,   Cfs,   Vp,    Cpar,
          Cdzcil, Cdzil, Cfsil, Cpar,  Vp
        ),
        5, 5,
        byrow = T
      )

      # phenotypic means
      me <- rep(0, 5)

      # true parameter
      true <- data.frame(
        est = c(SA1, SC1, ST1, SD1, SI1, SE1),
        component = c("A", "C", "T", "D", "N", "E")
      ) # note the difference due to the increase in phenotype variance given q as a consequence of AM

      # iterate simulation and estimate parameters
      est_ctd <- c()
      sel_vars <- c("stw1", "tw1", "s1", "tw2", "stw2")
      sel_vars_tw <- c("tw1", "tw2")

      for (i in 1) {
        ## exact data simulation ####
        zyg <- c(rep(1, Nmz), rep(2, Ndz))
        dat <- matrix(0, N, Nmem1)
        dat[, 1] <- zyg
        dat[zyg == 1, 2:Nmem1] <- mvrnorm(Nmz, mu = me, Sigma = Smz / cor_Nmz, emp = exact) # note here is sample parameters
        dat[zyg == 2, 2:Nmem1] <- mvrnorm(Ndz, mu = me, Sigma = Sdz / cor_Ndz, emp = exact) # note here is sample parameters
        colnames(dat) <- c("zyg", sel_vars)
        dat <- as.data.frame(dat)

        # data files
        mz <- dat[dat$zyg == 1, c("zyg", sel_vars_tw)][, -1]
        dz <- dat[dat$zyg == 2, c("zyg", sel_vars_tw)][, -1]
        mzfull <- dat[dat$zyg == 1, c("zyg", sel_vars)][, -1]
        dzfull <- dat[dat$zyg == 2, c("zyg", sel_vars)][, -1]

        # extract correlations
        rmz <- cor(mz)["tw1", "tw2"]
        rdz <- cor(dz)["tw1", "tw2"]
        rsib <- cor(mzfull)["tw1", "s1"]
        rmzil <- cor(mzfull)["tw1", "stw2"]
        rdzil <- cor(dzfull)["tw1", "stw2"]
        rsibil <- cor(mzfull)["s1", "stw2"]
        rmzcil <- cor(mzfull)["stw1", "stw2"]
        rdzcil <- cor(dzfull)["stw1", "stw2"]

        # ACE 2*rDZ>rMZ otherwise ADE
        model_ctd <- c()
        acde_mod <- if ((2 * rdz - rmz) > 0) {
          acde.fun(mz, dz, sel_vars_tw, cp.free = T)
        } else {
          acde.fun(mz, dz, sel_vars_tw, dp.free = T)
        }
        acde_fit <- mxRun(acde_mod, intervals = F)

        # extract results to compare
        est_ctd_i <- c()
        est_ctd_i <- acde_fit$US$result[, c("varP_A", "varP_C", "varP_D", "varP_E")] %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          rename(component = rowname, est = ".") %>%
          mutate(model = "ADCE")
        est_ctd_i <- rbind(est_ctd_i, data.frame(component = c("MZ", "DZ", "Full sib", "MZ in-law", "DZ in-law", "Full sib in-law", "MZ co-in-law", "Full sib co-in-law"), est = c(rmz, rdz, rsib, rmzil, rdzil, rsibil, rmzcil, rdzcil), model = "ctd"))
        est_ctd <- rbind(est_ctd, est_ctd_i)
      }

      # extended twin family design (iAM_DATE)
      est_estd <- c()
      for (i in 1) {
        # sampling data simulation ####  -----------------------
        zyg <- c(rep(1, Nmz), rep(2, Ndz))
        dat <- matrix(0, N, Nmem1)
        dat[, 1] <- zyg
        dat[zyg == 1, 2:Nmem1] <- mvrnorm(Nmz, mu = me, Sigma = Smz / cor_Nmz, emp = FALSE) # note here is population parameters
        dat[zyg == 2, 2:Nmem1] <- mvrnorm(Ndz, mu = me, Sigma = Sdz / cor_Ndz, emp = FALSE) # note here is population parameters
        colnames(dat) <- c("zyg", sel_vars)
        dat <- as.data.frame(dat)

        # data files
        mz <- dat[dat$zyg == 1, c("zyg", sel_vars)][, -1]
        dz <- dat[dat$zyg == 2, c("zyg", sel_vars)][, -1]

        # fit different NFTD models
        iAM_DATE_mod <- if (am == "dAM") {
          iam.date.fun(mz, dz, sel_vars, indirect.assortment = F, d1.free = T, n1.free = F)
        } else {
          iam.date.fun(mz, dz, sel_vars, indirect.assortment = T, d1.free = T, d1s.free = T, n1.free = F, n1s.free = F, e1.val = .5)
        }
        iAM_DATE_fit <- mxTryHard(iAM_DATE_mod, intervals = F)

        # extract results to compare
        est_estd_i <- c()
        est_estd_i <- iAM_DATE_fit$est_var$result[, c("varP_A", "varP_C", "varP_T", "varP_D", "varP_N", "varP_E")] %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          rename(component = rowname, est = ".") %>%
          mutate(model = "iamDATE")
        est_estd <- rbind(est_estd, est_estd_i)
      }
      # across values of mu
      est_mu <- rbind(est_mu, rbind(est_ctd %>% mutate(mu = mu), est_estd %>% mutate(mu = mu)))
    }
    # across values of h2
    est_h2 <- rbind(est_h2, est_mu %>% mutate(h2 = a2))
  }
  # across types of homogamy
  est_am <- rbind(est_am, est_h2 %>% mutate(assortment = am))
}

# PLOTTING --------------------------------------------------------------
# this is the starting h2 value we use as an illustrative example in text
a_example <- 0.30

## Figure 2a ####
# illustration of observed correlation under different causes of partner similarity
plot_a <- est_am %>%
  filter(component == "MZ" | component == "Full sib" | component == "MZ in-law" | component == "Full sib in-law") %>%
  filter(h2 == a_example) %>%
  mutate(
    component = factor(component, level = rev(c("MZ", "Full sib", "MZ in-law", "Full sib in-law"))),
    assortment = factor(assortment, level = c("dAM", "narrow-sense hom.", "broad-sense hom.", "social hom.", "idio. hom.", "mixed hom."))
  ) %>%
  ggplot(aes(assortment, est, fill = mu)) +
  geom_point(size = 3, shape = 21) +
  scale_fill_gradientn(colours = c(
    "#120246",
    "#050353",
    "#06175F",
    "#08306B",
    "#285C80",
    "#478495",
    "#68A7A9",
    "#89BDB5",
    "#AAD0C4",
    "#CCE3D7"
  )) +
  facet_grid(cols = vars(component), scales = "free") +
  labs(
    x = expression(paste("Causes of partner correlation")),
    y = expression(paste("Phenotypic correlation")),
    fill = expression(paste("Partner correlation"))
  ) +
  ylim(c(-.01, .5)) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

## Figure 2b ####
# differences in correlations with or without assortment
plot_b <- est_am %>%
  filter(component == "DZ") %>%
  # filter on the assortment that creates differences in DZ/sib correlations
  filter(assortment == "dAM" | assortment == "broad-sense hom." | assortment == "narrow-sense hom." | assortment == "mixed hom.") %>%
  mutate(mu = as.factor(mu)) %>%
  mutate(assortment = factor(assortment, levels = c("dAM", "narrow-sense hom.", "broad-sense hom.", "mixed hom."))) %>%
  pivot_wider(names_from = mu, values_from = est) %>%
  # we define delta as the difference between the sib correlation with or without assortment
  mutate(
    rm0.9 = `0.9` - `0`,
    rm0.8 = `0.8` - `0`,
    rm0.7 = `0.7` - `0`,
    rm0.6 = `0.6` - `0`,
    rm0.5 = `0.5` - `0`,
    rm0.4 = `0.4` - `0`,
    rm0.3 = `0.3` - `0`,
    rm0.2 = `0.2` - `0`,
    rm0.1 = `0.1` - `0`,
    rm0.0 = `0` - `0`
  ) %>%
  dplyr::select(h2, assortment, starts_with("rm")) %>%
  pivot_longer(rm0.9:rm0.0, names_to = "Partner_r", values_to = "dz_plus") %>%
  mutate(Partner_r = substr(Partner_r, 3, nchar(Partner_r))) %>%
  # plot delta curves
  ggplot(aes(h2, dz_plus, color = Partner_r)) +
  geom_vline(xintercept = a_example, linetype = "dashed") +
  facet_grid(cols = vars(assortment), scales = "free") +
  geom_line(size = 1) +
  theme_classic(base_size = 10) +
  scale_colour_manual(values = c(
    "#120246",
    "#050353",
    "#06175F",
    "#08306B",
    "#285C80",
    "#478495",
    "#68A7A9",
    "#89BDB5",
    "#AAD0C4",
    "#CCE3D7"
  )) +
  labs(
    x = expression(italic(h)^2),
    y = expression(paste(Delta, italic(r)[FS])),
    color = expression(paste(mu))
  ) +
  ylim(c(0, .35)) +
  scale_x_continuous(breaks = seq(a_example, 1, .15), limits = c(a_example, 0.85)) +
  theme(strip.background = element_blank())

# absolute bias in genetic effects
est_am_bias <- est_am %>%
  # select standardized components
  filter(grepl("var", component)) %>%
  mutate(
    d2 = .15, # we simulated dominance at .15
    component = case_match(
      component,
      "varP_E" ~ "E",
      "varP_T" ~ "T",
      "varP_D" ~ "D",
      "varP_C" ~ "C",
      "varP_E" ~ "E",
      "varP_D" ~ "D",
      "varP_C" ~ "C",
      "varP_N" ~ "N",
      "varP_N" ~ "N",
      "varP_A" ~ "A",
      "varP_A" ~ "A"
    )
  ) %>%
  filter(component %in% c("A", "D")) %>%
  pivot_wider(names_from = c("model", "component"), values_from = "est") %>%
  # compute bias if needed later
  mutate(
    biash2 = h2 - ADCE_A,
    biasd2 = d2 - ADCE_D
  )

## Figure 2c ####
# bias in h2twin estimates
plot_c <- est_am_bias %>%
  filter(assortment == "dAM" | assortment == "broad-sense hom." | assortment == "narrow-sense hom." | assortment == "mixed hom.") %>%
  mutate(assortment = factor(assortment, levels = c("dAM", "narrow-sense hom.", "broad-sense hom.", "mixed hom."))) %>%
  ggplot(aes(h2, ADCE_A, color = mu, group = mu)) +
  geom_point() +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = a_example, linetype = "dashed") +
  facet_grid(cols = vars(assortment), scales = "free") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = "dashed", color = "#B41F01") +
  theme_classic(base_size = 10) +
  ylim(c(.2, 1)) +
  scale_x_continuous(breaks = seq(a_example, 1, .15), limits = c(a_example, 0.85)) +
  scale_color_gradientn(colors = c(
    "#120246",
    "#050353",
    "#06175F",
    "#08306B",
    "#285C80",
    "#478495",
    "#68A7A9",
    "#89BDB5",
    "#AAD0C4",
    "#CCE3D7"
  )) +
  labs(
    x = expression(italic(h)^2),
    y = expression(hat(italic(h))[twin]^2),
    color = expression(paste("Partner correlation"))
  ) +
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )

## Figure 2d ####
# bias in d2twin estimates
plot_d <- est_am_bias %>%
  filter(assortment == "dAM" | assortment == "broad-sense hom." | assortment == "narrow-sense hom." | assortment == "mixed hom.") %>%
  mutate(assortment = factor(assortment, levels = c("dAM", "narrow-sense hom.", "broad-sense hom.", "mixed hom."))) %>%
  ggplot(aes(mu, ADCE_D, color = h2, group = h2)) +
  geom_hline(yintercept = .15, linetype = "dashed", color = "#FAA53C", linewidth = 1) +
  geom_point() +
  ylim(c(-0.01, 0.3)) +
  geom_line(linewidth = 1) +
  facet_grid(cols = vars(assortment), scales = "free") +
  theme_classic(base_size = 10) +
  scale_color_viridis_c(option = "plasma", breaks = seq(a_example, 1, .15)) +
  scale_x_continuous(limits = c(0, .95), breaks = seq(0, 1, .20)) +
  labs(
    y = expression(hat(italic(delta))[twin]^2),
    x = expression(mu),
    color = expression(italic(h)^2)
  ) +
  theme(strip.background = element_blank())

## Figure 2e ####
plot_e <- est_am_bias %>%
  filter(assortment == "dAM" | assortment == "broad-sense hom." | assortment == "narrow-sense hom." | assortment == "mixed hom.") %>%
  filter(mu != 0) %>%
  mutate(assortment = factor(assortment, levels = c("dAM", "narrow-sense hom.", "broad-sense hom.", "mixed hom."))) %>%
  ggplot(aes(h2, iamDATE_A, color = mu)) +
  geom_point() +
  geom_vline(xintercept = a_example, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = "dashed", color = "#B41F01") +
  scale_color_gradientn(colours = c(
    "#120246",
    "#050353",
    "#06175F",
    "#08306B",
    "#285C80",
    "#478495",
    "#68A7A9",
    "#89BDB5",
    "#AAD0C4",
    "#CCE3D7"
  )) +
  ylim(c(.2, 1)) +
  facet_grid(cols = vars(assortment), scales = "free") +
  theme_classic(base_size = 10) +
  scale_x_continuous(breaks = seq(a_example, 1, .15), limits = c(a_example, 0.85)) +
  labs(
    x = expression(italic(h)^2),
    y = expression(hat(italic(h))[twin - ped]^2),
    color = expression(paste("Partner correlation"))
  ) +
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )

## Figure 2f ####
plot_f <- est_am_bias %>%
  filter(assortment == "dAM" | assortment == "broad-sense hom." | assortment == "narrow-sense hom." | assortment == "mixed hom.") %>%
  mutate(assortment = factor(assortment, levels = c("dAM", "narrow-sense hom.", "broad-sense hom.", "mixed hom."))) %>%
  filter(mu != 0) %>%
  ggplot(aes(mu, iamDATE_D, color = h2)) +
  geom_point() +
  ylim(c(-0.01, 0.3)) +
  geom_hline(yintercept = .15, linetype = "dashed", color = "#FAA53C", linewidth = 1) +
  scale_color_viridis_c(option = "plasma", breaks = seq(a_example, 1, .15)) +
  scale_x_continuous(limits = c(0, .95), breaks = seq(0, 1, .2)) +
  facet_grid(cols = vars(assortment), scales = "free") +
  theme_classic(base_size = 10) +
  labs(
    x = expression(italic(mu)),
    y = expression(hat(italic(delta))[twin - ped]^2),
    color = expression(paste("Partner correlation"))
  ) +
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )

## Figure 2 ####
library(patchwork)
pdf(sprintf("%s/00_Figure1.pdf", wd_oa_figures),
  width = 10.75,
  height = 6
)
(plot_a + free(plot_b, side = "b") + plot_c + plot_d + plot_e + plot_f) +
  plot_layout(ncol = 2, widths = c(1, 1, 1), heights = c(1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a")
dev.off()