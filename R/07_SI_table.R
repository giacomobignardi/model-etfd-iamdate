#Author: Giacomo Bignardi
# adapted from https://github.com/giacomobignardi/h2_BMRQ/blob/main/R/04_r1_MRS_facets_correlations.R 
# prepare supplementary file
library(openxlsx)
library(tidyverse)
rm(list = ls())

# set open access working directories
wd_oa = getwd()

# load sup tables
st1 <- read.csv(sprintf("%s/SI/02_table1_dem.csv", wd_oa))
st2 <- read.csv(sprintf("%s/SI/03_table2_esat_mxcompare.csv", wd_oa))
st3 <- read.csv(sprintf("%s/SI/06_table3_mod_comp.csv", wd_oa))
st4 <- read.csv(sprintf("%s/SI/06_table4_est.csv", wd_oa))
st5 <- read.csv(sprintf("%s/SI/06_table5_est.csv", wd_oa))

# order of SI st1, st5, st2, st3, st4
# create a reference table
st0 <- data.frame(sheet = c(paste0("S", seq(0,5,1))),
                  content = c(
                    "Overview of sheet contents: provides information about the contents of all other sheets",
                    "Complete data per family member: displays the number of individuals with complete data for each family member",
                    "Alternative estimates dAM-ANE across rNA values: presents alternative parameter estimates from the dAM-ANE model across varying rNA values",
                    "Model comparison saturated model: contains results from the comparison against the saturated model",
                    "Model comparison iAM-DATE model: contains results from the comparison against the iAM-DATE",
                    "Parameter estimates iAM-DATE and dAM-ADE: shows parameter estimates derived from the iAM-DATE and the most parsimonious dAM-ADE"
                    ))
st_list = list(st0, st1, st5, st2, st3, st4)

# create a workbook
wb <- createWorkbook()

# create a list for sheet names
sn <- c("0.overview", 
        "1.n", 
        "2.dAMANE_alt_rNA", 
        "3.sat_mod_comp", 
        "4.iAMDATE_mod_comp", 
        "5.iAMDATE_parameter_est")

# add each correlation matrix as a sheet
for (i in 1:length(st_list)) {
  sheet_name <- sn[i]
  addWorksheet(wb, sheet_name)  # Add a new sheet
  writeData(wb, sheet = sheet_name, x = st_list[[i]]) 
}

# add note to sup 5
note_sup5 <- "Notes. The variance (var) is standardised; the parameter labels for the path coefficients to the sorting factor end in s (e.g., as);
varP = variance of the focal phenotype (i.e. self reported proneness to chills); varS = variance of the sorting factor"
writeData(wb, "5.iAMDATE_parameter_est", note_sup5, startCol = 1, startRow = 60)

# save Supplementary file
saveWorkbook(wb,sprintf("%s/SI/Supplementary_Table.xlsx", wd_oa), overwrite = T)
