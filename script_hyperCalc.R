library(data.table)
library(dplyr)
source("hyperCalc.R")
##
data <- fread("All_Mammals_productivity.csv")
#variables que pasaron el VIF. Solo como ejemplo
var_names <- c("PETDriestQuarter", "PETseasonality", "PETWettestQuarter")
#ver el numero de combinaciones pareadas
get_comb(var_names = var_names, n_comb = 2)
##cada core use como 600 Mb de RAM
hyperV_prod <- calc_hVol(data = data, cores = 20, var_names = var_names, n_occs = 100, n_comb = 2, samples_per_points = 10)
