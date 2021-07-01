library(data.table)
library(dplyr)
source("hyperCalc.R")
##
data <- fread("All_mammals_bio_produ_VIF.csv")
marine_sp <- read.csv("marines.csv")
#eliminar sp marinas
data_terr <- data %>% filter(!species %in% marine_sp$marines)
var_names <- var_names[8:length(var_names)]
#ver el numero de combinaciones pareadas
get_comb(var_names = var_names, n_comb = 2)
##cada core use como 600 Mb de RAM
hyperV_prod <- calc_hVol(data = data_terr, cores = 20, var_names = var_names, n_occs = 100, n_comb = 2, samples_per_points = 10)
write.csv(hyperV_prod, "hyperV_prod_biclim.csv", quote = F, row.names = F)
