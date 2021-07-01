#####función: entrega el hipervolumen para cada combinación pareada junto al promedio de todos los hipervolumenes por especie
#####args####
#data: datos brutos,
#cores: numero de cores que se quieran usar
#var_names: vector de los nombres de las variables a utilizar
#n_occs: numero de ocurrencias mínimo por especie
#n_comb: numero de variables (i.e. dimensiones de nicho) en cada combinación
#samples_per_points: numero de puntos aleatorios por ocurrencias empíricas (creo)
calc_hVol <- function(data, cores, var_names, n_occs, n_comb, samples_per_points) {
  suppressPackageStartupMessages({
    require(hypervolume)
    require(dplyr)
    require(purrr)
    require(parallel)
    require(foreach)
    require(doParallel)
    require(tidyr)
    require(tictoc)
    })
  #datos
  options(scipen = 999)
  set.seed(12345)
  data_file <- data %>% filter(complete.cases(.))
  #seleccion variables
  data_hv <- data_file %>% select("species", all_of(var_names))
  data_hv <- as.data.frame(data_hv)
  ##filtrado por numero de ocurrencias
  summ_data <- data_hv %>% count(species) %>% filter(n >= n_occs)
  data_hv <- data_hv %>% filter(species %in% summ_data$species)
  ##estandarización
  data_standarized <- data_hv %>% mutate_if(is.numeric, scale, center = TRUE, scale = TRUE)
  species_list <- split(data_standarized, f = data_standarized$species)
  #combinaciones de parámetros
  cat("creando todas las combinaciones de variables...\n")
  comb_pars <- as.data.frame(combn(var_names, n_comb))
  all_comb <- list()## pasar a foreach para consistencia de código
  for (i in 1:ncol(comb_pars)) {
    all_comb[[i]] <- map(species_list, ~select(., comb_pars[, i]))
    }
  #calcular el volumen para cada combinación en paralelo
  if (cores >= 2) {
    cat("calculando hipervolumenes en paralelo...\n")
    cl <- parallel::makeCluster(cores, setup_timeout = 0.5)
    registerDoParallel(cl)
    tic("tiempo de ejecución")
    vols <- foreach(c = 1:length(all_comb)) %:%
      foreach(n = 1:length(all_comb[[1]]), .packages = "hypervolume") %dopar% {
        get_volume(hypervolume(data = all_comb[[c]][[n]], method = "box", samples.per.point = samples_per_points))
        }
    toc(log = TRUE)
    stopCluster(cl)
    }
  if (cores == 1) {
    cat("calculando hipervolumnes en serie...\n")
    tic("tiempo de ejecución")
    vols <- foreach(c = 1:length(all_comb), .packages = "hypervolume") %:%
      foreach(n = 1:length(all_comb[[1]])) %do% {
        get_volume(hypervolume(data = all_comb[[c]][[n]], method = "box", samples.per.point = samples_per_points))
        }
    toc(log = TRUE)
    }
  cat("generando resultados...\n")
  suppressMessages(tmp <- map(vols, ~do.call(rbind, .)) %>% bind_cols() %>% mutate(species = names(species_list)) %>% relocate(species))
 names(tmp)[2:length(tmp)] <- paste0("comb_", seq(1:length(all_comb)))
 tmp_log <- pivot_longer(tmp, cols = comb_1:last_col() , names_to = "combinations", values_to = "hyperV")
 tmp_log_M <- tmp_log %>% group_by(species) %>% summarize(mean = mean(hyperV))
 vols_DF <- left_join(tmp, tmp_log_M, by = "species")
 vols_DF <- left_join(vols_DF, summ_data , by = "species")
 cat("!listo!...\n")
 return(vols_DF)
}
get_comb <- function(var_names, n_comb){
  comb_pars <- combn(var_names, n_comb)
  return(comb_pars)
}


