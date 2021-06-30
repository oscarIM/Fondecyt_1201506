#####funcion: entrega el hypervolumen para cada combinación pareada junto al promedio de todos los hypervolumnes por especie
#args: data: datos brutos, con una columna "species" (Quizas agregar argumento para ver si se quiere o no estandarizar y/o que método usar)
#cores: numero de cores que se quieran usar
#var_names: vector de los nombres de las variables a utilizar
#n_occs: numero de ocurrencias minimo por especie
#n_comb: numero de variables (i.e. dimensiones de nicho) en cada combinación (quizas agrgar un warning cuando hay mas dimensiones que occs).
#samples_per_points: numer de puntos aleatorios por ocurrencias empiricas (creo)
calc_hVol <- function(data, cores = NULL, var_names, n_occs, n_comb, samples_per_points, workers) {
  suppressPackageStartupMessages({
    require(hypervolume)
    require(dplyr)
    require(purrr)
    require(parallel)
    require(foreach)
    require(doParallel)
    require(tidyr)
    require(tictoc)
    require(furrr)
  })
  #data
  set.seed(12345)
  data_file <- data %>% filter(complete.cases(.))
  #vars to analyses
  data_hv <- data_file %>% select("species", all_of(var_names))
  data_hv <- as.data.frame(data_hv)
  ##filtrar >100
  summ_data <- data_hv %>% count(species) %>% filter(n >= n_occs)
  data_hv <- data_hv %>% filter(species %in% summ_data$species)
  ##estandarizar
  data_standarized <- data_hv %>% mutate_if(is.numeric, scale, center = TRUE, scale = TRUE)
  species_list <- split(data_standarized, f = data_standarized$species)
  #combinaciones de parametros
  cat("creating all variable combinations...\n")
  comb_pars <- as.data.frame(combn(var_names, n_comb))
  all_comb <- list()## pasar a foreach para consistencia de codigo
  for (i in 1:ncol(comb_pars)) {
    all_comb[[i]] <- map(species_list, ~select(., comb_pars[, i]))
  }
  #calcular los volumnes para cada combinacion en paralelo
  if (!is.null(cores) & cores >= 2) {
    cat("calculating hypervolume in parallel across combinations...\n")
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    plan(multicore, workers = workers)
    tic()
    vols <- foreach(i = 1:length(all_comb), .packages = c("future","furrr", "purrr","hypervolume")) %dopar% {
      furrr::future_map(all_comb[[i]], ~get_volume(hypervolume(data=., method = "box", samples.per.point = samples_per_points)), .progress = T)

    }
    exectime <- toc()
    exectime <- exectime$toc - exectime$tic
    #print(exectime)
    stopCluster(cl)
  }
  #calcular los volumnes de forma serial
  if (is.null(cores) || cores == 1) {
    cat("calculating hypervolume serially across combinations...\n")
    plan(multicore, workers = workers)
    tic()
    vols <- foreach(i = 1:length(all_comb), .packages = c("furrr","purrr","hypervolume")) %do% {
      furrr::future_map(all_comb[[i]], ~get_volume(hypervolume(data=., method = "box", kde.bandwidth = val_est_band, samples.per.point = samples_per_points)))
    }
    exectime <- toc()
    exectime <- exectime$toc - exectime$tic

  }
  cat("generating final dataframe...\n")
  vols_DF <- as.data.frame(t(bind_rows(vols)))
  vols_DF <- vols_DF %>% mutate(species = rownames(vols_DF)) %>% relocate(species)
  rownames(vols_DF) <- NULL
  df_log <- pivot_longer(vols_DF, cols = V1:last_col() , names_to = "combination", values_to = "hyperV")
  df_log_M <- df_log %>% group_by(species) %>% summarize(mean = mean(hyperV))
  vols_DF <- left_join(vols_DF, df_log_M, by = "species")
  vols_DF <- left_join(vols_DF, summ_data , by = "species")
  cat("ready...\n")
  return(vols_DF)

}

get_comb <- function(var_names, n_comb){
  comb_pars <- combn(var_names, n_comb)
  return(comb_pars)
}


