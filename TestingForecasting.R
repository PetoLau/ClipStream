source("forecasting.R")

pickwindow <- function(data, N.win){
  
  data_sub <- data[N.slid.win %in% N.win, lapply(Load, as.vector), by = N.K]
  data_sub[, N.K := NULL]
  data_sub <- data.matrix(data_sub)
  
  data_sub
}

ForecastClustersSimple <- function(dataset){
  
  n_day <- length(unique(dataset$N.slid.win))
  
  data_list <- lapply(0:(n_day-1), function(x) pickwindow(dataset, x))
  
  cl <- makeForkCluster(detectCores()-1, outfile = "")
  registerDoParallel(cl)
  pred_clusters <- parLapply(cl, 1:(length(data_list)),
                          function(i) lapply(1:dim(data_list[[i]])[1], 
                                             function(j) predSimAll(data_list[[i]][j,])))
  if(!is.null(cl)){
    parallel::stopCluster(cl)
    cl <- c()
  }

  res <- sapply(seq_len(length(pred_clusters[[1]][[1]])),
                function(k) as.vector(sapply(seq_len(length(pred_clusters)),
                                             function(i) rowSums(sapply(seq_len(length(pred_clusters[[i]])),
                                                                        function(j) as.vector(pred_clusters[[i]][[j]][[k]]))))))
  
  return(res)
}

ForecastAggregatedSimple <- function(dataset){
  
  win <- 21
  days_for <- length(dataset)/seas - win
  
  pred_sums <- lapply(0:(days_for-1),
                      function(i) predSimAll(dataset[((i*seas)+1):((seas*i)+(win*seas))]))
  
  predictions <- sapply(seq_len(length(pred_sums[[1]])),
                        function(k) sapply(seq_len(days_for),
                                           function(j) as.vector(pred_sums[[j]][[k]])))
  
  return(predictions)
}

computeMape <- function(real, predictions){
  
  win <- 21
  data_test <- real[-(1:(win*seas))]
  n_day <- length(data_test)/seas
  
  err_byday <- sapply(seq_len(dim(predictions)[2]),
                      function(i) sapply(0:(n_day-1),
                                         function(j) mape(as.vector(data_test[((j*seas)+1):(seas*(j+1))]), predictions[((j*seas)+1):(seas*(j+1)), i])))
  
  err_whole <- sapply(seq_len(dim(predictions)[2]),
                      function(i) mape(as.vector(data_test), predictions[, i]))
  
  return(list(ByDay = err_byday, Whole = matrix(err_whole)))
}
