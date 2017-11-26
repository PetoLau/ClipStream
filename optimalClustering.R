clusterOptimKmedoidsDB <- function(matrixOOM, k_low, k_high, ncores, criterium = "Davies_Bouldin"){

  mat <- data.matrix(matrixOOM)
  
  cl <- makeCluster(ncores)
  clusterExport(cl, varlist = c("mat", "pam", "k_low", "k_high"), envir = environment())
  clusterings <- parLapply(cl, c(k_low:k_high), function(x) pam(mat, x, metric = "euclidean"))
  if(!is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- c()
  }

  cl <- makeCluster(ncores)
  clusterExport(cl, varlist = c("clusterings", "mat", "intCriteria", "criterium"), envir = environment())
  DB_values <- parSapply(cl, seq_along(clusterings), function(x) intCriteria(mat, as.integer(clusterings[[x]]$clustering), criterium))
  if(!is.null(cl)){
    parallel::stopCluster(cl)
    cl <- c()
  }

  if(criterium == "Davies_Bouldin") {
    return(clusterings[[which.min(DB_values)]])
  } else return(clusterings[[which.max(DB_values)]])
  
}
