# Davies-Bouldin index for k-Shape

DB_kshape <- function(clustering) {
  
  n <- max(clustering@cluster)
  
  R <- as.data.frame(matrix(NA,n,n))
  M <- as.data.frame(matrix(NA,n,n))
  D <- as.data.frame(matrix(NA,n,1))
  
  for(i in 1:n) {
    for(j in 1:n) {
      M[i,j] <- SBD(clustering@centroids[[i]], clustering@centroids[[j]], znorm = FALSE, error.check = TRUE)$dist
      if (i != j) {
        R[i,j] <- sum(clustering@clusinfo[i,2],clustering@clusinfo[j,2])/M[i,j]
      }
    }
    D[i,1] <- max(R[i,], na.rm = TRUE)
  }
  
  DB <- sum(D)/n
  
  return(DB)
}

# k-Shape clustering

clusterKShape <- function(matrixOOM, k_low, k_high, ncores){
  
  # PAA representation computing, z-score normalisation of original long time series windows
  mat <- repr_matrix(matrixOOM, func = repr_paa,
                     args = list(q = 6, func = mean), normalise = TRUE,
                     func_norm = norm_z)
  
  mat_sd <- apply(mat, 1, sd)
  
  for(i in which(mat_sd == 0)) {
    mat[i,] <- rnorm(ncol(mat), 0, 0.01)
  }
  
  cl <- makeCluster(ncores)
  clusterExport(cl, varlist = c("mat", "k_low", "k_high"), envir = environment())
  clusterCall(cl, function() library(dtwclust))
  clusterCall(cl, function() library(proxy))
  clusterings <- parLapply(cl, c(k_low:k_high), function(x) tsclust(mat, k = x,
                                                                    distance = "sbd", centroid = "shape", type = "partitional", preproc = zscore,
                                                                    control = partitional_control(pam.precompute = FALSE), trace = FALSE))
  if(!is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- c()
  }
  
  cl <- makeCluster(ncores)
  clusterExport(cl, varlist = c("clusterings", "DB_kshape", "SBD"), envir = environment())
  DB_values <- parSapply(cl, seq_along(clusterings), function(x) DB_kshape(clusterings[[x]]))
  if(!is.null(cl)){
    parallel::stopCluster(cl)
    cl <- c()
  }
  
  return(clusterings[[which.min(DB_values)]]@cluster)

}
