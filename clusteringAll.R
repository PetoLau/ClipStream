source("DetectOutliers.R")
source("optimalClustering.R")
source("conceptDrift.R")
source("forecasting.R")

my_knn <- function(data, query) {
  as.integer(which.min(as.matrix(dist(rbind(data, query)))[nrow(data)+1,-(nrow(data)+1)]))
}

check.matrix <- function(mat) {
  if(is.matrix(mat) == TRUE){
    return(mat)
  } else
    return(as.matrix(mat))
}

check.matrix2 <- function(mat) {
  if(is.matrix(mat) == TRUE){
    return(mat)
  } else
    return(t(as.matrix(mat)))
}

clusterAll <- function(data, k_low, k_high, outlier_tresh, ncores = detectCores()-1, freq = 48){
  
  clip_oom <- repr_matrix(data, func = repr_feaclip,
                          windowing = T, win_size = freq)
  
  # detect outliers
  outliers <- detectOutliersIQR(clip_oom, treshold = outlier_tresh)
  
  # clustering
  clip_filtered <- clip_oom[-outliers$outliers,]
  clus_res <- clusterOptimKmedoidsDB(clip_filtered, k_low, k_high, ncores, criterium = "Davies_Bouldin")
  
  clustering <- outliers$class
  clustering[which(clustering == 1)] <- clus_res$clustering
  
  # Assign outliers to clusters (to nearest medoid)
  out_clust <- sapply(seq_along(outliers$outliers),
                      function(z) my_knn(clus_res$medoids, as.numeric(clip_oom[outliers$outliers[z],])))
  
  clustering[which(clustering == 0)] <- out_clust

  # Aggregating clusters
  final_ts_sums <- t(sapply(unique(clustering), function(x) colSums(check.matrix(data[clustering == x,]))))
  
  return(list(final_ts = final_ts_sums,
              outliers = outliers$class,
              clustering = clustering,
              representation = clip_oom))
}

clusterRepr <- function(data, data_rep, k_low, k_high, outlier_tresh, alpha, ncores = 4){
  
  # detect outliers
  outliers <- detectOutliersIQR(data_rep, treshold = outlier_tresh)
  
  # clustering
  clip_filtered <- data_rep[-outliers$outliers,]
  clus_res <- clusterOptimKmedoidsDB(clip_filtered, k_low, k_high, ncores, criterium = "Davies_Bouldin")
  
  clustering <- outliers$class
  clustering[which(clustering == 1)] <- clus_res$clustering
  
  # Assign outliers to clusters (to nearest medoid)
  out_clust <- sapply(seq_along(outliers$outliers),
                      function(z) my_knn(clus_res$medoids, as.numeric(data_rep[outliers$outliers[z],])))
  
  clustering[which(clustering == 0)] <- out_clust
  # table(clustering)
  
  # Aggregating clusters
  final_ts_sums <- t(sapply(unique(clustering), function(x) colSums(check.matrix(data[clustering == x,]))))
  
  return(list(final_ts = final_ts_sums,
              outliers = outliers$class,
              clustering = clustering))
}
