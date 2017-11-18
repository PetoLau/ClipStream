source("RepresentationsWeekly.R")
source("DetectOutliers.R")
source("optimalClustering.R")
source("mergeSums.R")
source("conceptDrift.R")
source("forecasting.R")

my_knn <- function(data, query) {
  as.integer(which.min(as.matrix(dist(rbind(data, query)))[nrow(data)+1,-(nrow(data)+1)]))
}

clusterAll <- function(data, k_low, k_high, ncores = detectCores()-1, freq = 48){
  
  # clip_oom <- FeaClipRep(data)
  clip_oom <- repr_matrix(data,
                          func = repr_feaclip, windowing = T, win_size = freq)
  
  # detect outliers
  outliers <- detectOutliersCustom(clip_oom)
  
  # clustering
  clip_filtered <- clip_oom[-outliers$outliers]
  clus_res <- clusterOptimKmedoidsDB(clip_filtered, k_low, k_high, ncores, criterium = "Davies_Bouldin")
  
  clustering <- outliers$class
  clustering[which(clustering == 1)] <- clus_res$clustering
  
  # Assign outliers to clusters (to nearest medoid)
  out_clust <- sapply(seq_along(outliers$outliers),
                      function(z) my_knn(clus_res$medoids, as.numeric(clip_oom[outliers$outliers[z],])))
  
  clustering[which(clustering == 0)] <- out_clust
  # table(clustering)
  
  # Merge clusters
  final_ts_sums <- mergeClusters(data, clustering)
  
  return(list(final_ts = final_ts_sums$final_ts,
              merge_ts = final_ts_sums$clustering,
              outliers = outliers$class,
              clustering = clustering,
              representation = clip_oom))
}

clusterRepr <- function(data, data_rep, k_low, k_high, ncores = 3){
  
  # detect outliers
  outliers <- detectOutliersCustom(data_rep)
  
  # clustering
  clip_filtered <- data_rep[-outliers$outliers]
  clus_res <- clusterOptimKmedoidsDB(clip_filtered, k_low, k_high, ncores, criterium = "Davies_Bouldin")
  
  clustering <- outliers$class
  clustering[which(clustering == 1)] <- clus_res$clustering
  
  # Assign outliers to clusters (to nearest medoid)
  out_clust <- sapply(seq_along(outliers$outliers),
                      function(z) my_knn(clus_res$medoids, as.numeric(data_rep[outliers$outliers[z],])))
  
  clustering[which(clustering == 0)] <- out_clust
  # table(clustering)
  
  # Merge clusters
  final_ts_sums <- mergeClusters(data, clustering)
  
  return(list(final_ts = final_ts_sums$final_ts,
              merge_ts = final_ts_sums$clustering,
              outliers = outliers$class,
              clustering = clustering))
}
