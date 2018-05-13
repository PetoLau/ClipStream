source("optimalClustering.R") # K-medoids
source("kshape.R")
source("TestingForecasting.R")

batchClust <- function(data, string, k.min, k.max, method) {
  
  data.clust <- data.frame(N.slid.win = 0, N.K = 0, Load = 0)
  
  data_oom <- data
  win <- 21
  days_for <- ncol(data_oom)/seas - win
  
  for(k in 0:(days_for-1)){
    
    oomDT.select <- data_oom[, ((k*seas)+1):((k+win)*seas), with = F]
    
    if (method == "lm") {
      
      reprs <- repr_matrix(oomDT.select, func = repr_lm,
                           args = list(freq = c(seas, seas*7), method = "lm"),
                           normalise = TRUE, func_norm = norm_min_max)
      
      clustering <- clusterOptimKmedoidsDB(reprs, k.min, k.max, ncores = 7, criterium = "Davies_Bouldin")$clustering
      
    }
    
    if (method == "dft") {
      
      reprs <- repr_matrix(oomDT.select, func = repr_dft, args = list(coef = 8*21),
                           normalise = TRUE, func_norm = norm_min_max)
      
      clustering <- clusterOptimKmedoidsDB(reprs, k.min, k.max, ncores = 7, criterium = "Davies_Bouldin")$clustering
      
    }
    
    if (method == "feaclip") {
      
      reprs <- repr_matrix(oomDT.select, func = repr_feaclip,
                           windowing = TRUE, win_size = seas)
      
      clustering <- clusterOptimKmedoidsDB(reprs, k.min, k.max, ncores = 7, criterium = "Davies_Bouldin")$clustering
      
    }
    
    if (method == "kshape") {
      
      clustering <- clusterKShape(oomDT.select, k_low = k.min, k_high = k.max, ncores = 7)
      
    }
    
    km_sums <- t(sapply(1:length(unique(clustering)), function(x) colSums(oomDT.select[clustering == x,])))
    
    for(l in 1:nrow(km_sums)){
      data.clust <- rbind(data.clust,
                          data.frame(N.slid.win = k, N.K = l, Load = km_sums[l,]))
    }
    
    # print(k)
  }
  
  data.clust <- data.clust[-1,]
  
  gc()
  
  write.table(data.clust, paste("Kmed", method, string, k.min, k.max, ".csv", sep = "_"), row.names = F, col.names = T, quote = F)
  
  res_clus_for <- ForecastClustersSimple(as.data.table(data.clust))
  
  gc()
  
  write.table(res_clus_for, paste("res", method, string, k.min, k.max, ".csv", sep = "_"), row.names = F, col.names = F, quote = F)
  
}
