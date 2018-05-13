source("clusteringAll.R")

streamClust <- function(dataset, k.min, k.max, outlier_tresh, alpha, freq = 48, win = 21) {
  
  data.clust <- data.frame(N.slid.win = 0, N.K = 0, Load = 0)
  clust.info <- data.frame(N.slid.win = 0, N.K = 0,
                           N.out = 0, N.cd = 0, CD = F)
  
  oomDT <- dataset
  n_days <- ncol(oomDT)/freq - win
  
  # First iteration
  i <- 0
  data_oom <- oomDT[, ((i*freq)+1):((i+win)*freq), with = F]
  
  results <- clusterAll(data_oom, k.min, k.max, outlier_tresh = outlier_tresh)
  clip_oom <- results$representation
  
  for(l in 1:nrow(results$final_ts)){
    data.clust <- rbind(data.clust,
                        data.frame(N.slid.win = i, N.K = l, Load = results$final_ts[l,]))
  }
  
  clust.info <- rbind(clust.info, data.frame(N.slid.win = i, N.K = max(results$clustering),
                                             N.out = length(which(results$outliers == 0)), N.cd = 0, CD = F))
  
  cd_prev <- list(cd = FALSE, n_cd = 0)
  
  # Main cycle
  for(i in 1:(n_days-1)){
    
    cd_info <- cdTest(results$final_ts, alpha = alpha)
    
    if(cd_prev$cd == FALSE & cd_info$cd == FALSE & cd_info$n_cd > cd_prev$n_cd){
      cd_TF <- TRUE
    } else if(cd_info$cd == TRUE){
      cd_TF <- TRUE
    } else {
      cd_TF <- FALSE
    }
    
    cd_prev <- cd_info
    
    if(cd_info$cd == FALSE & cd_TF == TRUE){
      cd_prev$cd <- TRUE
    }
    
    clust.info[i+1, c("N.cd", "CD")] <- data.frame(N.cd = cd_info$n_cd, CD = cd_TF)
    
    if(cd_TF){
      
      data_oom <- oomDT[, ((i*freq)+1):((i+win)*freq), with = F]
      data_new <- data_oom[, tail(seq_len(ncol(data_oom)), freq), with = F]
      
      clip_new <- repr_matrix(data_new, func = repr_feaclip)
      clip_oom <- cbind(clip_oom[, -(1:ncol(clip_new))], clip_new)
      
      if(max(results$clustering) >= k.min + 2){
        k_min <- max(results$clustering) - 2
      } else k_min <- k.min
      
      if(max(results$clustering) <= k.max - 2){
        k_max <- max(results$clustering) + 2
      } else k_max <- k.max
      
      results <- clusterRepr(data_oom, clip_oom, k_min, k_max, outlier_tresh = outlier_tresh)
      
    } else {
      
      data_oom <- oomDT[, ((i*freq)+1):((i+win)*freq), with = F]
      data_new <- data_oom[, tail(seq_len(ncol(data_oom)), freq), with = F]
      
      clip_new <- repr_matrix(data_new, func = repr_feaclip)
      results$representation <- cbind(clip_oom[, -(1:ncol(clip_new))], clip_new)
      
      results$final_ts <- t(sapply(unique(results$clustering), function(x) colSums(check.matrix(data_oom[results$clustering == x,]))))
      # results$final_ts <- t(sapply(unique(results$merge_ts), 
      #                              function(x) colSums(check.matrix2(sums_clustered[results$merge_ts == x,]))))
      
    }
    
    for(l in 1:nrow(results$final_ts)){
      data.clust <- rbind(data.clust,
                          data.frame(N.slid.win = i, N.K = l, Load = results$final_ts[l,]))
    }
    
    clust.info <- rbind(clust.info, data.frame(N.slid.win = i, N.K = max(results$clustering),
                                               N.out = length(which(results$outliers == 0)), N.cd = 0, CD = F))
    
  }
  
  data.clust <- data.clust[-1,]
  clust.info <- clust.info[-1,]
  
  return(list(data.clust = data.clust, info = clust.info))
}
