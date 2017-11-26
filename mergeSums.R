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

mergeClusters <- function(data, clustering, freq = 48) {
  
  period <- freq * 7
  slid.win <- ncol(data) / freq
  
  sums_clustered <- t(sapply(unique(clustering), function(x) colSums(check.matrix(data[clustering == x,]))))

  clip_sums <- repr_matrix(sums_clustered,
                           func = repr_feaclip, windowing = T, win_size = freq)
  
  clip_sums <- repr_matrix(clip_sums,
                           func = repr_seas_profile, args = list(freq = 8, func = meanC))
  
  clus_res_sums <- clusterOptimKmedoidsDB(clip_sums, 3, max(clustering)-1, ncores = 2, criterium = "Davies_Bouldin")

  final_ts_sums <- sapply(unique(clus_res_sums$clustering), 
                          function(x) colSums(check.matrix2(sums_clustered[clus_res_sums$clustering == x,])))
  
  return(list(final_ts = t(final_ts_sums),
              clustering = clus_res_sums$clustering))
}
