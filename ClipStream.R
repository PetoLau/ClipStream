source("ClusterStreamResults.R")
source("TestingForecasting.R")

testClipStream <- function(data, string, k.min, k.max, tresh, alpha, freq = seas, win = 21) {
  
  res_clipstream <- streamClust(data, k.min = k.min, k.max = k.max, outlier_tresh = tresh, alpha = alpha, freq = seas, win = win)
  
  gc()
  
  write.table(res_clipstream$data.clust, paste("Clip_Stream", string, win, "win", k.min, k.max, tresh, "outlier", alpha, "CD.csv", sep = "_"), row.names = F, col.names = T, quote = F)
  write.table(res_clipstream$info, paste("Clip_Stream_Info", string, win, "win", k.min, k.max, tresh, "outlier", alpha, "CD.csv", sep = "_"), row.names = F, col.names = T, quote = F)
  
  res_clus_clipstream <- ForecastClustersSimple(as.data.table(res_clipstream$data.clust))
  
  gc()
  
  write.table(res_clus_clipstream, paste("res_clipstream", string, win, "win", k.min, k.max, tresh, "outlier", alpha, "CD.csv", sep = "_"), row.names = F, col.names = F, quote = F)
}
