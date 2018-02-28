# Header ----
# rm(list=ls())
# gc()

# library(ggplot2)
library(data.table)
library(parallel)
library(cluster)
library(clusterCrit)
library(TSrepr) # install via devtools::install_github('petolau/TSrepr')
library(kSamples)
library(rpart)
library(party)
library(smooth)
library(forecast)
library(doParallel)
library(randomForest)

# Data reading ----
# data must be in the format where time series streams are in rows of a matrix or a data.frame
data <- fread("path")
seas <- 48

# Offline batch clustering - DFT ----
source("optimalClustering.R") # K-medoids

data.clust <- data.frame(N.slid.win = 0, N.K = 0, Load = 0)
# clust.res <- data.frame(N.slid.win = 0, Class = 0, ID = 0)

data_oom <- data
win <- 21
days_for <- ncol(data_oom)/seas - win

for(k in 0:(days_for-1)){
  
  oomDT.select <- data_oom[, ((k*seas)+1):((k+win)*seas), with = F]

  reprs <- repr_matrix(oomDT.select, func = repr_dft, args = list(coef = 6),
                       normalise = TRUE, func_norm = norm_min_max)
  
  km_res <- clusterOptimKmedoidsDB(reprs, 8, 18, 8)$clustering
  
  km_sums <- t(sapply(1:length(unique(km_res)), function(x) colSums(oomDT.select[km_res == x,])))
  
  for(l in 1:length(unique(km_res))){
    data.clust <- rbind(data.clust,
                        data.frame(N.slid.win = k, N.K = l, Load = km_sums[l,]))
  }
  
  # clust.res <- rbind(clust.res,
  #                    data.frame(N.slid.win = k, Class = km_res, ID = 1:nrow(oomDT.select)))
  
  print(k)
}

data.clust <- data.clust[-1,]
summary(data.clust)
str(data.clust)

write.table(data.clust, "DFT_Kmed.csv", row.names = F, col.names = T, quote = F)

# Multiple Data Streams Clustering - ClipStream with FeaClip ----
source("ClusterStreamResults.R")

# optional - subset time series data
# data <- trainSetHack[, 1:(seas*28), with = F]

res_clipstream <- streamClust(data, k.min = 8, k.max = 16, freq = seas)

summary(res_clipstream$data.clust)
summary(res_clipstream$info)

plot(ts(res_clipstream$info$N.K))
plot(ts(res_clipstream$info$N.K.final))
plot(ts(res_clipstream$info$N.out))
plot(ts(res_clipstream$info$N.cd))

write.table(res_clipstream$data.clust, "Clip_Stream_Kmed.csv", row.names = F, col.names = T, quote = F)
write.table(res_clipstream$info, "Clip_Info_Stream_Kmed.csv", row.names = F, col.names = T, quote = F)

# Forecasting on clusters and simple agg. ----
source("TestingForecasting.R")
files <- list.files(pattern = "*.csv")

data.clust <- fread(files[1])
data_sum <- colSums(data)
res_clus_clipstream <- ForecastClustersSimple(data.clust)
err_clus_clipstream <- computeMape(data_sum, res_clus_clipstream)
gc()

data.clust <- fread(files[2])
data_sum <- colSums(data)
res_clus_dft <- ForecastClustersSimple(data.clust)
err_clus_dft <- computeMape(data_sum, res_clus_dft)
gc()

write.table(res_clus_clipstream, "res_clipstream.csv", row.names = F, col.names = F, quote = F)
write.table(res_clus_dft, "res_dft.csv", row.names = F, col.names = F, quote = F)

# simple aggregate forecasting ----

data_sum <- colSums(data)
res_sim <- ForecastAggregatedSimple(data_sum)
err_sim <- computeMape(data_sum, res_sim)
gc()

write.table(res_sim_1, "res_sim.csv", row.names = F, col.names = F, quote = F)

# Evaluation of results ----
err_agg <- err_sim$ByDay
err_clip <- err_clus_clipstream$ByDay
err_dft <- err_clus_dft$ByDay

colMeans(err_agg)
colMeans(err_clip)
colMeans(err_dft)

wilcErr <- function(x, y){
  wilcox.test(x, y, paired = T, alternative = "less")$p.value
}

round(as.matrix(sapply(seq_len(dim(err_agg)[2]),
                       function(i) wilcErr(err_clip[,i], err_agg[,i]))), digits = 6)

round(as.matrix(sapply(seq_len(dim(err_agg)[2]),
                       function(i) wilcErr(err_dft[,i], err_agg[,i]))), digits = 6)
