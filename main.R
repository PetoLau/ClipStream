# Header ----
# library(ggplot2)
library(data.table)
library(parallel)
library(cluster)
library(clusterCrit)
library(TSrepr)
library(kSamples)
library(rpart)
library(party)
library(smooth)
library(forecast)
library(doParallel)
library(randomForest)
library(dtwclust)
library(foreign)

# Data reading ----
# data must be in the format where time series streams are in rows of a matrix or a data.frame
data <- as.data.table(read.arff("London_5months_5066ID.ARFF"))
seas <- 48

# Offline batch clustering - benchmarks ----
source("batchClust.R")

batchClust(data, string = "London", k.min = 20, k.max = 30, method = "lm")
batchClust(data, string = "London", k.min = 20, k.max = 30, method = "dft")
batchClust(data, string = "London", k.min = 20, k.max = 30, method = "kshape")

# ClipStream ----
source("ClipStream.R")

testClipStream(data, string = "London", k.min = 20, k.max = 30, tresh = 1.5, alpha = 0.05)

# Simple aggregate forecasting ----
source("TestingForecasting.R")

data_sum <- colSums(data)
res_sim <- ForecastAggregatedSimple(data_sum)
err_sim <- computeMape(data_sum, res_sim)
gc()

write.table(res_sim_london, "res_sim.csv", row.names = F, col.names = F, quote = F)

# Evaluation of results ----
err_agg <- err_sim$ByDay
res_clipstream <- fread("result_of_clipstream.csv")
err_clip <- computeMape(data_sum, res_clipstream)

err_clip <- err_clip$ByDay
err_agg <- err_agg$ByDay

colMeans(err_agg)
colMeans(err_clip)

wilcErr <- function(x, y) {
  wilcox.test(x, y, paired = T, alternative = "less")$p.value
}

round(as.matrix(sapply(seq_len(dim(err_agg)[2]),
                       function(i) wilcErr(err_clip[,i], err_agg[,i]))), digits = 6)
