# Outlier Detection ----
detectOutliersIQR <- function(data, treshold = 1.5) {
  
  p <- ncol(data)
  
  # subset sum_1 and crossings
  clip_out_10 <- data[, seq(2, p, by = 8)]
  clip_out_jumps <- data[, seq(4, p, by = 8)]
  
  # compute quartiles and find outliers - Q_0.75 + treshold*IQR or Q_0.25 - treshold*IQR
  # usually treshold = 1.5
  jumps <- rowMeans(clip_out_jumps)
  jumps <- which(jumps > quantile(jumps, probs = 0.75) + treshold * (quantile(jumps, probs = 0.75) - quantile(jumps, probs = 0.25)))
  zero_one <- rowMeans(clip_out_10)
  under <- which(zero_one > quantile(zero_one, probs = 0.75) + treshold * (quantile(zero_one, probs = 0.75) - quantile(zero_one, probs = 0.25)))
  upper <- which(zero_one < quantile(zero_one, probs = 0.25) - treshold * (quantile(zero_one, probs = 0.75) - quantile(zero_one, probs = 0.25)))
  
  class <- rep(1, nrow(data))
  class[c(jumps, under, upper)] <- 0
  
  return(list(outliers = as.vector(which(class == 0)),
              class = class))
}
