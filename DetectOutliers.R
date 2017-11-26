detectOutliersCustom <- function(data) {
  
  p <- ncol(data)
  
  # subset sum_1 and jumps
  clip_out_10 <- data[, seq(2, p, by = 8)]
  clip_out_jumps <- data[, seq(4, p, by = 8)]
  
  # compute quartiles and find outliers
  jumps <- rowMeans(clip_out_jumps)
  jumps <- which(jumps > 1.5*summary(jumps)[5])
  zero_one <- rowMeans(clip_out_10)
  under <- which(zero_one > 1.25*summary(zero_one)[5])
  upper <- which(zero_one < 0.75*summary(zero_one)[2])
  
  class <- rep(1, nrow(data))
  class[c(jumps, under, upper)] <- 0
  
  return(list(outliers = as.vector(which(class == 0)),
              class = class))
}
