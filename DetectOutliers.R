detectOutliersCustom <- function(data) {
  
  # data <- results$representation
  
  p <- ncol(data)
  
  clip_out_10 <- data[, seq(2, p, by = 8), with = F]
  # clip_out3 <- data[, seq(4, p, by = 9), with = F]
  clip_out_jumps <- data[, seq(4, p, by = 8), with = F]
  
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
