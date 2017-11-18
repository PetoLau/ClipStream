## Representations of Time Series ----
# DFT - FFT
fftinv <- function(x){
  fft(x, inverse = TRUE) / length(x)
  # fft(x, inverse = TRUE)
}

makeRepFourier <- function(y, koef){
  fourier_fft <- fft(y)
  inv_fft <- fftinv(fourier_fft[1:koef])

  return(as.vector(Re(inv_fft)))
}

DFTRep <- function(matTS){
  mat <- data.matrix(matTS)

  cl <- makeCluster(detectCores())
  clusterExport(cl, varlist = c("mat", "fftinv", "makeRepFourier"), envir = environment())
  DFT <- parSapply(cl, 1:dim(mat)[1], function(i) makeRepFourier(mat[i,], 6))
  if(!is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- c()
  }

  return(as.data.table(t(DFT)))
}
