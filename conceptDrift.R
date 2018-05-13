ad_test_edit <- function(x) {
  
  if (length(unique(x)) == 1) {
    return(1)
  } else {
    return(ad.test(x)$ad[1, 3])
  }
}

cdTest <- function(data, alpha = 0.05, freq = 48) {
  
  p <- ncol(data)
  
  freq_win <- (p/freq)/7
  win <- p/freq_win
  
  # subset weeks and normalise them
  ts_test <- lapply(1:nrow(data),
                    function(i) lapply(0:(freq_win-1),
                                       function(j) norm_min_max(data[i,((win*j)+1):((win*j)+win)])))

  # Anderson-Darling test
  ad_pvalues <- sapply(seq_len(length(ts_test)),
                       function(x) ad_test_edit(ts_test[[x]]))
  
  # EDF change detected or not
  n_concept_drifts <- length(which(ad_pvalues < alpha))
  
  if(n_concept_drifts >= ceiling(nrow(data)/2)){ 
    res_cd <- TRUE
  } else res_cd <- FALSE
  
  return(list(n_cd = n_concept_drifts,
              cd = res_cd))
}
