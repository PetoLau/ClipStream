cdTest <- function(data, alpha = 0.05, freq = 48){
  
  # Normalise time series
  ts_news <- sapply(seq_len(nrow(data)),
                    function(x) norm_min_max(tail(data[x,], 3*freq)))
  
  ts_olds <- sapply(seq_len(nrow(data)),
                    function(x) norm_min_max(head(tail(data[x,], 10*freq), 3*freq)))
  
  # Anderson-Darling test
  ad_pvalues <- sapply(seq_len(nrow(data)),
                       function(x) ad.test(ts_olds[,x], ts_news[,x])$ad[1, 3])
  
  # EDF change detected or not
  n_concept_drifts <- length(which(ad_pvalues < alpha))
  
  if(n_concept_drifts >= ceiling(nrow(data)/2)){ 
    res_cd <- TRUE
  } else res_cd <- FALSE
  
  return(list(n_cd = n_concept_drifts,
              cd = res_cd))
}
