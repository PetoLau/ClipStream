simRpart <- function(Y, K, freq = 48, h = 48){
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K*2))
  fuur_test <- as.data.frame(fourier(data_ts, K = c(K, K*2), h = h))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             fuur)
  
  tree_1 <- rpart(Load ~ ., data = matrix_train,
                  control = rpart.control(minsplit = 2,
                                          maxdepth = 30,
                                          cp = 0.000001))
  
  # new data and prediction
  pred_tree <- predict(tree_1, fuur_test) + trend_for
  
  return(as.vector(pred_tree))
}

simCtreeLag <- function(Y, freq = 48, h = 48){
  
  N <- length(Y)
  window <- (N / freq) - 1
  
  K <- 2
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K*2))
  fuur_test <- as.data.frame(fourier(data_ts, K = c(K, K*2), h = h))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  train_s <- decom_1$time.series[1:(freq*window), 1]
  
  matrix_trigonom <- data.table(Load = tail(new_load, window*freq),
                                fuur[(freq+1):N,],
                                Seas.Ts = train_s)
  
  test_s <- decom_1$time.series[((freq*window)+1):(freq*(window+1)), 1]
  
  matrix_test <- data.table(fuur_test,
                            Seas.Ts = test_s)
  
  tree_2 <- party::ctree(Load ~ ., data = matrix_trigonom,
                         controls = party::ctree_control(teststat = c("quad"),
                                                         testtype = c("Teststatistic"),
                                                         mincriterion = 0.95,
                                                         minsplit = 1,
                                                         minbucket = 1,
                                                         mtry = 0, maxdepth = 0))
  
  pred_tree <- predict(tree_2, matrix_test) + trend_for
  
  return(as.vector(pred_tree))
}

simRFLag <- function(Y, freq = 48, h = 48){
  
  N <- length(Y)
  window <- (N / freq) - 1
  
  K <- 2
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K*2))
  fuur_test <- fourier(data_ts, K = c(K, K*2), h = h)
  
  colnames(fuur) <- NULL
  colnames(fuur_test) <- NULL
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  train_s <- decom_1$time.series[1:(freq*window), 1]
  
  matrix_trigonom <- data.table(Load = tail(new_load, window*freq),
                                fuur[(freq+1):N,],
                                Seas.Ts = train_s)
  
  test_s <- decom_1$time.series[((freq*window)+1):(freq*(window+1)), 1]
  
  matrix_test <- data.table(fuur_test,
                            Seas.Ts = test_s)
  
  rf_model <- randomForest(Load ~ ., data = matrix_trigonom,
                           ntree = 1100, mtry = 3, nodesize = 3, importance = TRUE)
  
  pred_rf <- predict(rf_model, matrix_test) + trend_for
  
  return(as.vector(pred_rf))
}

simCtreeFur <- function(Y, K, freq = 48, h = 48){
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K*2))
  fuur_test <- as.data.table(fourier(data_ts, K = c(K, K*2), h = h))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             fuur)
  
  tree_2 <- party::ctree(Load ~ ., data = matrix_train,
                         controls = party::ctree_control(teststat = c("quad"),
                                                         testtype = c("Teststatistic"),
                                                         mincriterion = 0.95,
                                                         minsplit = 1,
                                                         minbucket = 1,
                                                         mtry = 0, maxdepth = 0))
  
  pred_tree <- predict(tree_2, fuur_test) + trend_for
  
  return(as.vector(pred_tree))
}

simRFFur <- function(Y, K, freq = 48, h = 48){
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K*2))
  fuur_test <- fourier(data_ts, K = c(K, K*2), h = h)
  
  colnames(fuur) <- NULL
  colnames(fuur_test) <- NULL
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             fuur)
  
  rf_model <- randomForest(Load ~ ., data = matrix_train,
                           ntree = 1100, mtry = 4,
                           nodesize = 3, importance = TRUE)
  
  pred_tree <- predict(rf_model, data.table(fuur_test)) + trend_for
  
  return(as.vector(pred_tree))
}

simEsExtra <- function(Y, freq = 48, h = 48){
  
  data_ts <- ts(Y, freq = freq*7)
  
  pred <- es(data_ts, model = c("AAA", "ANA", "AAdA"), ic = c("AIC"), h = h, cfType = "MAE", intervals = "none", silent = "all")$forecast
  
  return(as.vector(pred))
}

# All methods in one FUN ----

predSimAll <- function(Y, freq = 48, h = 48){
  
  pred_rpart_2 <- simRpart(Y, 2, freq = freq, h = h)
  pred_ctree_lag <- simCtreeLag(Y, freq = freq, h = h)
  pred_ctree_fur <- simCtreeFur(Y, 5, freq = freq, h = h)
  pred_rf_lag <- simRFLag(Y, freq = freq, h = h)
  pred_rf_fur <- simRFFur(Y, 3, freq = freq, h = h)
  pred_ese <- simEsExtra(Y, freq = freq, h = h)

  return(list(RPART1 = pred_rpart_2,
              CTREE_LAG = pred_ctree_lag,
              CTREE_FUR = pred_ctree_fur,
              RF_LAG = pred_rf_lag,
              RF_FUR = pred_rf_fur,
              ESE = pred_ese))
}
