#Example code for the predictive performance of Dataset 1 models: SSMEC3

#Required libraries
library(arfima)
library(KFAS)

#Dividing data into training and testing sets
y <- read.csv('C:/Users/justj/Documents/Thesis files/realised_volatility_data.csv', row.names=1)
y <- as.matrix(y)
y_test <- as.matrix(y[241:298,])
y <- as.matrix(y[1:240,]) 
h <- nrow(y_test)

#Theta estimated using only the training set
new_theta <- readRDS('C:/Users/justj/Documents/Thesis files/Thetas/Thetas prediction/theta_training_commonRW_UCUS.rds')

k <- ncol(y)
ar_order <- rep(2, k)

#Two functions necessary for the model building, functions are from Hartl (2025)
#ToComp function: stability check
toComp <- function(x)
{
  if (is.null(dim(x))) x <- matrix(x,1)
  k <- nrow(x)
  p <- ncol(x)/k
  if (p>1) x <- rbind(x, cbind(diag((p-1)*k),matrix(0,(p-1)*k,k)))
  eigv <- eigen(x,only.values=TRUE)$values
  return(list(CompMat=x,eigv=eigv,stable=all(round(abs(eigv), digits = 2)<1)))
}

#Bdiag: building block diagonal matrices
bdiag <- function (...) 
{
  if (nargs() == 1) 
    x <- as.list(...)
  else x <- list(...)
  n <- length(x)
  if (n == 0) 
    return(NULL)
  x <- lapply(x, function(y) if (length(y)) 
    as.matrix(y)
    else stop("Zero-length component in x"))
  d <- array(unlist(lapply(x, dim)), c(2, n))
  rr <- d[1, ]
  cc <- d[2, ]
  rsum <- sum(rr)
  csum <- sum(cc)
  out <- array(0, c(rsum, csum))
  ind <- array(0, c(4, n))
  rcum <- cumsum(rr)
  ccum <- cumsum(cc)
  ind[1, -1] <- rcum[-n]
  ind[2, ] <- rcum
  ind[3, -1] <- ccum[-n]
  ind[4, ] <- ccum
  imat <- array(1:(rsum * csum), c(rsum, csum))
  iuse <- apply(ind, 2, function(y, imat) imat[(y[1] + 1):y[2], 
                                               (y[3] + 1):y[4]], imat = imat)
  iuse <- as.vector(unlist(iuse))
  out[iuse] <- unlist(x)
  return(out)
}

#Model building function
buildTC_3_mv_season <- function(theta, y, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- length(ar_order)
  
  q_trend <- theta[1:k]
  theta <- theta[-(1:k)]
  q_cycle <- theta[1:ncycle]
  theta <- theta[-(1:ncycle)]
  
  ar <- theta
  
  T_trend <- matrix(c(1), 1, 1) 
  R_trend <- matrix(c(1), ncol = 1)                         
  Z_trend <- matrix(0, nrow = k, ncol = 1)   
  Z_trend[1, 1] <- 1                         
  
  for (i in 2:k) {
    T_trend <- bdiag(T_trend, matrix(c(1), 1, 1))
    R_trend <- bdiag(R_trend, matrix(c(1), ncol = 1))
    Z_trend_block <- matrix(0, nrow = k, ncol = 1)
    Z_trend_block[i, 1] <- 1
    Z_trend <- cbind(Z_trend, Z_trend_block)
  }
  
  genT_ar <- function(ar) toComp(ar)$CompMat 
  
  p_ar <- ar_order[1] 
  ar_order <- ar_order[-1]
  ar_coef <- ar[1:p_ar]
  ar <- ar[-(1:p_ar)]
  
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, rep(0, p_ar -1))
  Z_ar <- matrix(0, nrow = k, ncol = p_ar)
  Z_ar[1, 1] <- 1
  
  for (i in 2:k) {
    p_ar <- ar_order[1] 
    ar_order <- ar_order[-1]
    ar_coef <- ar[1:p_ar]
    ar <- ar[-(1:p_ar)]
    
    T_ar <- bdiag(T_ar, genT_ar(ar_coef))
    R_ar <- bdiag(R_ar, c(1, rep(0, p_ar -1)))
    Z_ar_1 <- matrix(0, nrow = k, ncol = p_ar)
    Z_ar_1[i, 1] <- 1
    Z_ar <- cbind(Z_ar, Z_ar_1)
  }
  
  T_season <- matrix(0, 3, 3)
  T_season[cbind(2:3, 1:2)] <- 1
  T_season[1, 1:3] <- -1
  
  T_seas <- T_season
  Z_seas <- matrix(0, nrow = k, ncol = 3)
  Z_seas[1,1] <- 1
  
  for (i in 2:k) {
    T_seas <- bdiag(T_seas, T_season) 
    Z_seas_1 <- matrix(0, nrow = k, ncol = 3) 
    Z_seas_1[i, 1] <- 1 
    Z_seas <- cbind(Z_seas, Z_seas_1)
  }
  
  Tt <- bdiag(T_trend, T_ar, T_seas) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))
  Rt <- rbind(bdiag(R_trend, R_ar), matrix(0, nrow = 3*k, ncol = ncol(Qt))) 
  Zt <- cbind(Z_trend, Z_ar, Z_seas)
  

  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  
  SSM$P1inf[-(1: (k + NROW(T_trend) + NROW(T_ar))) , -(1: (k + NROW(T_trend) + NROW(T_ar)))] <- diag(NROW(T_seas)) 
  
  if (ret.ll) {
    ll <- logLik(SSM)
    cat('Log likelihood: ', ll, '\n')
    return(-ll)
  } else {
    return(SSM)
  }
}


model_train <- buildTC_3_mv_season(
  theta = new_theta,
  y = y,
  ret.ll = FALSE,
  ar_order
)

y_future <- matrix(NA, nrow = h, ncol = ncol(y))
model_future <- buildTC_3_mv_season(
  theta = new_theta,
  y = y_future,
  ret.ll = FALSE,
  ar_order
)

y_pred <- predict(model_train, newdata = model_future)

city_names <- names(y_pred)
y_pred_mat <- do.call(cbind, lapply(city_names, function(city) {
  y_pred[[city]][, 'fit']
}))
colnames(y_pred_mat) <- city_names

rmse <- sqrt(colMeans((y_test - y_pred_mat)^2, na.rm = TRUE)) 
mae <- colMeans(abs(y_test - y_pred_mat), na.rm = TRUE)  
overall_rmse <- sqrt(mean((y_test - y_pred_mat)^2, na.rm = TRUE))
overall_mae <- mean(abs(y_test - y_pred_mat), na.rm = TRUE)
  
