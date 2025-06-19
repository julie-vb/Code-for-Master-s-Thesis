#Code to build the SSMs for the European city dataset

#Required libraries
library(KFAS)
library(arfima)

y <- read.csv('C:/Users/justj/Documents/Thesis files/realised_volatility_data.csv', row.names=1)
y <- as.matrix(y)

k <- ncol(y)
ar_order <- rep(2, k)

#The following two functions used to build the SSMs are from Hartl (2025).
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
                
#Step 1 - Unique trend, cyclical, and seasonal components; varying trend types.

#SSMEC1: Quadratic trends
buildTC_1_mv_season <- function(theta, y, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- length(ar_order)
  
  q_trend <- theta[1:k]
  theta <- theta[-(1:k)]
  q_cycle <- theta[1:ncycle]
  theta <- theta[-(1:ncycle)]
  
  ar <- theta
  
  T_trend <- t(matrix(c(2, -1, 1, 0), 2, 2))  
  R_trend <- matrix(c(1, 0), ncol = 1)                         
  Z_trend <- matrix(0, nrow = k, ncol = 2)   
  Z_trend[1, 1] <- 1                         
  
  for (i in 2:k) {
    T_trend <- bdiag(T_trend, t(matrix(c(2, -1, 1, 0), 2, 2)))
    R_trend <- bdiag(R_trend, matrix(c(1, 0), ncol = 1))
    Z_trend_block <- matrix(0, nrow = k, ncol = 2)
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

#SSMEC2: Random walk with deterministic drift trends
buildTC_2_mv_season <- function(theta, y, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- length(ar_order)
  
  q_trend <- theta[1:k]
  theta <- theta[-(1:k)]
  q_cycle <- theta[1:ncycle]
  theta <- theta[-(1:ncycle)]
  
  ar <- theta
  
  T_trend <- t(matrix(c(1, 1, 0, 1), 2, 2))  
  R_trend <- matrix(c(1, 0), ncol = 1)                         
  Z_trend <- matrix(0, nrow = k, ncol = 2)   
  Z_trend[1, 1] <- 1                         
  
  for (i in 2:k) {
    T_trend <- bdiag(T_trend, t(matrix(c(1, 1, 0, 1), 2, 2)))
    R_trend <- bdiag(R_trend, matrix(c(1, 0), ncol = 1))
    Z_trend_block <- matrix(0, nrow = k, ncol = 2)
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
  
  SSM$P1inf[(k+2) , (k+2)] <- 1
  SSM$P1inf[(k+4) , (k+4)] <- 1
  SSM$P1inf[(k+6) , (k+6)] <- 1
  SSM$P1inf[(k+8) , (k+8)] <- 1
  SSM$P1inf[(k+10) , (k+10)] <- 1
  SSM$P1inf[(k+12) , (k+12)] <- 1
  SSM$P1inf[(k+14) , (k+14)] <- 1
  SSM$P1inf[(k+16) , (k+16)] <- 1
  SSM$P1inf[(k+18) , (k+18)] <- 1
  SSM$P1inf[(k+20) , (k+20)] <- 1
  SSM$P1inf[(k+22) , (k+22)] <- 1
  SSM$P1inf[(k+24) , (k+24)] <- 1
  SSM$P1inf[(k+26) , (k+26)] <- 1
  SSM$P1inf[(k+28) , (k+28)] <- 1
  SSM$P1inf[(k+30) , (k+30)] <- 1
  SSM$P1inf[(k+32) , (k+32)] <- 1
  SSM$P1inf[(k+34) , (k+34)] <- 1
  SSM$P1inf[(k+36) , (k+36)] <- 1
  SSM$P1inf[(k+38) , (k+38)] <- 1
  SSM$P1inf[(k+40) , (k+40)] <- 1
  SSM$P1inf[(k+42) , (k+42)] <- 1
  SSM$P1inf[(k+44) , (k+44)] <- 1
  SSM$P1inf[(k+46) , (k+46)] <- 1
  SSM$P1inf[(k+48) , (k+48)] <- 1
  SSM$P1inf[(k+50) , (k+50)] <- 1
  SSM$P1inf[(k+52) , (k+52)] <- 1
  SSM$P1inf[(k+54) , (k+54)] <- 1
  SSM$P1inf[(k+56) , (k+56)] <- 1
  SSM$P1inf[(k+58) , (k+58)] <- 1
  SSM$P1inf[(k+60) , (k+60)] <- 1
  SSM$P1inf[(k+62) , (k+62)] <- 1
  SSM$P1inf[(k+64) , (k+64)] <- 1
  SSM$P1inf[(k+66) , (k+66)] <- 1
  SSM$P1inf[(k+68) , (k+68)] <- 1
  SSM$P1inf[(k+70) , (k+70)] <- 1
  SSM$P1inf[(k+72) , (k+72)] <- 1
  SSM$P1inf[(k+74) , (k+74)] <- 1
  
  SSM$P1inf[-(1: (k + NROW(T_trend) + NROW(T_ar))) , -(1: (k + NROW(T_trend) + NROW(T_ar)))] <- diag(NROW(T_seas)) 
  
  if (ret.ll) {
    ll <- logLik(SSM)
    cat('Log likelihood: ', ll, '\n')
    return(-ll)
  } else {
    return(SSM)
  }
}

#SSMEC3: Random walk without drift trends
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
  
  
  print(dim(Tt))
  print(dim(Zt))
  print(dim(Rt))
  print(dim(Qt))
  
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

#SSMEC4: No trend
buildTC_4_mv_season <- function(theta, y, ret.ll = FALSE, ar_order, nseas=3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- length(ar_order)
  
  q_cycle <- theta[1:ncycle]
  theta <- theta[-(1:ncycle)]
  
  ar <- theta
  
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
  
  Tt <- T_ar
  Qt <- diag(exp(q_cycle))
  Rt <- R_ar
  Zt <- Z_ar
  
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
  
  Tt <- bdiag(T_ar, T_seas) 
  Qt <- diag(exp(q_cycle))
  Rt <- rbind(R_ar, matrix(0, nrow = 3*k, ncol = ncol(Qt))) 
  
  Zt <- cbind(Z_ar, Z_seas)
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  SSM$P1inf[-(1: (k + NROW(T_ar))) , -(1: (k + NROW(T_ar)))] <- diag(NROW(T_seas)) 
  
  if (ret.ll) {
    ll <- logLik(SSM)
    cat('Log likelihood: ', ll, '\n')
    return(-ll)
  } else {
    return(SSM)
  }
}


#Step 2: Cycles: unique vs common
#SSMEC5: Unique random walk without drift trends, common cycle, unique seasonals
buildTC_5_mv_season <- function(theta, y, ret.ll = FALSE, ar_order=2, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- length(ar_order)
  
  loadings_cycle <- theta[1:((k - 1) * ncycle)]
  theta <- theta[-(1:((k - 1) * ncycle))]
  
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
  
  ar_coef <- ar  
  if (!toComp(ar_coef)$stable) {
    if (ret.ll) return(1e10)
    else return(NULL)
  }
  T_ar <- genT_ar(ar_coef)
  R_ar <- matrix(c(1,0),2,1)
  Z_ar <- matrix(0, nrow = k, ncol = 2)
  Z_ar[, 1] <- c(1, loadings_cycle)
  
  
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


#Step 3: Trend: common vs unique
#SSMEC6: Common random walk without drift trend, unique cycles, unique seasonals
buildTC_6_mv_season <- function(theta, y, ntrend, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- length(ar_order)
  
  loadings_trend <- theta[1:((k - 1) * ntrend)]
  theta <- theta[-(1:((k - 1) * ntrend))]
  
  q_trend <- theta[1:ntrend]
  theta <- theta[-(1:ntrend)]
  q_cycle <- theta[1:ncycle]
  theta <- theta[-(1:ncycle)]
  
  ar <- theta
  
  T_trend <- matrix(c(1), 1, 1)
  
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
  Rt <- rbind(bdiag(c(1), R_ar), matrix(0, nrow = 3*k, ncol = ncol(Qt))) 
  
  
  Z_trend <- cbind(c(1, loadings_trend), matrix(0, ncol = 2 * (ntrend - 1), nrow = k))
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


#Step 7: Seasonals: common vs unique vs none
#SSMEC7: Common random walk without drift trend, unique cycles, common seasonal
buildTC_7_mv_season <- function(theta, y, ntrend, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- length(ar_order)
  
  loadings_trend <- theta[1:((k - 1) * ntrend)]
  theta <- theta[-(1:((k - 1) * ntrend))]
  
  loadings_season <- theta[1:((k - 1))]
  theta <- theta[-(1:((k - 1)))]
  
  q_trend <- theta[1:ntrend]
  theta <- theta[-(1:ntrend)]
  q_cycle <- theta[1:ncycle]
  theta <- theta[-(1:ncycle)]
  
  ar <- theta
  
  T_trend <- matrix(c(1), 1, 1)
  
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
  Z_seas[-1, 1] <- loadings_season
  
  
  Tt <- bdiag(T_trend, T_ar, T_seas) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle))) 
  Rt <- rbind(bdiag(c(1), R_ar), matrix(0, nrow = 3, ncol = ncol(Qt))) 
  
  
  Z_trend <- cbind(c(1, loadings_trend), matrix(0, ncol = 2 * (ntrend - 1), nrow = k))
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

#SSMEC8: Common random walk without drift trend, unique cycles, no seasonal
  buildTC_8_mv_season <- function(theta, y, ntrend, ret.ll = FALSE, ar_order) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- length(ar_order)
  
  loadings_trend <- theta[1:((k - 1) * ntrend)]
  theta <- theta[-(1:((k - 1) * ntrend))]
  
  
  q_trend <- theta[1:ntrend]
  theta <- theta[-(1:ntrend)]
  q_cycle <- theta[1:ncycle]
  theta <- theta[-(1:ncycle)]
  
  ar <- theta
  
  T_trend <- matrix(c(1), 1, 1)
  
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
  
  
  Tt <- bdiag(T_trend, T_ar) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle))) 
  Rt <- rbind(bdiag(c(1), R_ar)) 

  
  Z_trend <- cbind(c(1, loadings_trend), matrix(0, ncol = 2 * (ntrend - 1), nrow = k))
  Zt <- cbind(Z_trend, Z_ar)
  print(dim(Tt)) #mxm
  print(dim(Zt)) #kxm
  print(dim(Rt)) #mxr
  print(dim(Qt)) #rxr
  
  
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  
  if (ret.ll) {
    system.time(ll <- logLik(SSM))
    cat("Log likelihood: ", ll, "\n")
    return(-ll)
  } else {
    return(SSM)
  }
}
