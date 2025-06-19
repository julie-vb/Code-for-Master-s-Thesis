#Code to build the models for Dataset 2: RV across space

#Required libraries
library(KFAS)
library(arfima)

y <- read.csv('C:/Users/justj/Documents/Thesis files/global_rv_quarterly_outliersaveraged.csv')
y <- as.matrix(y[ , 2])

k <- 1
ar_order <- 2

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

#SSMG1: Quadratic trend, AR2 cycle, seasonal included
buildTC_g1_uv_season <- function(theta, y, ntrend, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- t(matrix(c(2, -1, 1, 0), 2, 2))
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:2]
  ar <- ar[-(1:2)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 2)
  Z_ar[1, 1] <- 1
  
  
  T_season <- matrix(0, 3, 3)
  T_season[cbind(2:3, 1:2)] <- 1
  T_season[1, 1:3] <- -1
  
  T_seas <- T_season
  R_seas <- diag(3)
  Z_seas <- matrix(0, nrow = k, ncol = 3)
  Z_seas[1,1] <- 1
  
  
  Tt <- bdiag(T_trend, T_ar, T_seas) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))  
  Rt <- rbind(bdiag(c(1, 0), R_ar), matrix(0, nrow = 3*k, ncol = ncol(Qt))) 
  
  
  Z_trend <- matrix(c(1, 0), nrow = k, ncol = 2)
  Zt <- cbind(Z_trend, Z_ar, Z_seas)
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  
  SSM$P1inf[-(1: (k + NROW(T_trend) + NROW(T_ar))) , -(1: (k + NROW(T_trend) + NROW(T_ar)))] <- diag(NROW(T_seas)) 
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}               

#SSMG2: Random walk with deterministic drift trend, AR2 cycle, seasonal included
buildTC_g2_uv_season <- function(theta, y, ntrend, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- t(matrix(c(1, 1, 0, 1), 2, 2))
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:2]
  ar <- ar[-(1:2)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 2)
  Z_ar[1, 1] <- 1
  
  
  T_season <- matrix(0, 3, 3)
  T_season[cbind(2:3, 1:2)] <- 1
  T_season[1, 1:3] <- -1
  
  T_seas <- T_season
  R_seas <- diag(3)
  Z_seas <- matrix(0, nrow = k, ncol = 3)
  Z_seas[1,1] <- 1
  
  
  Tt <- bdiag(T_trend, T_ar, T_seas) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))  
  Rt <- rbind(bdiag(c(1, 0), R_ar), matrix(0, nrow = 3*k, ncol = ncol(Qt))) 
  
  
  Z_trend <- matrix(c(1, 0), nrow = k, ncol = 2)
  Zt <- cbind(Z_trend, Z_ar, Z_seas)
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  SSM$P1inf[(k+2) , (k+2)] <- 1
  SSM$P1inf[-(1: (k + NROW(T_trend) + NROW(T_ar))) , -(1: (k + NROW(T_trend) + NROW(T_ar)))] <- diag(NROW(T_seas)) 
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG3: Random walk without drift, AR2 cycle, seasonal included
buildTC_g3_uv_season <- function(theta, y, ntrend, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- matrix(c(1), 1, 1)
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:2]
  ar <- ar[-(1:2)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 2)
  Z_ar[1, 1] <- 1
  
  
  T_season <- matrix(0, 3, 3)
  T_season[cbind(2:3, 1:2)] <- 1
  T_season[1, 1:3] <- -1
  
  T_seas <- T_season
  R_seas <- diag(3)
  Z_seas <- matrix(0, nrow = k, ncol = 3)
  Z_seas[1,1] <- 1
  
  
  Tt <- bdiag(T_trend, T_ar, T_seas) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))  
  Rt <- rbind(bdiag(c(1), R_ar), matrix(0, nrow = 3*k, ncol = ncol(Qt))) 
  
  
  Z_trend <- matrix(c(1), nrow = k, ncol = 1)
  Zt <- cbind(Z_trend, Z_ar, Z_seas)
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  
  SSM$P1inf[-(1: (k + NROW(T_trend) + NROW(T_ar))) , -(1: (k + NROW(T_trend) + NROW(T_ar)))] <- diag(NROW(T_seas)) 
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG4: No trend, AR2 cycle, seasonal included
buildTC_g4_uv_season <- function(theta, y, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:2]
  if (any(abs(toComp(ar_coef)$eigv) >= 0.99)) return(1e10)
  ar <- ar[-(1:2)]
  
  #Added following line because the optimisation gave an error of non-finite values
  #being produced
  if (!is.finite(q_cycle) || exp(q_cycle) > 1e6 || exp(q_cycle) < 1e-8) return(1e10)
  
  T_ar <- genT_ar(ar_coef)
  R_ar <- matrix(c(1, 0),ncol=1)
  Z_ar <- matrix(0, nrow = k, ncol = 2)
  Z_ar[1, 1] <- 1
  
  
  T_season <- matrix(0, 3, 3)
  T_season[cbind(2:3, 1:2)] <- 1
  T_season[1, 1:3] <- -1
  
  T_seas <- T_season
  R_seas <- diag(3)
  Z_seas <- matrix(0, nrow = k, ncol = 3)
  Z_seas[1,1] <- 1
  
  
  Tt <- bdiag(T_ar, T_seas)
  Qt <- matrix(c(exp(q_cycle)),nrow=1,ncol=1)
  Rt <- rbind(R_ar, matrix(0, nrow = 3*k, ncol = ncol(Qt))) 
  
  Zt <- cbind(Z_ar, Z_seas)
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  
  SSM$P1inf[-(1: (k + NROW(T_ar))) , -(1: (k + NROW(T_ar)))] <- diag(NROW(T_seas)) 
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}               

#SSMG5: Quadratic trend, AR2 cycle, no seasonal
buildTC_g5_uv_noseason <- function(theta, y, ntrend, ret.ll = FALSE, ar_order) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- t(matrix(c(2, -1, 1, 0), 2, 2))
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:2]
  ar <- ar[-(1:2)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 2)
  Z_ar[1, 1] <- 1
  
  
  Tt <- bdiag(T_trend, T_ar) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))  
  Rt <- rbind(bdiag(c(1, 0), R_ar)) 
            
              
  Z_trend <- matrix(c(1, 0), nrow = k, ncol = 2)
  Zt <- cbind(Z_trend, Z_ar)
              
              
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
              
              
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG6: Random walk with deterministic drift trend, AR2 cycle, no seasonal
buildTC_g6_uv_noseason <- function(theta, y, ntrend, ret.ll = FALSE, ar_order) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- t(matrix(c(1, 1, 0, 1), 2, 2))
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:2]
  ar <- ar[-(1:2)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 2)
  Z_ar[1, 1] <- 1
  
  
  Tt <- bdiag(T_trend, T_ar) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))  
  Rt <- rbind(bdiag(c(1, 0), R_ar)) 
  
  
  Z_trend <- matrix(c(1, 0), nrow = k, ncol = 2)
  Zt <- cbind(Z_trend, Z_ar)
  

  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  SSM$P1inf[(k+2) , (k+2)] <- 1
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

SSMG7: Random walk without drift trend, AR2 cycle, no seasonal
buildTC_g7_uv_noseason <- function(theta, y, ntrend, ret.ll = FALSE, ar_order) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- t(matrix(c(1), 1, 1))
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:2]
  ar <- ar[-(1:2)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 2)
  Z_ar[1, 1] <- 1
  
  
  Tt <- bdiag(T_trend, T_ar) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))  
  Rt <- rbind(bdiag(c(1), R_ar)) 
  
  
  Z_trend <- matrix(c(1), nrow = k, ncol = 1)
  Zt <- cbind(Z_trend, Z_ar)
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG8: No trend, AR2 cycle, no seasonal
buildTC_g8_uv_noseason <- function(theta, y, ret.ll = FALSE, ar_order) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:2]
  #if (!toComp(ar_coef)$stable) return(1e10)
  ar_coef <- tryCatch(arfima::PacfToAR(ar[1:2]), error = function(e) return(rep(NA, 2)))
  #ar_coef <- ar[1:2]
  ar <- ar[-(1:2)]
  
  #Added following line because the optimisation gave an error of non-finite values
  #being produced
  #if (!is.finite(q_cycle) || exp(q_cycle) > 1e6 || exp(q_cycle) < 1e-8) return(1e10)
  
  T_ar <- genT_ar(ar_coef)
  R_ar <- matrix(c(1, 0),ncol=1)
  Q_ar <- matrix(c(exp(q_cycle)),nrow=1,ncol=1)
  Z_ar <- matrix(0, nrow = k, ncol = 2)
  Z_ar[1, 1] <- 1
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Z_ar, T = T_ar, R = R_ar, Q = Q_ar), H = diag(k) * 0)
  
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG9: Quadratic trend, AR4 cycle, seasonal included
buildTC_9_uv_season <- function(theta, y, ntrend, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- t(matrix(c(2, -1, 1, 0), 2, 2))
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:4]
  #Added because ar coefs became instable during optimisation:
  if (!toComp(ar_coef)$stable) return(1e10) 
  ar <- ar[-(1:4)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0, 0, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 4)
  Z_ar[1, 1] <- 1
  
  
  T_season <- matrix(0, 3, 3)
  T_season[cbind(2:3, 1:2)] <- 1
  T_season[1, 1:3] <- -1
  
  T_seas <- T_season
  R_seas <- diag(3)
  Z_seas <- matrix(0, nrow = k, ncol = 3)
  Z_seas[1,1] <- 1
  
  
  Tt <- bdiag(T_trend, T_ar, T_seas) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))  
  Rt <- rbind(bdiag(c(1, 0), R_ar), matrix(0, nrow = 3*k, ncol = ncol(Qt))) 
  
  
  Z_trend <- matrix(c(1, 0), nrow = k, ncol = 2)
  Zt <- cbind(Z_trend, Z_ar, Z_seas)
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  
  SSM$P1inf[-(1: (k + NROW(T_trend) + NROW(T_ar))) , -(1: (k + NROW(T_trend) + NROW(T_ar)))] <- diag(NROW(T_seas)) 
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG10: Random walk with deterministic drift trend, AR4 cycle, seasonal included
buildTC_g10_uv_season <- function(theta, y, ntrend, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- t(matrix(c(1, 1, 0, 1), 2, 2))
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:4]
  #Added because ar coefs became instable during optimisation:
  if (!toComp(ar_coef)$stable) return(1e10) 
  ar <- ar[-(1:4)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0, 0, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 4)
  Z_ar[1, 1] <- 1
  
  
  T_season <- matrix(0, 3, 3)
  T_season[cbind(2:3, 1:2)] <- 1
  T_season[1, 1:3] <- -1
  
  T_seas <- T_season
  R_seas <- diag(3)
  Z_seas <- matrix(0, nrow = k, ncol = 3)
  Z_seas[1,1] <- 1
  
  
  Tt <- bdiag(T_trend, T_ar, T_seas) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))  
  Rt <- rbind(bdiag(c(1, 0), R_ar), matrix(0, nrow = 3*k, ncol = ncol(Qt))) 
  
  
  Z_trend <- matrix(c(1, 0), nrow = k, ncol = 2)
  Zt <- cbind(Z_trend, Z_ar, Z_seas)
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  SSM$P1inf[(k+2) , (k+2)] <- 1
  SSM$P1inf[-(1: (k + NROW(T_trend) + NROW(T_ar))) , -(1: (k + NROW(T_trend) + NROW(T_ar)))] <- diag(NROW(T_seas)) 
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG11: Random walk without drift trend, AR4 cycle, seasonal included
buildTC_g11_uv_season <- function(theta, y, ntrend, ret.ll = FALSE, ar_order) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- matrix(c(1), 1, 1)
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:4]
  #Added because ar coefs became instable during optimisation:
  if (!toComp(ar_coef)$stable) return(1e10) 
  ar <- ar[-(1:4)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0, 0, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 4)
  Z_ar[1, 1] <- 1
  
  Tt <- bdiag(T_trend, T_ar) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))  
  Rt <- bdiag(c(1), R_ar) 
  
  
  Z_trend <- matrix(c(1), nrow = k, ncol = 1)
  Zt <- cbind(Z_trend, Z_ar)
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG12: No trend, AR4 cycle, seasonal included
buildTC_g12_uv_season <- function(theta, y, ret.ll = FALSE, ar_order, nseas = 3) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:4]
  ar <- ar[-(1:4)]
  
  #Added following line because the optimisation gave an error of non-finite values
  #being produced
  if (!is.finite(q_cycle) || exp(q_cycle) > 1e6 || exp(q_cycle) < 1e-8) return(1e10)
  
  T_ar <- genT_ar(ar_coef)
  R_ar <- matrix(c(1, 0, 0, 0),nrow=4,ncol=1)
  Z_ar <- matrix(0, nrow = k, ncol = 4)
  Z_ar[1, 1] <- 1
  
  
  
  T_season <- matrix(0, 3, 3)
  T_season[cbind(2:3, 1:2)] <- 1
  T_season[1, 1:3] <- -1
  
  T_seas <- T_season
  R_seas <- diag(3)
  Z_seas <- matrix(0, nrow = k, ncol = 3)
  Z_seas[1,1] <- 1
  
  
  Tt <- bdiag(T_ar, T_seas)
  Qt <- matrix(c(exp(q_cycle)),nrow=1,ncol=1)
  Rt <- rbind(R_ar, matrix(0, nrow = 3*k, ncol = ncol(Qt))) 
  
  Zt <- cbind(Z_ar, Z_seas)

  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  
  SSM$P1inf[-(1: (k + NROW(T_ar))) , -(1: (k + NROW(T_ar)))] <- diag(NROW(T_seas)) 
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG13: Quadratic trend, AR4 cycle, no seasonal
buildTC_13_uv_noseason <- function(theta, y, ntrend, ret.ll = FALSE, ar_order) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- t(matrix(c(2, -1, 1, 0), 2, 2))
  Z_trend <- matrix(c(1,0),nrow=k,ncol=2)
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:4]
  #Added because ar coefs became instable during optimisation:
  if (!toComp(ar_coef)$stable) return(1e10) 
  ar <- ar[-(1:4)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- matrix(c(1, 0, 0, 0),nrow=4,ncol=1)
  Z_ar <- matrix(0, nrow = k, ncol = 4)
  Z_ar[1, 1] <- 1
  
  Tt <- bdiag(T_trend, T_ar)
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))
  Rt <- bdiag(c(1,0),R_ar)
  Zt <- cbind(Z_trend,Z_ar)
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  

  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG14: Random walk with deterministic drift, AR4 cycle, no seasonal
buildTC_14_uv_noseason <- function(theta, y, ntrend, ret.ll = FALSE, ar_order) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- t(matrix(c(1, 1, 0, 1), 2, 2))
  Z_trend <- matrix(c(1,0),nrow=k,ncol=2)
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:4]
  #Added because ar coefs became instable during optimisation:
  if (!toComp(ar_coef)$stable) return(1e10) 
  ar <- ar[-(1:4)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0, 0, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 4)
  Z_ar[1, 1] <- 1
  
  Tt <- bdiag(T_trend, T_ar)
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))
  Rt <- bdiag(c(1,0),R_ar)
  Zt <- cbind(Z_trend,Z_ar)
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  SSM$P1inf[(k+2) , (k+2)] <- 1
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG15: Random walk without drift trend, AR4 cycle, no seasonal
buildTC_g15_uv_noseason <- function(theta, y, ntrend, ret.ll = FALSE, ar_order) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_trend <- theta[1]
  theta <- theta[-1]
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  T_trend <- matrix(c(1), 1, 1)
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:4]
  #Added because ar coefs became instable during optimisation:
  if (!toComp(ar_coef)$stable) return(1e10) 
  ar <- ar[-(1:4)]
  T_ar <- genT_ar(ar_coef)
  R_ar <- c(1, 0, 0, 0)
  Z_ar <- matrix(0, nrow = k, ncol = 4)
  Z_ar[1, 1] <- 1
  
  Tt <- bdiag(T_trend, T_ar) 
  Qt <- diag(c(exp(q_trend), exp(q_cycle)))  
  Rt <- bdiag(c(1), R_ar) 
  
  
  Z_trend <- matrix(c(1), nrow = k, ncol = 1)
  Zt <- cbind(Z_trend, Z_ar)
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt), H = diag(k) * 0)
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}

#SSMG16: No trend, AR4 cycle, no seasonal
buildTC_g16_uv_noseason <- function(theta, y, ret.ll = FALSE, ar_order) {
  n <- nrow(y)
  k <- ncol(y)
  ncycle <- 1
  
  q_cycle <- theta[1]
  theta <- theta[-1]
  
  ar <- theta
  
  genT_ar <- function(ar) toComp(ar)$CompMat
  
  ar_coef <- ar[1:4]
  ar <- ar[-(1:4)]
  
  #Added following line because the optimisation gave an error of non-finite values
  #being produced
  #if (!is.finite(q_cycle) || exp(q_cycle) > 1e6 || exp(q_cycle) < 1e-8) return(1e10)
  if (!toComp(ar_coef)$stable) return(1e10)
  
  T_ar <- genT_ar(ar_coef)
  R_ar <- matrix(c(1, 0, 0, 0),nrow=4,ncol=1)
  Q_ar <- matrix(c(exp(q_cycle)),nrow=1,ncol=1)
  Z_ar <- matrix(0, nrow = k, ncol = 4)
  Z_ar[1, 1] <- 1
  
  
  SSM <- SSModel(y ~ SSMcustom(Z = Z_ar, T = T_ar, R = R_ar, Q = Q_ar), H = diag(k) * 0)
  
  
  if (ret.ll) {
    return(-logLik(SSM))
  } else {
    return(SSM)
  }
}
