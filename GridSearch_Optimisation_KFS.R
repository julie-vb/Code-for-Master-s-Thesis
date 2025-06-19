#Example code for the grid search and optimisation: SSMEC3 (unique RW trends, unique cycles, unique seasonals).
#Subsequent application of the Kalman Filter and Smoother.

#Required libraries
library(arfima)
library(KFAS)

#Grid search for good starting values for the parameter estimates
set.seed(42)
R <- 100  #Number of random samples
p <- 2   #AR order

START.val <- matrix(NA, nrow = R, ncol = k + k + p * k)
ll_vals <- rep(Inf, R)
cat('Dimensions of START.val: ', dim(START.val), '\n') 

for (i in 1:R) {
  
  log_q_trend <- runif(1, -25, -5)
  log_q_cycle <- runif(1, -10, 0)
  
  as <- matrix(NA, k, p)  
  as[, 1] <- runif(k, 0.5, 1)   
  
  if (p > 1) {
    as[, 2:p] <- runif(k * (p - 1), -0.7, 0.7)  
  }
  
  ar_mat <- matrix(NA, k, p)  
  for (city in 1:k) {
    ar_mat[city, ] <- arfima::PacfToAR(as[city, ]) 
  }
  apply(ar_mat, 1, function(x) max(abs(toComp(x)$eigv)))
  
  theta_i <- c(                
    rep(log_q_trend, k),                   
    rep(log_q_cycle, k),         
    as.vector(t(ar_mat))              
  )
  
  START.val[i, ] <- theta_i  
}
        
thetas <- vector('list', R)
ll_vals <- rep(Inf, R)

pb <- txtProgressBar(min = 0, max = R, style = 3)

for (i in 1:R) {
  theta_try <- START.val[i, ]
  
  ll_vals[i] <- tryCatch(
    buildTC_2_mv_season(theta_try, y, ar_order = rep(p, k), ret.ll = TRUE),
    error = function(e) Inf
  )
  
  thetas[[i]] <- theta_try  
  setTxtProgressBar(pb, i)
}
close(pb)

#Extracting the best initial values
best_index <- which.min(ll_vals)
theta_best <- thetas[[best_index]]

#True theta
buildTC_2_mv_season(theta = theta_best, y, ar_order = rep(p, k), ret.ll = TRUE)
ll_vals[best_index]

#True theta
buildTC_2_mv_season(theta = theta_best, y, ar_order = rep(p, k), ret.ll = TRUE)
ll_vals[best_index]

#Optimisation
est <- optim(
  theta_best,
  buildTC_2_mv_season,
  y = y,
  ret.ll=TRUE, 
  ar_order = ar_order,
  method = 'BFGS',
  control = list(maxit = 1000)
)

new_theta <- est$par
final_ssm <- buildTC_2_mv_season(
  theta = new_theta,
  y = y,          
  ret.ll = FALSE,       
  ar_order = ar_order,  
  nseas = 3             
)

#Applying KFS
filtered_ssm <- KFAS::KFS(final_ssm)
filtered_states <- filtered_ssm$alphahat
smoothed <- KFS(model_opt, smoothing = 'state', filtering = 'state')

