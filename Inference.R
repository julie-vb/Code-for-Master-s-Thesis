#Example code for inference performed on SSMEC3

library(arfima)
library(KFAS)
library(numDeriv)

#Optimised parameters
theta <- readRDS('new_theta_uRW_uc_us.rds')
negll_group1 <- function(theta) {
  buildTC_2_mv_season(
    theta = theta,
    y = y,           
    ret.ll = TRUE,
    ar_order = rep(2, k),     
    nseas = 3
  )
}
H <- numDeriv::hessian(func = negll_group1, x = theta)
J <- numDeriv::jacobian(function(x) (exp(x)), c(q_trend, q_cycle))
cov0 <- solve(H)
Jacobian <- bdiag(J,diag(length(ar_coefs)))
Cov <- Jacobian %*% cov0 %*% t(Jacobian)
theta.trans <- c(exp(q_trend),exp(q_cycle),ar_coefs)
se.trans <- sqrt(abs(diag(Cov)))
EST <- cbind(theta.trans, se.trans) 
colnames(EST) <- c('par', 'se') 
alpha <- 0.05  
z <- qnorm(1 - alpha/2)  
lower <- theta.trans - z*se.trans
upper <- theta.trans + z*se.trans

z_scores <- theta.trans/se.trans
p_values <- 2*pnorm(-abs(z_scores))
p_adj <- p.adjust(p_values, method = 'BH')

CI <- cbind(lower, upper)
EST <- cbind(
  Estimate = theta.trans,
  SE = se.trans,
  Lower95 = lower,
  Upper95 = upper
)

m <- 148
alpha <- 0.05
alpha_bonf <- alpha/m
z_bonf <- qnorm(1 - alpha_bonf/2)

lower_bonf <- theta.trans - z_bonf*se.trans
upper_bonf <- theta.trans + z_bonf*se.trans

CI_bonf <- cbind(lower_bonf, upper_bonf)
EST_bonf <- cbind(
  Estimate = theta.trans,
  SE = se.trans,
  Lower95_Bonf = lower_bonf,
  Upper95_Bonf = upper_bonf
)
