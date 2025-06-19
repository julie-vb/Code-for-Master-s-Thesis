#Code for preliminary analysis Dataset 2

library(ggplot2)
library(zoo)
library(dplyr)

data <- read.csv('C:/Users/justj/Documents/Thesis files/global_rv_quarterly_outliersaveraged.csv')
data$YearQuarter <- as.yearqtr(data$YearQuarter, format = '%Y Q%q')
ggplot(data = data, aes(x = YearQuarter, y = realised_volatility)) +
  geom_line() +
  labs(title = 'Realised Volatility Across Space (Quarterly)', x = 'Quarter', y = 'Realised Volatility') +
  theme_minimal()

ggplot(data = data, aes(x = YearQuarter, y = log(realised_volatility))) +
  geom_line() +
  labs(title = 'Log of Realised Volatility Across Space', x = 'Quarter', y = 'Log of Realised Volatility') +
  theme_minimal()

summary(data$realised_volatility)
sd(data$realised_volatility)

top_quarters <- data %>%
  arrange(desc(realised_volatility)) %>%
  head(25)  
print(top_quarters)
#Q1:4; Q2:13; Q3:3; Q4:5.

bottom_quarters <- data %>%
  arrange(realised_volatility) %>%
  head(50)  
print(bottom_quarters)
#Q1:0; Q2:0; Q3:25; Q4:0
#Top 50 is only Q1, Q3, and one time Q4.

#Stationarity tests
library(urca)
data <- read.csv('C:/Users/justj/Documents/Thesis files/State-Space Models/Global dataset/Data/global_rv_quarterly_outliersaveraged.csv')
data <- as.matrix(data[ , 2])
adf_drift <- ur.df(data, type = 'drift')  
adf_trend <- ur.df(data, type = 'trend')
summary(adf_drift)
summary(adf_trend)
kpss_level <- ur.kpss(data, type = "mu") 
summary(kpss_level)
kpss_trend <- ur.kpss(data, type = "tau") 
summary(kpss_trend)

kpss_pval_interp <- function(statistic, type = 'mu') {
  if (type == 'mu') {
    crit_vals <- c(0.347, 0.463, 0.574, 0.739)
    sig_levels <- c(0.10, 0.05, 0.025, 0.01)
  } else if (type == 'tau') {
    crit_vals <- c(0.119, 0.146, 0.176, 0.216)
    sig_levels <- c(0.10, 0.05, 0.025, 0.01)
  } else {
    stop('Invalid KPSS type. Use 'mu' or 'tau'.')
  }
  
  if (statistic < crit_vals[1]) return(1.0)
  if (statistic > crit_vals[length(crit_vals)]) return(0.001)

  approx(x = crit_vals, y = sig_levels, xout = statistic, rule = 2)$y
}

kpss_mu_stats <- numeric(ncol(data))
kpss_tau_stats <- numeric(ncol(data))

for (i in seq_len(ncol(data))) {
  kpss_mu_stats[i] <- ur.kpss(data[, i], type = 'mu')@teststat
  kpss_tau_stats[i] <- ur.kpss(data[, i], type = 'tau')@teststat
}

# Get interpolated p-values
pvals_mu <- sapply(kpss_mu_stats, kpss_pval_interp, type = 'mu')
pvals_tau <- sapply(kpss_tau_stats, kpss_pval_interp, type = 'tau')

pvals <- c(pvals_mu,pvals_tau,2.2e-16,2.2e-16)
p.adjust(pvals,method='BH')
