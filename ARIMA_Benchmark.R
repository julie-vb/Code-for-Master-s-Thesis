#Code for making the ARIMA(2,1,2) and evaluating its predictive performance

library(forecast)
y <- read.csv('C:/Users/justj/Documents/Thesis files/State-Space Models/Global dataset/Data/global_rv_quarterly_outliersaveraged.csv')
y <- as.matrix(y[ , 2])
fit2 <- Arima(y, order = c(2,1,2),method='ML')
checkresiduals(fit2)
AIC(fit2)
BIC(fit2)

#Predictions
y <- read.csv('C:/Users/justj/Documents/Thesis files/global_rv_quarterly_outliersaveraged.csv')
y <- as.matrix(y[ , 2])
y_test <- as.matrix(y[125:172])
y <- as.matrix(y[1:124]) 
h <- nrow(y_test)
arima_212 <- Arima(y, order = c(2,1,2), method = 'ML')
forecast_arima <- forecast(arima_212, h = h)
y_pred <- forecast_arima$mean
rmse <- sqrt(mean((y_test - y_pred)^2))
mae <- mean(abs(y_test - y_pred))
