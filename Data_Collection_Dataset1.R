library(readr)
library(raster)
library(ncdf4)
library(ggplot2)
library(dplyr)

#Data from the E-OBS dataset (v30.0e) from Copernicus Climate Change Service
nc <- nc_open('C:/Users/justj/Downloads/rr_ens_mean_0.1deg_reg_v30.0e.nc')
lon <- ncvar_get(nc, 'longitude')  
lat <- ncvar_get(nc, 'latitude')   
time <- ncvar_get(nc, 'time')  

time_units <- ncatt_get(nc, 'time', 'units')$value  
time_origin <- sub('.*since ', '', time_units)      
time_readable <- as.POSIXct(time * 86400, origin = time_origin, tz = 'UTC')

#Obtaining daily precipitation values for every city & computing the realised volatility over time using a global mean

#Athens
desired_lon <- 23.72  
desired_lat <- 37.98 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Athens', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#per quarter

desired_lon <- 23.72  
desired_lat <- 37.98  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$Year <- format(prec_df$time, '%Y')
prec_df$Quarter <- paste0('Q', ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_volatility <- prec_df %>%
  group_by(Year, Quarter) %>%
  summarize(realized_volatility = mean((precipitation - global_mean)^2, na.rm = TRUE))
ggplot(quarterly_volatility, aes(x = as.numeric(Year), y = log(realized_volatility), 
                                 color = Quarter, group = Quarter)) +
  geom_line(size = 1) +        
  scale_color_manual(values = c('Q1' = 'blue', 'Q2' = 'red', 'Q3' = 'green', 'Q4' = 'gold')) + 
  labs(title = 'Realized Volatility of Precipitation per Quarter in Athens',
       x = 'Year',
       y = 'Log of Realized Volatility',
       color = 'Quarter') +  
  theme_minimal()



#Barcelona
desired_lon <- 2.16  
desired_lat <- 41.39  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Barcelona', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Belgrade
desired_lon <- 20.47  
desired_lat <- 44.80 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Belgrade', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Berlin
desired_lon <- 13.4000  
desired_lat <- 52.4700 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Berlin', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Bern
desired_lon <- 7.45  
desired_lat <- 46.95  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Bern', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Brno-Turany
desired_lon <- 16.7000  
desired_lat <- 49.1600
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Brno-Turany', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Bucharest
desired_lon <- 26.10  
desired_lat <- 44.44  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Bucharest', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Budapest
desired_lon <- 19.03  
desired_lat <- 47.52  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Budapest', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Copenhagen
desired_lon <- 12.57  
desired_lat <- 55.68 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Copenhagen', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#De_Bilt
desired_lon <- 5.30  
desired_lat <- 52.17 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in De Bilt', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#per quarter
desired_lon <- 5.30  
desired_lat <- 52.17
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$Year <- format(prec_df$time, '%Y')
prec_df$Quarter <- paste0('Q', ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_volatility <- prec_df %>%
  group_by(Year, Quarter) %>%
  summarize(realized_volatility = mean((precipitation - global_mean)^2, na.rm = TRUE))
ggplot(quarterly_volatility, aes(x = as.numeric(Year), y = log(realized_volatility), 
                                 color = Quarter, group = Quarter)) +
  geom_line(size = 1) +        
  scale_color_manual(values = c('Q1' = 'blue', 'Q2' = 'red', 'Q3' = 'green', 'Q4' = 'gold')) + 
  labs(title = 'Realized Volatility of Precipitation per Quarter in De Bilt',
       x = 'Year',
       y = 'Log of Realized Volatility',
       color = 'Quarter') +  
  theme_minimal()


#Dublin
desired_lon <- -6.27 
desired_lat <- 53.35
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Dublin', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Helsinki
desired_lon <- 24.94  
desired_lat <- 60.21 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Helsinki', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Hohenpeissenberg
desired_lon <- 11.02  
desired_lat <- 47.80
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Hohenpeissenberg', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Innsbruck
desired_lon <- 11.385  
desired_lat <- 47.26 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Innsbruck', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Karlsruhe
desired_lon <- 8.33 
desired_lat <- 48.97  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Karlsruhe', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Kremsmünster
desired_lon <- 14.13  
desired_lat <- 48.055 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Kremsmünster', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#per quarter
desired_lon <- 14.13 
desired_lat <- 48.055
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$Year <- format(prec_df$time, '%Y')
prec_df$Quarter <- paste0('Q', ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_volatility <- prec_df %>%
  group_by(Year, Quarter) %>%
  summarize(realized_volatility = mean((precipitation - global_mean)^2, na.rm = TRUE))
ggplot(quarterly_volatility, aes(x = as.numeric(Year), y = log(realized_volatility), 
                                 color = Quarter, group = Quarter)) +
  geom_line(size = 1) +        
  scale_color_manual(values = c('Q1' = 'blue', 'Q2' = 'red', 'Q3' = 'green', 'Q4' = 'gold')) + 
  labs(title = 'Realized Volatility of Precipitation per Quarter in Kremsmünster',
       x = 'Year',
       y = 'Log of Realized Volatility',
       color = 'Quarter') +  
  theme_minimal()


#Lisbon
desired_lon <- -9.14
desired_lat <- 38.71 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Lisbon', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#London
desired_lon <- -0.12  
desired_lat <- 51.51
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in London', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Lviv
desired_lon <- 24.03  
desired_lat <- 49.84  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Lviv', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Madrid
desired_lon <- -3.70  
desired_lat <- 40.42 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Madrid', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Marseille
desired_lon <- 5.37  
desired_lat <- 43.30
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Marseille', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Milan
desired_lon <- 9.19  
desired_lat <- 45.47 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Milan', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Munich
desired_lon <- 11.55  
desired_lat <- 48.17 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Munich', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Oslo
desired_lon <- 10.76  
desired_lat <- 59.91 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Oslo', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Paris
desired_lon <- 2.50  
desired_lat <- 48.80 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Paris', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Regensburg
desired_lon <- 12.10  
desired_lat <- 49.05 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Regensburg', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Rome
desired_lon <- 12.50  
desired_lat <- 41.90  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Rome', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Sofia
desired_lon <- 23.32  
desired_lat <- 42.70  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Sofia', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Stockholm
desired_lon <- 18.07  
desired_lat <- 59.32 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Stockholm', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Stuttgart
desired_lon <- 9.20  
desired_lat <- 48.83
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Stuttgart', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Tirana
desired_lon <- 19.82 
desired_lat <- 41.33  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Tirana', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Trondheim
desired_lon <- 10.35 
desired_lat <- 63.35
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Trondheim', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Uppsala
desired_lon <- 17.63  
desired_lat <- 59.87
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Uppsala', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Vienna
desired_lon <- 16.36  
desired_lat <- 48.25 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Vienna', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Vilnius
desired_lon <- 25.10 
desired_lat <- 54.63  
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Vilnius', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Warsaw
desired_lon <- 20.97 
desired_lat <- 52.17 
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Warsaw', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')

#Zagreb
desired_lon <- 15.97  
desired_lat <- 45.82
lon_index <- which.min(abs(lon - desired_lon))
lat_index <- which.min(abs(lat - desired_lat))
precipitation <- ncvar_get(nc, 'rr', start = c(lon_index, lat_index, 1), count = c(1, 1, -1))
prec_df <- data.frame(
  time = time_readable,
  precipitation = precipitation
)
global_mean <- mean(prec_df$precipitation, na.rm = TRUE)
prec_df$YearQuarter <- paste0(format(prec_df$time, '%Y-Q'), ceiling(as.numeric(format(prec_df$time, '%m')) / 3))
quarterly_df_global <- prec_df %>%
  group_by(YearQuarter) %>%
  summarize(across(-time, function(x) mean((x - global_mean)^2, na.rm = TRUE)))
plot(log(quarterly_df_global$precipitation), main = 'Realised Volatility of Precipitation in Zagreb', xlab = 'Quarter', ylab = 'Log of realised volatility', type ='l')


