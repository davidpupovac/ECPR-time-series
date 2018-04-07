###################################################################
###################################################################
#
#                 Times Series Analysis - Day 2
#
###################################################################
###################################################################


###################################################################
#
#                          DECOMPOSITION
#
###################################################################

setwd("D:/My Documents/Desktop/ECPR data") 


# ==============
# Additive decomposition - by hand!!!

# Load data beer.csv
ts_beer<-read.csv(file.choose(),header=TRUE)
ts_beer <- ts(ts_beer, frequency = 4, start=c(1956,1))
plot(ts_beer)

# First, we will estimate a trend by calculating annual moving average of the series
# Have in mind this is quarterly data, so we use order = 4
library(forecast)
trend_beer = ma(ts_beer, order = 4, centre = T)
plot(as.ts(ts_beer))
lines(trend_beer, col="red",lwd=2)

# Second, we detrend the time series by subtracting 
detrend_beer = ts_beer - trend_beer
plot(as.ts(detrend_beer))  # so this is leftover seasonality and error

# Third, calculate average seasonality
# Here we have quarterly data so we have quarterly seasonality
# We create a matrix of 4 rows, calculate mean for each quarter
m_beer = t(matrix(data = detrend_beer, nrow = 4))
seasonal_beer = colMeans(m_beer, na.rm = T) 
plot(as.ts(rep(seasonal_beer,16)))

# Fourth, calculate random noise
random_beer = ts_beer - trend_beer - seasonal_beer
plot(as.ts(random_beer))

# Of course we can reconstruct the original signal by adding components
recomposed_beer = trend_beer+seasonal_beer+random_beer
plot(as.ts(recomposed_beer))

#------------------------------------------------------------------
# Additive decomposition using the DECOMPOSE() function:

decompose_beer = decompose(ts_beer, "additive")

plot(as.ts(decompose_beer$seasonal))
plot(as.ts(decompose_beer$trend))
plot(as.ts(decompose_beer$random))
plot(decompose_beer)

#------------------------------------------------------------------
# Multiplicative decomposition - by hand!!!

data(AirPassengers)
AP <- AirPassengers
plot(as.ts(AP))

# First, similarly to additive decomposition, we calculate annual moving average
trend_air <- ma(AP, order = 12, centre = T)  
plot(as.ts(AP))
lines(trend_air, col="red",lwd=2)

# Second, we detrend the time series by dividing 
detrend_air <- AP/trend_air
plot(as.ts(detrend_air))  # so this is again leftover seasonality and error

# Third, calculate average seasonality
m_air = t(matrix(data = detrend_air, nrow = 12))
seasonal_air = colMeans(m_air, na.rm = T)  # calculate mean for each quarter
plot(as.ts(rep(seasonal_air,12)))

# Fourth, calculate random noise
random_air = AP/(trend_air*seasonal_air)
plot(as.ts(random_air))

# Again, we can reconstruct the original signal by multplying components
recomposed_air = trend_air*seasonal_air*random_air
plot(as.ts(recomposed_air))

#------------------------------------------------------------------
# Multiplicative decomposition using the DECOMPOSE( ) function:

decompose_air <- decompose(AP, "multiplicative")

plot(as.ts(decompose_air$seasonal))
plot(as.ts(decompose_air$trend))
plot(as.ts(decompose_air$random))
plot(decompose_air)

#------------------------------------------------------------------
#------------------------------------------------------------------
# Additional options for decomposition

# There are many other options of decomposing time series. Often they differ with respect to filtering. 
# The function stl() performs a seasonal decomposition of a given time series by determining 
# the trend using 'loess' regression and then calculating the seasonal component 
# (and the residuals) from the differences

fit <- stl(AP, s.window=12)
plot(fit)


###################################################################
#
#                       SMOOTHERS/FILTERS
#
###################################################################


# ==============
# Exponential smoothing - models level component

# Exponential smoothing is a special case of the Holt-Winters algorithm with the additional parameters omitted.
# If we do not specify a value for alpha, R will find the value that minimizes the one-step-ahead prediction error.
# Coefficients - estimated value of the mean

# load data: monthly.csv
mydata <- read.csv(file.choose(),header=TRUE)
mydata <- mydata[1:840,]
ts_mydata <- ts(mydata[8], frequency = 12)

exp_smooth <- HoltWinters(na.omit(ts_mydata), beta = FALSE, gamma = FALSE)
# this is not same as exp_smooth <- HoltWinters(na.omit(ts_mydata), beta = 0, gamma = 0)
# this is forcing parameters to be zero
plot(exp_smooth)
exp_smooth$SSE
exp_smooth

exp_smooth <- HoltWinters(na.omit(ts_mydata), alpha=0.7, beta = FALSE, gamma = FALSE)
plot(exp_smooth)
exp_smooth$SSE
exp_smooth

predict <- predict(exp_smooth, n.ahead = 30) # 30 is too long
ts.plot(na.omit(ts_mydata), predict, lty = 1:2, col = c("black","red"))

# ==============
# Holt-Winters - level and trend components

exp_smooth <- HoltWinters(na.omit(ts_mydata), gamma=FALSE)
plot(exp_smooth)
exp_smooth$SSE

predict <- predict(exp_smooth, n.ahead = 30)
ts.plot(na.omit(ts_mydata), predict, lty = 1:2, col = c("black","red"))

# ==============
# Holt-Winters - level, trend, and seasonal components

# Here we are letting R choose the alpha, beta and gamma but in general you should set these values yourself!!!
# Typically values of beta and gamma are set about 0.2
# The default in R is to use values obtained from the decompose procedure

# =====
# additive

# Load data beer.csv
ts_beer<-read.csv(file.choose(),header=TRUE)
ts_beer <- ts(ts_beer, frequency = 4, start=c(1956,1))
plot(ts_beer)

AP.hw <- HoltWinters(ts_beer, seasonal = "additive")
plot(AP.hw)

AP.predict <- predict(AP.hw, n.ahead = 4 * 4) # predict for 4 years
ts.plot(ts_beer, AP.predict, lty = 1:2, col = c("black","red"))

# =====
# multiplicative

data(AirPassengers)
AP = AirPassengers
plot(as.ts(AP))

AP.hw <- HoltWinters(AP, seasonal = "multiplicative")
plot(AP.hw)

AP.predict <- predict(AP.hw, n.ahead = 4 * 12) # predict for 4 years
ts.plot(AP, AP.predict, lty = 1:2, col = c("black","red"))

# ==============
# Comparison of smoothers

ma_smooth <- filter(na.omit(ts_mydata), rep(1/4, 4), sides=1) # backward looking MA
exp_smooth <- HoltWinters(na.omit(ts_mydata), beta = FALSE, gamma = FALSE)
double_exp_smooth <- HoltWinters(na.omit(ts_mydata),gamma = FALSE)
hw_smooth <- HoltWinters(na.omit(ts_mydata))


plot(na.omit(ts_mydata)[time(na.omit(ts_mydata))>= 70],type="l",lwd=2,
     main="Comparison of smoothers",xlab="Time (months)",ylab="y") 
grid()
legend("topleft", inset=.02,
       c("Data", "Moving Average","Exponential","Double exponential", "Holt Winters"), 
       fill=c("black","gold3", "red", "chartreuse4", "darkmagenta"), horiz=F)

lines(ma_smooth[time(ma_smooth) >= 70], col = "gold3",lwd=2, lty = 6)
lines(exp_smooth$fitted[,1][time(exp_smooth$fitted) >= 70], col = "red",lwd=2, lty = 6)
lines(double_exp_smooth$fitted[,1][time(double_exp_smooth$fitted) >= 70], col = "chartreuse4", lwd=2, lty = 6)
lines(hw_smooth$fitted[,1][time(hw_smooth$fitted) >= 70], col = "darkmagenta",lwd=2, lty = 6)


###################################################################
#
#        BUILDING BLOCKS OF STOCHASTIC MODELS
#
###################################################################


# =====
# Gaussian white noise series 
set.seed(1)
w <- rnorm(200)
plot(w, type = "l")

# histogram
hist(w)


# =====
# simulate random walk
x <- w <- rnorm(200)
for (t in 2:200) x[t] <- x[t - 1] + w[t]
plot(x, type = "l")

acf(x)
acf(diff(x)) # differenced random walk

plot(diff(x), type = "l") # difference stationary


# =====
# simulate random walk with drift
x <- w <- rnorm(200)
d <- 0.15
for (t in 2:200) x[t] <- x[t - 1] + d + w[t]
plot(x, type = "l")



# =====
# simulate deterministic trend
x <- w <- rnorm(200)
dt <-seq(1:200) # time
beta <-0.015
for (t in 2:200) x[t] <- beta*dt[t] + w[t]
plot(x, type = "l")

acf(x, lag.max = 100)

time <- 1:length(x)
plot(resid(lm(x~time)), type = "l") # trend stationary


###################################################################
#
#                          ARMA
#
###################################################################


# ======
# Simulate AR(1) process
# this is one possible realisation of the model!

set.seed(12345)
y <- w <- rnorm(1000)
for (t in 2:1000) y[t] <- 0.7 * y[t - 1] + w[t]
layout(1:1)
plot(y, type = "l")

# ACF starts at lag 0
acf(y) 
# the acf will not tell us anything about autoregressive order

# PACF starts at lag 1
pacf(y) 
# the pacf indicates autoregressive order. Here we have significant correlation at lag 1

# ======
# autoregressive model fitted to the simulated series with approximately  95% confidence interval of parameter
# ar function by default selects the complexity (order) by AIC
# the estimation can be based on maximum likelihood estimation, yule-walker or OLS

y.ar <- ar(y, method = "mle")

y.ar$order
y.ar$ar                                       # estimate of autoregressive parameter 
y.ar$ar + c(-2, 2) * sqrt(y.ar$asy.var.coef)  # variance is extracted using y.ar$asy.var.coef by AIC
y.ar$ar + c(-2, 2) * as.vector(sqrt(y.ar$asy.var.coef))  # without Warning message...
mean(y)                                       # mean should be 0, but due to chance variation it will never be exactly 0

# ======
# Let's compare various estimates of AR(1)

# OLS estimate
z<-embed(y,2) # creating dataset with lagged dependent variable

summary(lm(z[,1]~z[,2])) # simple regression

# default, Yule-Walker estimate
y.ar <- ar(y)                  
y.ar$ar 
sqrt(y.ar$asy.var)

# maximum likelihood estimate
y.ar <- ar(y, method = "mle") 
y.ar$ar 
sqrt(y.ar$asy.var)

# ======
# predict using AR model

x <- predict(y.ar,  n.ahead = 10, se.fit = TRUE)
z <- c(y, x$pred[1:10])

plot(z[980:1010], type = "l", ylim = c(-5,5))
abline(h=mean(y), lty = 2, col='darkgray')  # plot mean 
abline(v=21,  col = "red",lwd=2, lty = 6)   # indicate the last in-sample point 

# ======
# Simulate AR(2) process

set.seed(1234)
y <- w <- rnorm(1000)
for (t in 3:1000) y[t] <- 0.4 * y[t - 1] +  0.2 * y[t - 2] + w[t]

acf(y)
pacf(y)

#------------------------------------------------------------------
# Moving average models
#------------------------------------------------------------------

# An MA(q) model can be fitted to data in R using the arima function with the order function parameter set to c(0,0,q). 
# Unlike the function ar, the function arima does not subtract the mean by default and estimates an intercept term. 


# ======
# Code used to simulate the MA(1)
set.seed(12345)
x <- w <- rnorm(1000)
for (t in 2:1000) x[t] <- x[t] + 0.6 * w[t - 1]
plot(x, type = "l")

acf(x)  # remember, ACF starts at lag 0
pacf(x) # PACF starts at lag 1

arima(x, order = c(0, 0, 1))


# ======
# Code used to simulate the MA(3)

set.seed(1)
b <- c(0.8, 0.6, 0.4)
x <- w <- rnorm(1000)
for (t in 4:1000) {
  for (j in 1:3) x[t] <- x[t] + b[j] * w[t - j]
}
plot(x, type = "l")

acf(x)
pacf(x)

# Moving average model is fitted to the simulated series
# For this we have to use function arima

arima(x, order = c(0, 0, 3))

