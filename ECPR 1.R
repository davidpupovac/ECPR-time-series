###################################################################
###################################################################
#
#                 Times Series Analysis - Day 1
# 
###################################################################
###################################################################

# let us install the packages we will use during the course

doInstall <- TRUE  
toInstall<-c("lmtest","sandwich","car","nlme","orcutt",
             "forecast","Ecdat","fArma","urca","topicmodels",
             "tseries","uroot","MASS","FinTS","fGarch","mvtnorm",
             "devtools","vars","dplyr","urca","quantmod", "pdR")
if(doInstall){install.packages(toInstall, repos = "http://cran.r-project.org")}


# package ‘FinTS’ is not available (for R version 3.4.1)

###################################################################
#
#                          BASICS
#
###################################################################


#  How to declare time series objects?
ts(1:100, frequency = 4, start = c(1959, 2)) # quarterly
ts(1:100, start=c(2009, 1), end=c(2014, 12), frequency=12) # monthly 
ts(1:100, frequency = 1, start = c(1993)) # annually

# How to create multiple time series?
working_data <- ts(matrix(rnorm(300), 100, 3), start = c(1961, 1), frequency = 12)

# Let's check how this data set looks and works 
working_data

class(working_data)
methods(class = "ts")

# Most of the typical subsetting functions work with ts data
dim(working_data)
colnames(working_data)
head(working_data)
working_data[3:6,2:3]
head(working_data[,3])

# However, there are ts subsetting options
ob<-window(working_data, start=c(1965,1), end=c(1967,12))

#------------------------------------------------------------------

# Options for plotting time series data

# Plotting multiple time series on separate panels
plot(working_data)
plot.ts(working_data)

# Plotting multiple time series on a single panel
ts.plot(working_data)

# Additional options for plotting
plot(working_data, plot.type = "single", lty = 4:6, lwd=3, col=c("gold3", "red", "chartreuse4"))
plot(working_data, plot.type = "multiple", lwd=1, col=c("red"))

#------------------------------------------------------------------

# Examples with real time series data 

data(AirPassengers)
AP <- AirPassengers
AP

plot(AP, ylab = "Passengers (1000's)")

# Extract time from a time series object
time(AP) 

# Return the cycles for each value in a series
cycle(AP) 

# Extract star and end dates, and the frequency  
start(AP); end(AP); frequency(AP)

layout(1:3)
plot(AP)   # note "aggregate" is not base aggregate function requiring "by"
plot(aggregate(AP, FUN=mean)) 
plot(aggregate(AP, FUN=sd))  

# Some usual graphics functions can be applied to the data
layout(1:1)
boxplot(AP ~ cycle(AP))

#------------------------------------------------------------------

# Additional functions for time series data 

setwd("D:/My Documents/Desktop/ECPR data") 

# Unemployment in the United States January 1996-October 2006.
# Load data monthly.csv
ts_mydata<-read.csv(file.choose(),header=TRUE)
ts_mydata <- ts(ts_mydata, frequency = 12, start=c(1913,1))

# Let us look at this data. Here we have some lack of overlap in data. 
plot(ts_mydata)

# Let's make lack of overlap more obvious for two variables. First select the window for which we have data.
unemp.ts <- window(ts_mydata[, 8], start = c(1950,1), end=c(2000,12) ,freq = 12) # unemployment 
pmmms30.ts <- window(ts_mydata[, 7], start = c(1980,1), freq = 12) #  mortgage rate

plot(cbind(unemp.ts, pmmms30.ts))

# Let us try to combine these data series without missing  
Overlap <- ts.intersect(unemp.ts, pmmms30.ts)
plot(Overlap)

start(Overlap)
end(Overlap)

#------------------------------------------------------------------
# Lagging variables 

set.seed(1234)
y <- rnorm(10)

# this is the best way to lag variables (IMHO)
# install.packages("dplyr")
library("dplyr")
lag_y <- dplyr::lag(y,1) 
lead_y <- dplyr::lead(y,1) 

# other option 1
stats::lag(y,-1) 
stats::lag(y,1) 

# other option 2
embed(y,dimension=2) # you essentially lose a row here
head(y)

###################################################################
#
#                          Autocorrelation
#
###################################################################

#------------------------------------------------------------------
# Generalised least squares

# ======
# Simulate data

set.seed(12345)
e <- u <- rnorm(100, sd = 2.5)
for (t in 2:100) e[t] <- 0.8 * e[t - 1] + u[t]

x1 <- rnorm(100, mean=15, sd = 2)
x2 <- rnorm(100, mean=30, sd = 3)

y <- 20 + 4 * x1 + 6 * x2 + e
gls_data<- data.frame(y, x1, x2)

plot(y, xlab = "time", type = "l")
plot.ts(gls_data)

# ------ 
# Let's fit OLS

lm_fit<- lm(y ~ x1 + x2, data = gls_data) 
summary(lm_fit)

# ======
# Generalized least squares

library(nlme)

# lag 1 autocorrelation of 0.8 is used above because this value was used to simulate the data
gls_fit <- gls(y ~ x1 + x2, cor = corAR1(0.8))

# ols estimate
summary(lm_fit)

# compare to gls coefficients - quite different
coef(gls_fit)

# standard errors of the parameters are considerably greater than those obtained from lm 
sqrt(diag(vcov(gls_fit)))

# confidence intervals
confint(gls_fit)

# ======
# Feasible generalized least squares

# There is a package that estimates FGLS 
# install.packages("lmtest")
library(lmtest)
library(orcutt)

lm_fit <- lm(y ~ x1 + x2, data=gls_data)
fgls_fit <- cochrane.orcutt(lm_fit)

fgls_fit

summary.orcutt(fgls_fit)
fgls_fit$std.error
summary(lm_fit)$coefficients[,2]

#------------------------------------------------------------------
# Feasible generalized least squares

# from John Fox - McMaster University 

cochrane.orcutt.lm <- function(mod){
  X <- model.matrix(mod)
  y <- model.response(model.frame(mod))
  e <- residuals(mod)
  n <- length(e)
  names <- colnames(X)
  rho <- sum(e[1:(n-1)]*e[2:n])/sum(e^2)
  y <- y[2:n] - rho * y[1:(n-1)]
  X <- X[2:n,] - rho * X[1:(n-1),]
  mod <- lm(y ~ X - 1)
  result <- list()
  result$coefficients <- coef(mod)
  names(result$coefficients) <- names
  summary <- summary(mod, corr = F)
  result$cov <- (summary$sigma^2) * summary$cov.unscaled
  dimnames(result$cov) <- list(names, names)
  result$sigma <- summary$sigma
  result$rho <- rho
  class(result) <- 'cochrane.orcutt'
  result
}

# coefficients
cochrane.orcutt.lm(lm_fit)$coefficients

# standard errors
sqrt(diag(cochrane.orcutt.lm(lm_fit)$cov))

# ======

prais.winsten.lm <- function(mod){
  X <- model.matrix(mod)
  y <- model.response(model.frame(mod))
  e <- residuals(mod)
  n <- length(e)
  names <- colnames(X)
  rho <- sum(e[1:(n-1)]*e[2:n])/sum(e^2)
  y <- c(y[1] * (1 - rho^2)^0.5, y[2:n] - rho * y[1:(n-1)])
  X <- rbind(X[1,] * (1 - rho^2)^0.5, X[2:n,] - rho * X[1:(n-1),])
  mod <- lm(y ~ X - 1)
  result <- list()
  result$coefficients <- coef(mod)
  names(result$coefficients) <- names
  summary <- summary(mod, corr = F)
  result$cov <- (summary$sigma^2) * summary$cov.unscaled
  dimnames(result$cov) <- list(names, names)
  result$sigma <- summary$sigma
  result$rho <- rho
  class(result) <- 'prais.winsten'
  result
}

# coefficients
prais.winsten.lm(lm_fit)$coefficients

# standard errors
sqrt(diag(prais.winsten.lm(lm_fit)$cov))

#------------------------------------------------------------------
#------------------------------------------------------------------
# let's segue a bit...

# Heteroscedastic data

# ======
# Define some heteroscedastic data
set.seed(1234)

# Generate uniformly distributed variable (no negative values like in rnorm)
x <- runif(1000, 0, 1)  
# For each successive iteration of random number generation we have a larger standard deviation (x)
y <- rnorm(1000, 0, x)  

#------------------------------------------------------------------
# Detecting heteroscedastic data

# ======
# Plot the data
plot(x,y)

het_fit <- lm(y~x)
plot(het_fit)

# ======
# Non-robust standard error
summary(het_fit)

# ======
# Breusch-Pagan test of homoscedasticity
library(lmtest)
bptest(het_fit)  # (H0) there is no heteroscedastic among residuals

#------------------------------------------------------------------
# Getting covariance matrices

# ======
library(sandwich)

vcov(het_fit) 

vcovHC(het_fit, "HC0") # White (1980)
vcovHC(het_fit, "HC1") # MacKinnon and White (1985) - for small samples

# ======
# Getting corrected standard errors


summary(het_fit)$coefficients[,2] # original estimate

# library(lmtest)

coeftest(het_fit, vcov = NULL)  # original estimate with package lmtest

coeftest(het_fit, vcov = vcovHC(het_fit, "HC0"))    # robust; HC0 
coeftest(het_fit, vcov = vcovHC(het_fit, "HC1"))    # robust; HC1 (Stata default)

#------------------------------------------------------------------
#------------------------------------------------------------------

# Autocorrelated data

# ======
# Simulate process with autocorrelated errors
set.seed(1234)
e <- u <- rnorm(100)

# create autocorrelated errors; rho=0.3
for (t in 2:100) e[t] <- 0.60 * e[t - 1] + u[t] 
hist(e)

# create independent variable
x <- rnorm(100,mean=50,sd=2)

# create dependent variable
y = 30 + 3*x + e
ac_errors <- ts(cbind(y,x))

# estimate the model
fit_ac<- lm(y~x, data=ac_errors)

# save residuals
resid <- residuals(fit_ac)

# compare errors
cor(e,resid)

#------------------------------------------------------------------
# Detection of autocorrelation

# =====
# detection via visualisation 
ts.plot(resid)

# =====
# Durbin-Watson test 
# library(lmtest)
dwtest(fit_ac) # (H0) is there is no autocorrelation among residuals

library(car)
durbinWatsonTest(fit_ac) # here we are also getting estiamte of rho and lag

# =====
# Portmanteau - Ljung-Box and Box-Pierce tests

# (H0) is that there is no correlation among residuals
Box.test(resid, lag = 1, type = "Ljung-Box", fitdf = 0) 
Box.test(resid, lag = 1, type = "Box-Pierce", fitdf = 0)

# =====
# Breusch-Godfrey 

# library(lmtest)
bgtest(y~x,order=1,type="Chisq") # (H0) is that there is no correlation among residuals
bgtest(y~x,order=1,type="F")

# auxiliary regression that includes the lagged residuals

coeftest(bgtest(fit_ac,order=1,type="Chisq"))

#------------------------------------------------------------------
# Getting Newey-West standard errors

library(sandwich)
library(lmtest)

fit_ac<- lm(y~x, data=ac_errors)
summary(fit_ac)

coeftest(fit_ac, vcov = NeweyWest) # Newey-West standard errors

#-------------------------
# ??????????????????????
#-------------------------

# exogenity

lag_y <- dplyr::lag(y,1) 

data <- data.frame(ac_errors,cbind(lag_y,resid))
cor(data$lag_y,data$resid, use="complete.obs")

# ---
fit_ac<- lm(y~x+lag_y, data=data)
resid <- residuals(fit_ac)

# (H0) is that there is no correlation among residuals
Box.test(resid, lag = 1, type = "Ljung-Box", fitdf = 0) 
Box.test(resid, lag = 1, type = "Box-Pierce", fitdf = 0)

cor(cbind(na.omit(data$lag_y),resid), use="complete.obs")
