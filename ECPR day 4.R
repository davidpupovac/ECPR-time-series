###################################################################
###################################################################
#
#             Times Series Analysis - Day 4
#
###################################################################
###################################################################

setwd("D:/My Documents/Desktop/ECPR data") 

###################################################################
#
#                Seasonal components
#
###################################################################

#------------------------------------------------------------------
# Linear seasonal model
#------------------------------------------------------------------

# Seasonal model for the temperature series
# load data: global.csv
Global.ts<-read.csv(file.choose(),header=TRUE)

# whole serie
Global.ts <- ts(Global.ts, st = c(1856, 1), end = c(2005,12), fr = 12)
temp <- window(Global.ts, start = 1970)
ts.plot(temp)

Seas <- cycle(temp) # seasonality 
Time <- time(temp)  # time

# -----
temp.lm <- lm(temp ~ Time)
coef(temp.lm)

ts.plot(predict(temp.lm), temp, col=c("red", "black")) # line indicates trend

# If it the equation is estimated by OLS, a zero is used within the formula to ensure 
# that the model does not have an intercept. If the intercept is included in the formula,
# one of the seasonal terms will be dropped and an estimate for the intercept will appear 
# in the output. With or without an intercept, would be equivalent, as can be easily 
# verified by rerunning the algorithm without the zero in the formula. 

# OLS
temp.lm <- lm(temp ~ 0 + Time + factor(Seas))
coef(temp.lm)

# seasonal variation of global temperatures
ts.plot(coef(temp.lm)[2:13])

# prediction from the model
ts.plot(predict(temp.lm))

# lets compare to the original series
ts.plot(predict(temp.lm), temp, col=c("red", "black"))

# ======

# forecast based on OLS
new.t <- seq(2006, len = 2 * 12, by = 1/12)
trend <- coef(temp.lm)[1]
seasonal <- rep(coef(temp.lm)[2:13], 2)

ts.plot((trend * new.t + seasonal)[1:24]) # forecast for two years

# alternatively you can use predict

new.dat <- data.frame(Time = new.t, Seas = rep(1:12, 2))
predict(temp.lm, new.dat)[1:24]
ts.plot(predict(temp.lm, new.dat)[1:24])

#------------------------------------------------------------------
# Seasonal ARIMA models
#------------------------------------------------------------------


# ======
# Additive SARIMA
# ======


# load data wpi.csv
# U.S. Wholesale Price Index

# This is replication of analysis from Enders (2004) - Applied Econometric Time Series model.

wpi <-read.csv(file.choose(),header=TRUE)
wpi <- ts(wpi, start = 1960, freq = 4)
plot(wpi)

plot(diff(wpi))
plot(diff(log(wpi)))

diff_wpi <- diff(log(wpi))

# Enders identified appropriate model as ARMA(1,1) . However, an MA(4) term is also included to 
# account for a remaining quarterly effect. This is our starting point. First we will estimate full ARMA(1,4)

arima(diff_wpi, order = c(1, 0, 4)) # MA(2) and MA(3) coefficients are insignificant
           
# to estimate an ARIMA model with AR(1), MA(1) and MA(4) components
# you have to put constrains on parameters. NA indicates free elements, 
# while O (zero) indicates constrained elements. 
# With command fixed=c(NA,NA,0,0,NA,NA) we are specifying: 
# AR(1) varies, MA(1) varies, MA(2)=0, MA(3)=0, MA(4) varies, intercept varies

arima(diff_wpi, order = c(1, 0, 4), 
      fixed=c(NA,NA,0,0,NA,NA))

# let's compare models using AIC:
AIC (arima(diff_wpi, order = c(1, 0, 4), 
           fixed=c(NA,NA,0,0,NA,NA)))

AIC (arima(diff_wpi, order = c(1, 0, 4), # here we also allow MA(2) to vary
           fixed=c(NA,NA,NA,0,NA,NA)))

AIC (arima(diff_wpi, order = c(1, 0, 4))) # here we allow everyting to varies

# ======
# check the residuals
res_diff_wpi<- resid(arima(diff_wpi, order = c(1, 0, 4), 
      fixed=c(NA,NA,0,0,NA,NA)))

# check autocorrelation function
library(forecast)
Acf(res_diff_wpi)

# The null hypothesis of uncorrelatedness 
Box.test(res_diff_wpi, lag=20, fitdf=4, type="Ljung")

# The hypothesis of normally distributed errors
shapiro.test(res_diff_wpi)

library(tseries)
jarque.bera.test(res_diff_wpi)

# The null normality hypothesis 
qqnorm(res_diff_wpi)

hist(res_diff_wpi)


# ======
# Multiplicative SARIMA
# ======


data(AirPassengers)
AP <- AirPassengers

plot(as.ts(AP))

# Frist, let's take a log of series to enforce constant variance
plot(log(AP))

# Second, I difference to induce stationarity 
# (you would probably use unit roots test to check if this is adequate)
plot(diff(log(AP)))

diff_AP <- diff(log(AP))

Acf(diff_AP) 
Pacf(diff_AP) # this does not tell us much, increase lag length

Acf(diff_AP, lag=60)
Pacf(diff_AP, lag=60)
# obviously we have here a non-stationary process in seasonal component

diff_diff_AP <-diff(diff_AP,12) # seasonal differencing

Acf(diff_diff_AP, lag=60)
Pacf(diff_diff_AP, lag=60)
# we still have some significant correlations at seasonal component,
# but maybe we get rid of it with modeling

# ======
# comparison of  models

AIC (arima(diff_diff_AP, order = c(1,0,0),
           seas = list(order = c(1,0,0), 12)))
AIC (arima(diff_diff_AP, order = c(1,0,0),
           seas = list(order = c(0,0,1), 12)))
AIC (arima(diff_diff_AP, order = c(0,0,1),
           seas = list(order = c(1,0,0), 12)))
AIC (arima(diff_diff_AP, order = c(0,0,1),
           seas = list(order = c(0,0,1), 12)))
AIC (arima(diff_diff_AP, order = c(1,0,1),
           seas = list(order = c(0,0,1), 12)))
AIC (arima(diff_diff_AP, order = c(0,0,2),
           seas = list(order = c(0,0,1), 12)))

# ======
# final model

arima(diff_diff_AP, order = c(0,0,1),
      seas = list(order = c(0,0,1), 12))

# =====
res <- resid(arima(diff_diff_AP, order = c(0,0,1),
                   seas = list(order = c(0,0,1), 12)))

# check autocorrelation function
acf(res)

Box.test(res, lag=20, fitdf=4, type="Ljung")
# The null hypothesis of uncorrelatedness 

shapiro.test(res)
# the null normality hypothesis 

qqnorm(res)

hist(res)

# ===============================
# there is automated version as well

auto.arima(log(AP))


# ===============================
# we can forecast from this model

# notice that we do not use function 'arima' but function 'Arima' from forecast package
sarima011_011 <- Arima(log(AP), order = c(0,1,1),
             seas = list(order = c(0,1,1), 12))

plot(forecast(sarima011_011, h=12))
plot(forecast(sarima011_011, h=12), xlim=c(1958,1962), ylim=c(5.5,7))

#------------------------------------------------------------------
#   Seasonal Unit Roots
#------------------------------------------------------------------

library(urca)

# you can also use library(uroot)
# install.packages("https://cran.r-project.org/src/contrib/Archive/uroot/uroot_1.4.tar.gz", repos = NULL, type="source")

library(pdR) 

# this is replication analysis of Hylleberg et al. [1990].

data(UKconinc) # real disposable income in the United Kingdom
incl<-ts(UKconinc$incl,start=c(1955,1),
         end=c(1984,4),frequency=4)

plot(incl)

# The authors have chosen the lags 1, 4, 5 in the augmented test regression 
# These are differenced lags as in augmented Dickey Fuller test

# Here we are testing if there are seasonal unit root if  deterministic terms (intercept,  trend, 
# seasonal dummy variables) or as well as lagged seasonal differences, are added to the test regression

HEGY.test(wts=incl,itsd=c(0,0,c(0)),      # no intercept, no linear trend in equation
          selectlags=list(mode=c(1,4,5)))

HEGY.test(wts=incl,itsd=c(1,0,c(0)),      # intercept, no linear trend in equation
          selectlags=list(mode=c(1,4,5)))

HEGY.test(wts=incl,itsd=c(1,1,c(0)),      # intercept, linear trend in equation
          selectlags=list(mode=c(1,4,5)))

HEGY.test(wts=incl,itsd=c(1,0,c(1,2,3)),  # intercept, no linear trend and 3 dummies in equation
          selectlags=list(mode=c(1,4,5)))

HEGY.test(wts=incl,itsd=c(1,1,c(1,2,3)),  # intercept, no linear trend and 3 dummies in equation
          selectlags=list(mode=c(1,4,5)))

# tpi_1      if pi significant =  no usual (non-seasonal) stochastic stationary component
# tpi_1      if pi significant =  no evidence of a biannual (twice a year) cycle unit root
# Fpi_3:4    if pi significant (on basis of F test) =  no evidence of annual cycle unit root

# -----------
# The specification of the test regression is determined by itsd, a three-element vector:
# If the first element is set to one, a constant is included, 
# If the second element of itsd is set to one, a linear trend is included, 
# The inclusion of seasonal dummy variables is controlled by the third element

# The inclusion of lagged seasonal differences is set by the argument selectlags
# Additional regressors can be included by the argument regvar (allows to model 
# structural breaks in the seasonal means)

###################################################################
#
#                      ARCH/GARCH models
#
###################################################################


library(MASS)
data(SP500)
plot(SP500, type = 'l')

library(forecast)
Acf(SP500, lag=60)
# Note that the correlogram of a volatile series does not differ significantly from white noise
# but the series is non-stationary since the variance is different at different times.

# If a correlogram appears to be white noise then volatility can be detected by looking at 
# the correlogram of the squared values since the squared values are equivalent to 
# the variance (provided the series is adjusted to have a mean of zero).

Acf((SP500 - mean(SP500))^2, lag=60)
Pacf((SP500 - mean(SP500))^2, lag=60)

# there is evidence of serial correlation in the squared values, 
# so there is evidence of conditional heteroskedastic behaviour and volatility

# ======
# ARCH

# =====
# simulate arch(1) process

n=1000
apha1=0.5
apha0=0.05

set.seed(12345)
w=rnorm(n)
epsilon=rnorm(n)
sigma2=w

for(t in 2:n){
  sigma2[t]=apha0+apha1*epsilon[t-1]^2
  epsilon[t]=w[t]*sqrt(sigma2[t])
}

plot(epsilon,type="l")

# simulate AR(1) with ARCH errors
y <- rnorm(1000)
for (t in 2:1000) y[t] <- 0.7*y[t - 1]  + epsilon[t]

ts.plot(y)

# ======
Acf(y)
Pacf(y)

Acf(epsilon, lag.max = 100) # this seems like a white noise proces

Acf(epsilon^2)
Pacf(epsilon^2)

# ======
arima(y,order=c(1,0,0))

auto.arima(y) # just a comparision with auto arima function, lesson: do not use abuse automated version

# ======
# check residuals

res_arma10 <- resid(arima(y,order=c(1,0,0)))

plot(res_arma10)


# ======
# determine the lag order

install.packages("https://cran.r-project.org/src/contrib/Archive/FinTS/FinTS_0.4-5.tar.gz", repos = NULL, type="source")
library(FinTS) 

# Lagrange Multiplier (LM) test for ARCH
ArchTest(res_arma10,lags=3) # at least one component is significant
ArchTest(res_arma10,lags=2) # at least one component is significant
ArchTest(res_arma10,lags=1) 

Acf(res_arma10^2)
Pacf(res_arma10^2)

# ======
# estimate ARCH 

library(tseries)

arch_1 <- garch(res_arma10,order=c(0,1),trace=F) # to set ARCH model order=c(0,p)
summary(arch_1)

arch_2 <- garch(res_arma10,order=c(0,2),trace=F) # overfit to check if there is need for additional term
summary(arch_2)

# ======
# additional residual analysis of ARCH

res_arch_1<- resid(garch(res_arma10,order=c(0,1),trace=F))

shapiro.test(res_arch_1)
# the normality hypothesis cannot be rejected

qqnorm(res_arch_1)
# quantiles plot seems straight

#------------------------------------------------------------------

# ======
# Simulation of GARCH errors

set.seed(1)
alpha0 <- 0.1
alpha1 <- 0.4
beta1 <- 0.2
w <- rnorm(10000)
h <- rep(0, 10000)
epsilon <- rep(0, 10000)
for (i in 2:10000) {
  h[i] <- alpha0 + alpha1 * (epsilon[i - 1]^2) + beta1 * h[i - 1]
  epsilon[i] <- w[i] * sqrt(h[i])
}

ts.plot(epsilon)

Acf(epsilon)
Acf(epsilon^2)
Pacf(epsilon^2)

# ======
# fitting garch model

library(tseries)

# The default is GARCH(1,1), higher-order models can be specified with the parameter order=c(p,q) 
epsilon.garch <- garch(epsilon, grad = "numerical", trace = FALSE)

summary(epsilon.garch)
confint(epsilon.garch)

###################################################################
#
#                   ARMAX/ARIMAX models
#
###################################################################

# =====
# ARMAX

# load data icecream.csv
# sales of icecream during a month

Icecream <- read.csv(file.choose(),header=TRUE)
Icecream<-data.frame(Icecream)

plot.ts(Icecream[,2]) 

auto.arima(Icecream$cons)

arima(Icecream$cons,order=c(2,0,1), xreg=Icecream[3:5])
arima(Icecream$cons,order=c(1,0,1), xreg=Icecream[3:5]) 
arima(Icecream$cons,order=c(1,0,0), xreg=Icecream[3:5]) 


# =====
# ARIMAX

# be careful with this, it may lead to seriously wrong inference!!!  
# we will address these problems in the next class

set.seed(12345)

y <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = 200)
x <- arima.sim(list(order = c(0,1,0)), n = 200)

plot(y, x, pch=1, cex=0.5)
cor(y, x)    

ts.plot(y)
ts.plot(x)

ts.plot(diff(y))
ts.plot(diff(x))

cor(diff(y),diff(x))

ar(x = diff(y))

arima(diff(y),order=c(1,0,0), xreg=diff(x)) 

# =====
# If you are interested in transfer function models you can experiment with package TSA

library(TSA)
arimax(y, order = c(1, 1, 0), xreg=x)

