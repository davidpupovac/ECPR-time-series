###################################################################
###################################################################
#
#                 Times Series Analysis - Day 3 
#
###################################################################
###################################################################


setwd("D:/My Documents/Desktop/ECPR data") 


###################################################################
#
# Example of unit roots within the circle
# 
###################################################################

# -----------------------------------------------------------------
# demostration of stationarity AR  

set.seed(1234)
y <- w <- rnorm(100)
for (t in 3:100) y[t] <- -0.4 * y[t - 1] -  0.8 * y[t - 2] + w[t]

ts.plot(y)

arma<-arima(y,c(2,0,0),include.mean=FALSE, # include.mean=FALSE - we do not want intercept
            transform.pars=FALSE,method="ML")
arma

# checking if the estimated values satisfy the stability condition.
x<-seq(-1,1,length=100)
y1<-sqrt(1-x^2)
y2<--sqrt(1-x^2)
plot(c(x,x),c(y1,y2),xlab='Real part', col="red",
     ylab='Complex part',type='l',
     main='Unit Circle',ylim=c(-2,2),xlim=c(-2,2))
abline(h=0)
abline(v=0)

points(Re(polyroot(c(1,-arma$coef))), # this commands will be discussed later
       Im(polyroot(c(1,-arma$coef))),pch=19) 
legend(-1.5,-1.5,legend="Roots of AR(2)",pch=19)

# there is a function that will allows us to assess stability 
# of AR and MA components. For every component you plug relevant coefficients

library(fArma)
armaRoots(c(-0.4039, -0.7878)) # in this fuction you have to put estimated coefficients

# ------
# explosive processes 1 

y <- w <- rnorm(200)
for (t in 3:200) y[t] <- 0.5 * y[t - 1] + 0.9 * y[t - 2] + w[t]

ts.plot(y)

x<-seq(-1,1,length=200)
y1<-sqrt(1-x^2)
y2<--sqrt(1-x^2)
plot(c(x,x),c(y1,y2),xlab='Real part', col="red",
     ylab='Complex part',type='l',
     main='Unit Circle',ylim=c(-2,2),xlim=c(-2,2))
abline(h=0)
abline(v=0)

points(Re(polyroot(c(1,-0.5, -0.9))),
       Im(polyroot(c(1,-0.5, -0.9))),pch=19)
legend(-1.8,-1.8,legend="Roots of the process",pch=19)

armaRoots(c(0.5, 0.9))

# we could try to estimate this model using arima command, 
# but the extreme nature of series will create problems

# ------
# explosive processes 2

set.seed(1234)
y <- w <- rnorm(100)
for (t in 3:100) y[t] <- -0.4 * y[t - 1] +  0.7 * y[t - 2] + w[t]

ts.plot(y)

armaRoots(c(-0.4, 0.7))

# -------
# unit roots process

y <- armaSim(n=100, model=list(ar = c(-0.5, 0.9, -0.1, -0.5)))
# here we are using armaSim command instead of our usual simulation 

ts.plot(y)

coefficients = c(-0.5, 0.9, -0.1, -0.5)

# The moduli of the characteristic polynomial are retrieved with Mod() and the real 
# and complex parts with the functions Re() and Im(), respectively.

Re(polyroot(c(1,0.5, -0.9, 0.1, 0.5))) 
Im(polyroot(c(1,0.5, -0.9, 0.1, 0.5)))
# notice that coefficient must be negative
# notice 1 at the beginning 

Mod(polyroot(c(1,0.5, -0.9, 0.1, 0.5)))

armaRoots(coefficients)


###################################################################
#
#       Autocorrelations and partial autocorrelations
#
###################################################################

# ---------------------------------
# white noise

set.seed(12345)
w <- rnorm(100)

plot.ts(w,ylab='')
acf(w,main='Autocorrelations',ylab='',
    ylim=c(-1,1),ci.col="red")
pacf(w,main='Partial Autocorrelations',ylab='',
     ylim=c(-1,1),ci.col="red")

# ---------------------------------
# random walk

set.seed(12345)
y <- w <- rnorm(100)
for (t in 2:100) y[t] <- y[t - 1]  + w[t]

plot.ts(y,ylab='')
acf(y,main='Autocorrelations',ylab='',
    ylim=c(-1,1),ci.col="red")
pacf(y,main='Partial Autocorrelations',ylab='',
     ylim=c(-1,1),ci.col="red")

# ---------------------------------
# AR(1)

set.seed(12345)
y <- arima.sim(n=100,list(ar=0.7),innov=rnorm(100))

plot.ts(y,ylab='')
acf(y,main='Autocorrelations',ylab='',
    ylim=c(-1,1),ci.col="red")
pacf(y,main='Partial Autocorrelations',ylab='',
     ylim=c(-1,1),ci.col="red")

# ---------------------------------
# MA(1)

set.seed(123456)
y <- arima.sim(n=100,list(ma=0.7),innov=rnorm(100))

plot.ts(y,ylab='')
acf(y,main='Autocorrelations',ylab='',
    ylim=c(-1,1),ci.col="red")
pacf(y,main='Partial Autocorrelations',ylab='',
     ylim=c(-1,1),ci.col="red")

# ---------------------------------
# ARMA(2,1)

set.seed(12345)
y <- arima.sim(n=100,list(ar = c(0.4, 0.2), ma=0.1),innov=rnorm(100))

acf(y,main='Autocorrelations',ylab='',
    ylim=c(-1,1),ci.col="red")
pacf(y,main='Partial Autocorrelations',ylab='',
     ylim=c(-1,1),ci.col="red")


###################################################################
#
#               Box-Jenkins approach 
#
###################################################################

# -----------------------------------------------------------------
#                   Example 1
# -----------------------------------------------------------------

library(urca)

# unemployment rate in the United States
data(npext)
y<-ts(na.omit(npext$unemploy),start=1909,end=1988,
      frequency=1) 

# ======

plot(y,ylab="Unemployment rate") # no trending behavior

# examine autocorrelation functions
acf(y,main='Autocorrelations',ylab='',ylim=c(-1,1))
# if you want to lose 0 lag do the following:
acf(y,main='Autocorrelations',ylab='',ylim=c(-1,1), xlim=c(1,20))

pacf(y,main='Partial Autocorrelations',ylab='', ylim=c(-1,1))
# autocorrelation function tapers off
# partial autocorrelation function has two significant correlations

# ======
# tentative ARMA(2,0) 

arma20<-arima(y,order=c(2,0,0)) # function arima() contained in the package stats
arma20 # check significance (i.e. standard error * 2, standard error * -2)

# ------
# checking if the estimated values satisfy the stability condition (unit roots).

polyroot(c(1,-arma20$coef[1:2]))
Re(polyroot(c(1,-arma20$coef[1:2])))
Im(polyroot(c(1,-arma20$coef[1:2])))

Mod(polyroot(c(1,-arma20$coef[1:2])))

armaRoots(arma20$coef[1:2]) # or we can just plot it on unit circle

# ------
# information criteria and log likelihood

ll20<- logLik(arma20) # saving logLik
aic20<- arma20$aic    # saving AIC

# ------
# residuals are retrieved and stored
res20<-residuals(arma20)

# residuals can be inspected visually,as can their autocorrelation functions (ACF)
# and partial autocorrelation functions (PACF)
plot(res20,ylab="Residual")
acf(res20,main="Residual acf",ylab="",ylim=c(-1,1), xlim=c(1,20))
pacf(res20,main="Residual pacf",ylab="", ylim=c(-1,1))

# ------
# checking uncorrelatedness

Box.test(res20,lag=20,type="Ljung-Box") # Ljung-Box Portmanteau test
# The null hypothesis of uncorrelatedness up to order 20 cannot be rejected

# ------
# The hypothesis of normally distributed errors

shapiro.test(res20)
# the normality hypothesis cannot be rejected

# or you can use Jarque-Bera test
library(tseries)
jarque.bera.test(res20)
# the normality hypothesis cannot be rejected

# or you can use normal quantiles plot
qqnorm(res20)

# ====================
# tentative ARMA (3,0) - overfitting

arma30<-arima(y,order=c(3,0,0))
arma30
# coefficient for the third lag is not significantly different from zero
# estimates for the first- and second-order AR-coefficients are almost unchanged.

# ------
armaRoots(arma30$coef[1:3])

# ------
# information criteria and log likelihood

ll30 <- logLik(arma30)
aic30 <- arma30$aic

# ------
# checking residuals

res30<-residuals(arma30)

plot(res30,ylab="Residual")

Box.test(res30,lag=20,type="Ljung-Box") 
# The null hypothesis of uncorrelatedness up to order 20 cannot be rejected

shapiro.test(res30)
# the normality hypothesis cannot be rejected

jarque.bera.test(res30)
qqnorm(res30)

# ------
# likelihood ratio test
lrtest<-as.numeric(2*(ll30-ll20))
pchisq(lrtest,df=1,lower.tail=FALSE)
# the improvement in the log-likelihood is not significant given a p-value of the likelihood-ratio test

# ====================
##ARMA (1,1)

arma11<-arima(y,order=c(1,0,1))
arma11

# ------
armaRoots(arma11$coef[1]) # notice that we are ignoring ma component

# ------
# information criteria and log likelihood

ll11<-logLik(arma11)
aic11<-arma11$aic

# ------
# checking residuals

res11 <-residuals(arma11)

plot(res11,ylab="Residual")

Box.test(res11,lag=20,type="Ljung-Box") 
shapiro.test(res11)
jarque.bera.test(res11)
qqnorm(res11)

# ======
# Compare models 

ll20; ll30; ll11
aic20; aic30; aic11

# ARMA(1,1)-model should be favored 

# ======
# Using auto.arima()

# install.packages("forecast")

# Let us compare models using different information criteria.

library(forecast)
auto.arima(y,max.p=3,max.q=3,start.p=1,
           start.q=1,ic="aic")
auto.arima(y,max.p=3,max.q=3,start.p=1,
           start.q=1,ic="bic")
auto.arima(y,max.p=3,max.q=3,start.p=1,
           start.q=1,ic="aicc")

#------------------------------------------------------------------
# Forecasts based on ARMA  
#------------------------------------------------------------------

# We can use the function forecast() and its associated
# plot method contained in the package forecast.

# Forecasts
arma11.pred<-predict(arma11,n.ahead=10)

predict<-ts(c(rep(NA,length(y)-1),y[length(y)],
              arma11.pred$pred),start=1909,
            frequency=1)

# Forecasts confidence intervals  
upper<-ts(c(rep(NA,length(y)-1),y[length(y)],
            arma11.pred$pred+2*arma11.pred$se),
          start=1909,frequency=1)
lower<-ts(c(rep(NA,length(y)-1),y[length(y)],
            arma11.pred$pred-2*arma11.pred$se),
          start=1909,frequency=1)

# Observed
observed<-ts(c(y,rep(NA,10)),start=1909,
             frequency=1)


# Plot of actual and forecasted values
plot(observed,type="l",
     ylab="Actual and predicted values",xlab="")
lines(predict,col="blue",lty=2)
lines(lower,col="red",lty=5)
lines(upper,col="red",lty=5)
abline(v=1988,col="gray",lty=3) 
abline(h=mean(y),col="green",lty=3)

# -----------------------------------------------------------------
#                   Example 2
# -----------------------------------------------------------------

# load data: quarterly.csv
#  United States real GDP collected quarterly

y <- read.csv(file.choose(),header=TRUE)
y <- ts(y, frequency = 4,start = 1947)

# ------
plot(y)

log_y <- log(y)
plot(log_y)

# ------
acf(log_y)
acf(diff(log_y))

# differencing
growth<-log_y-lag(log_y,-1) # why is this called growth? Differencing implies the change of gdp
growth1<- diff(log_y)

cor(growth, growth1)
plot(growth)

length(growth)
length(log_y)

# ------
# inital assement of lag
acf(growth, xlim=c(0.25,6)) # data is quarterly so I have to adjust xlim
pacf(growth)

# ------
# I have differenced the series so I do not need to specify integration order, d in ARIMA(p,d,q)
arma<-arima(growth,order=c(1,0,2)) 
arma
mean(log_y) # mean of original data
mean(growth) # mean of the differenced data - this is drift

# ----------
# When fitting ARIMA models with R,constant term is NOT included in the model if there is 
# any integrations order. The following two commands should produce same AR and MA coefficients

arima(log_y,order=c(1,1,2))  # notice, I use log_y, not growth
arma  # the model from above

# to include constant include xreg=1:length(log_y)
# this estimates drift - notice, this is the mean of diff(log_y)
arma112<-arima(log_y,order=c(1,1,2), xreg=1:length(log_y)) 
arma112
round((1-pnorm(abs(arma112$coef)/sqrt(diag(arma112$var.coef))))*2,4) # significance

arma111<-arima(log_y,order=c(1,1,1), xreg=1:length(log_y))
arma111
round((1-pnorm(abs(arma111$coef)/sqrt(diag(arma111$var.coef))))*2,4) # significance

# ----------
# Our AR coeficient is statistically significant. But it is also twice as big, 
# which implies that the data are not providing the precise estimate of the 
# dynamic responses od GDP growth to economic shocks

# you can compare the dynamic responses of these two process 

irf1<-arma112$coef[1] + arma112$coef[2]
irf2<-arma112$coef[1] * irf1 + arma112$coef[3]
irf3 <- arma112$coef[1] *irf2
irf_arma112<-cbind(irf1, irf2, irf3)
rownames(irf_arma112)<-NULL

irf1<-arma111$coef[1] + arma111$coef[2]
irf2<-arma111$coef[1] * irf1 
irf3 <- arma111$coef[1] *irf2
irf_arma111<-cbind(irf1, irf2, irf3)
rownames(irf_arma111)<-NULL

irf_arma112; irf_arma111

# ----------
# I am trying two additional specifications

arma012<-arima(log_y,order=c(0,1,2), xreg=1:length(log_y))

irf1<-arma012$coef[1] 
irf2<-irf2+ arma012$coef[2]
irf3 <- 0
irf_arma012<-cbind(irf1, irf2, irf3)
rownames(irf_arma012)<-NULL

arma110<-arima(log_y,order=c(1,1,0), xreg=1:length(log_y))

irf1<-arma110$coef[1]
irf2<-arma110$coef[1] * irf1 
irf3 <- arma110$coef[1] *irf2
irf_arma110<-cbind(irf1, irf2, irf3)
rownames(irf_arma110)<-NULL

irf_arma012;irf_arma110

# ----------
arma112$aic
arma111$aic
arma012$aic
arma110$aic

# There is little change in AIC, so I select the most parsimonious model
# ----------
# checking residuals

res110 <-residuals(arma110)

plot(res110,ylab="Residual")

Box.test(res110,lag=20,type="Ljung-Box") 
acf(res110)
qqnorm(res110)
hist(res110)

# However, these tests reject the null of normality
shapiro.test(res110)
jarque.bera.test(res110)

# Why is this happening? 
# Jarque-Bera test is based on if skewness and kurtosis differ from their expected values under normality,
# so we may examine why skewness and kurtosis may be off. On the other hand Shapiro-Wilk test is biased 
# when sample sizes are samll; see this
set.seed(34567)
shapiro.test(runif(9)) 

# Lets try other models:
shapiro.test(residuals(arma112))
jarque.bera.test(residuals(arma112)) # so it is happening all the time.

# Also we have a very unusual behavior at the end of the series corresponding to the economic crisis of 2008
plot(y)
hist(res110)
# Most likely we would need to account for this regime change

# -------------------------

# Forecast 

arma110<-arima(log_y,order=c(1,1,0), xreg=1:length(log_y))
arma110

library(forecast)

arma110 <- Arima(log_y, order = c(1,1,0), include.drift=TRUE)
arma110

plot(forecast(arma110, h=90), main = "Predicted values")

###################################################################
#
#               Tests for unit roots
#
###################################################################

library(urca)

set.seed(1234)

y <- arima.sim(list(order = c(2,1,0),ar = c(0.557, -0.4858)), n = 500)
plot(y,type="l")


df=ur.df(y,type="none",lags=0)
df
summary(df) # so I am not rejecting null of process containing unit root

# ----------
# identify lag order

df=ur.df(y,type="none",lags=4)
summary(df)
df=ur.df(y,type="none",lags=3)
summary(df)
df=ur.df(y,type="none",lags=2) # we have evidence of two lags
summary(df)

# ----------
# Augmented Dickey Fuller

# so let us start with procedure

summary(ur.df(y, type = "trend",lags=2))
# test value for tau3 is -2.248, we are failing to reject presence of unit roots 
# test value for phi3 is 2.5924, we are failing to reject deterministic trend=0

summary(ur.df(y, type = "drift",lags=2))
# test value for tau2 is -1.2689, we are failing to reject presence of unit roots
# test value for phi1 is 0.8143, we are failing to reject drift=0

summary(ur.df(y, type = "none",lags=2))
# test value for tau2 is -1.1941, we are failing to reject presence of unit roots
# so we can say there is unit root in the series and series is random walk

# ----------
# KPSS test: Kwiatkowski-Phillips-Schmidt-Shin Test

ir.kpss<-ur.kpss(y,type="mu",use.lag=2) # this is level-stationary version 
summary(ir.kpss)
# we are rejecting hypothesis of stationarity around a constant mean

wg.kpss<-ur.kpss(y,type="tau",use.lag=2) # this is trend-stationary version
summary(wg.kpss)
# we are rejecting hypothesis of stationarity around deterministic trend

# ----------
# Phillips-Perron Test

y.ct<-ur.pp(y,type='Z-tau',model='trend', lags='long') # with intercept and trend
summary(y.ct)

y.co<-ur.pp(y,type='Z-tau',model='constant',lags='long') # with constant only, the trend is dropped
summary(y.co)

# ----------
# Schmidt-Phillips Test

y.tau.sp<-ur.sp(y,type="tau",pol.deg=2,signif=0.05)
summary(y.tau.sp)

y.rho.sp<-ur.sp(y,type="rho",pol.deg=2,signif=0.05)
summary(y.tau.sp)

# ----------
# Elliott-Rothenberg-Stock Test

y.df<-ur.ers(y,type="DF-GLS",model="trend",lag.max=4)
summary(y.df)

