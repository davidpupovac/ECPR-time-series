###################################################################
###################################################################
#
#                 Times Series Analysis - Day 5 
#
###################################################################
###################################################################

setwd("D:/My Documents/Desktop/ECPR data") 

#------------------------------------------------------------------
#   Bivariate and multivariate white noise
#------------------------------------------------------------------

library(mvtnorm)

# define the covariance matrix of sigma 1 and correlation of 0.8 
# (as multivariate white noise can be correlated)
set.seed(12345)
cov.mat <- matrix(c(1, 0.8, 0.8, 1), nr = 2)
w <- rmvnorm(1000, sigma = cov.mat)
cov(w)

wx <- w[, 1]
wy <- w[, 2]

plot(wx, wy, main = "")

# the cross-correlations (we will discuss this function more later)
c_cor<-ccf(wx, wy)
c_cor

#------------------------------------------------------------------
#                  Simulated VAR
#------------------------------------------------------------------
# =====
# simulated series with AR proces 

set.seed(12345)

x <- y <- rep(0, 1000)
x[1] <- wx[1]
y[1] <- wy[1]
for (i in 2:1000) {
  x[i] <- 0.6 * x[i - 1] + 0.5 * y[i - 1] + wx[i]
  y[i] <- 0.2 * x[i - 1] + 0.1 * y[i - 1] + wy[i]
}

# estimate models and select one on the basis of AIC
xy.ar <- ar(cbind(x, y)) # notice I use function ar to estimate var
xy.ar$order

xy.ar$ar 
# This output option is not the most elegant. The letters above coefficients indicate dependent variable

# To check for stationarity we check characteristic equation 

#  On the basis of the output above our matrix looks like:

# 1 - 0.5705L  -0.5362L
# -0.1932L     1 - 0.0952L

# so determinant is:
# (1-0.5705L) * (1-0.0952L) - (-0.5362 * -0.1932)L^2 = 1 - 0.6657L - 0.0492L^2

Mod(polyroot(c(1,-0.6657,-0.0492)))

###################################################################
#
#                        VAR model
# 
###################################################################

# A best fitting VAR model is fitted to the (mean-adjusted) gross national product (GNP) and real money (M1)

# load data: var_data.csv
var_data<-read.csv(file.choose(),header=TRUE)

ts.plot(var_data)

ts.plot(var_data$y1)
ts.plot(var_data$y2)

# we do not use acf and pacf in diagnostics of VAR, but we want to get some feel of the data
acf(var_data$y1)
pacf(var_data$y1)

acf(var_data$y2)
pacf(var_data$y2)

# similarly, we asses the stationarity of VAR only after estimation, but we want to understand the data
library(tseries)
adf.test(var_data$y1)
pp.test(var_data$y1)

adf.test(var_data$y2)
pp.test(var_data$y2)

# =====

library(vars)

# Determine lag order 

# This funcition may be influcenced by lag.max, so try different values
VARselect(var_data,lag.max=4)

# ======
# Estimating the model

# For estimation of the model we will use VAR() function from library(vars)
# Function estimates a VAR by utilising OLS per equation.
# You can include deterministic regressors, seasonal dummies, and exogenous variables

var_model<-VAR(var_data,p=2,type="const", season=NULL,exogen=NULL)
var_model
summary(var_model)


# Lets see if the linear regression will provide same solution (this is just one equation)
library(dplyr)
lag1_y1 <- dplyr::lag(var_data$y1,1)
lag2_y1 <- dplyr::lag(var_data$y1,2) 

lag1_y2 <- dplyr::lag(var_data$y2,1) 
lag2_y2 <- dplyr::lag(var_data$y2,2) 

summary(lm(var_data$y1~lag1_y1 + lag1_y2 + lag2_y1 + lag2_y2))

# ======
# check for stationarity

# the following function returns a vector of the eigenvalues of the companion coefficient matrix
# /these should be inside unite circle
vars::roots(var_model)
vars::roots(var_model, modulus = F) # vector of complex numbers is returned

#------------------------------------------------------------------
#   Diagnostic Tests of VAR-process
#------------------------------------------------------------------

# you can check acf for autocorrelation
acf(resid(var_model)[, 1])
acf(resid(var_model)[, 2])

# =====
# Test for autocorrelation

# The following tests have a null hypothesis of no residual autocorrelation

# Breusch-Godfrey LM test
serial.test(var_model,lags.pt=16, type="BG")
# Edgerton-Shukur F test - small sample modification of LM test
serial.test(var_model,lags.pt=16, type="ES")
# Portmanteau-Test
serial.test(var_model,lags.pt=16, type="PT.asymptotic")

# =====
# Testing for heteroscedasticity

# Null hypothesis of no heteroscedasticity

arch.test(var_model,lags.multi=5, multivariate.only=TRUE)

# tests for heteroscedasticity are the univariate and multivariate ARCH tests
# If multivariate.only = FALSE, the univariate tests are computed

# =====
# Testing for normality

# The Jarque-Bera normality tests applied to the residuals
# as well as separate tests for multivariate skewness and kurtosis

# null hypothesis of normal distribution

normality.test(var_model, multivariate.only=TRUE)

var.norm <- normality.test(var_model, multivariate.only=TRUE)
plot(var.norm)

# =====
# Teting for structural stability

# Structural stability can be tested by investigating the empirical fluctuation process

reccusum<-vars::stability(var_model, type="OLS-CUSUM")
plot(reccusum)

# A CUSUM chart is a time-weighted control chart that displays the cumulative sums (CUSUMs) of the deviations
# of each sample value from the target value. Because it is cumulative, even minor drifting in the process
# mean will lead to steadily increasing or decreasing cumulative deviation values. Therefore, this chart is
# especially useful in detecting slow shifts away from the target value

#------------------------------------------------------------------
#   Forecasting based on VAR
#------------------------------------------------------------------

# Below we give the predicted values for the next ten periods
VAR_forecast <- predict(var_model, n.ahead = 10)
VAR_forecast

# Ploting forecasts
args(vars:::plot.varprd)

plot(VAR_forecast,names="y1")
plot(VAR_forecast,names="y2", xlim=c(490, 510))


# Fancharts
args(fanchart)
# if colors and cis are left NULL, then as defaults a gray color scheme is used and
# the critical values are set from 0.1 to 0.9 with a step size of 0.1.

fanchart(VAR_forecast,names="y1")

fanchart(VAR_forecast,names="y2", 
         col.y = "gray38", main = "Fanchart",
         xlim=c(490, 510))

###################################################################
#
#             Recursive VAR - Causality analysis
#
###################################################################

#------------------------------------------------------------------
#   Cross correlations
#------------------------------------------------------------------

# load data us_macro.csv

# Quarterly data on federal funds rate, inflation and unemployment rate starting with 1960
us_macro <- read.csv(file.choose(),header=TRUE)
View(us_macro)

ccf(us_macro$unrate,us_macro$inflation)

# Negative value for lag represents a correlation between the x-series at a time before t and the y-series at time t.
# Second variable is y-series

# When lag < 0 (left side of plot), x leads y 
# significant and large ccf values imply that past values of x are correlated to the present values of  y

# When lag > 0 (right side of plot), x lags y 
# significant and large ccf values imply that future values of x are correlated to the present values of y

ccf(us_macro$unrate,us_macro$inflation,20, main = "") 
# future values of unrate are correlated to the present values of inflation

ccf(us_macro$inflation,us_macro$unrate,20, main = "") 
# past values of inflation are correlated to the present values of unrate


ccf(us_macro$inflation,us_macro$unrate,plot=FALSE)


# but the real world is seldom this clear
ccf(us_macro$inflation,us_macro$unrate,60, main = "")

#------------------------------------------------------------------
#   Granger causality
#------------------------------------------------------------------

library(vars)

plot.ts(us_macro)

# let us estimate a provisional VAR model
VARselect(us_macro,lag.max=4) 

var_us_macro<-VAR(us_macro,p=3,type="const", season=NULL,exogen=NULL)
vars::roots(var_us_macro)

summary(var_us_macro)

# =====
# Granger and instantaneous causality

causality(var_us_macro,cause="inflation")
# null hypothesis of no Granger causality has to be rejected
# null hypothesis of no Wald-type instantaneous causality cannot be rejected

causality(var_us_macro,cause="unrate")

causality(var_us_macro,cause="fedfunds")

# It would be prudent to estimate  VARs by excluding some of the variables and their lags
# and see if Granger causality can be asserted only jointly or for particular pairs of variables

#------------------------------------------------------------------
#   Impulse Response Functions
#------------------------------------------------------------------

# Impulse response analysis

# As our analysis has three variables (excluding lags) VAR consists of three equations, 
# In the code bellow impulse indicates the equation which is shocked. 
# For each equation we can trace the response all three variables. 
# Consequently, in this example we will have 9 sets of impulse responses. 


# Here IRF traces the response of unempolyment rate over time due to shocks to the inflation equation.
irf.infl1 <-irf(var_us_macro,impulse="inflation",
            response="unrate",n.ahead=40,
            ortho=FALSE,cumulative=FALSE,
            boot=TRUE,seed=12345)
plot(irf.infl1)
# Here IRF traces the response of fedfunds rate over time due to shocks to the inflation equation.
irf.infl2 <-irf(var_us_macro,impulse="inflation",
               response="fedfunds",n.ahead=40,
               ortho=FALSE,cumulative=FALSE,
               boot=TRUE,seed=12345)
plot(irf.infl2)
# Here IRF traces the response of inflation  over time due to shocks to the inflation equation.
irf.infl3 <-irf(var_us_macro,impulse="inflation",
               response="inflation",n.ahead=40,
               ortho=FALSE,cumulative=FALSE,
               boot=TRUE,seed=12345)
plot(irf.infl3)

# To get a good understanding of dynamic indicated by the model you should review all IRFs

# We can also have cumulative shocks (again there are 9 sets)
irf.infl.cum <- irf(var_us_macro,impulse="inflation",
                     response="unrate",n.ahead=40,
                     ortho=FALSE,cumulative=TRUE,
                     boot=TRUE,seed=12345)
plot(irf.infl.cum)

# The default length of the impulse responses is set to 10 via argument n.ahead.
# The computation of orthogonal and/or cumulated impulse responses is controlled by  ortho and cumulative
# Confidence bands can be returned by setting boot = TRUE (default)

#------------------------------------------------------------------
# Orthogonalized impulse response function
#------------------------------------------------------------------

# Let us interpret some  of OIRFs 

# OIRFs for federal funds rate:
oirf.fed1<-irf(var_us_macro,impulse="fedfunds",
             response="inflation",n.ahead=40,
             ortho=TRUE,cumulative=FALSE,
             boot=TRUE,seed=12345)
plot(oirf.fed1)

oirf.fed2<-irf(var_us_macro,impulse="fedfunds",
              response="unrate",n.ahead=40,
              ortho=TRUE,cumulative=FALSE,
              boot=TRUE,seed=12345)
plot(oirf.fed2)

oirf.fed3<-irf(var_us_macro,impulse="fedfunds",
               response="fedfunds",n.ahead=40,
               ortho=TRUE,cumulative=FALSE,
               boot=TRUE,seed=12345)
plot(oirf.fed3)

# According to the these graphs there is almost no effect of monetary policy (federal funds) 
# on unemployment or inflation; zero is typically within the confidence interval. However, there 
# is some initial drop in unemployment rate as a response to increase to federal funds rate
# A shock to fedfunds leads to initial increase of fedfunds rates, but subsequently leads to easing of monetary policy 

# OIRFs for inflation:
oirf.infl1<-irf(var_us_macro,impulse="inflation",
               response="inflation",n.ahead=40,
               ortho=TRUE,cumulative=FALSE,
               boot=TRUE,seed=12345)
plot(oirf.infl1)

oirf.infl2<-irf(var_us_macro,impulse="inflation",
               response="unrate",n.ahead=40,
               ortho=TRUE,cumulative=FALSE,
               boot=TRUE,seed=12345)
plot(oirf.infl2)

oirf.infl3<-irf(var_us_macro,impulse="inflation",
               response="fedfunds",n.ahead=40,
               ortho=TRUE,cumulative=FALSE,
               boot=TRUE,seed=12345)
plot(oirf.infl3)

# A shock to inflation leads to prolong high inflation period which is slowly decreasing, 
# although there is sharp drop in inflation in first few quarters. Increase in inflation generates 
# increase in unemployment rate (however, these results are significant only about 12 quarters later)
# A shock to inflation leads to relatively long increase in federal funds rate (so there is an attempt to control inflation)

#------------------------------------------------------------------
# Orthogonalized impulse response function - dependence on ordering
#------------------------------------------------------------------

irf(var_us_macro,impulse="inflation",
    response="fedfunds",n.ahead=10,
    ortho=TRUE,cumulative=FALSE,
    boot=TRUE,seed=12345)

# Notice that the first value is 0. This is because the ordering . Inflation is second variable/equation in our system 
# Therefore, shocks to the second equation have contemporaneous impact on all equations except the first one

# =====
# Now we are going to change the order of the variables in the data set
# This implies that we are chaning the order of equations in VAR 
# (actually this ordering is more in line with Granger tests above)

us_macro_2<- data.frame(us_macro$inflation, us_macro$unrate, us_macro$fedfunds)
colnames(us_macro_2) <- c("inflation","unrate","fedfunds")

# We are estimating the model with same specifications, but different order 
var_us_macro_2<-VAR(us_macro_2,p=3,type="const", season=NULL,exogen=NULL)

# =====

oirf.fed1_2<-irf(var_us_macro_2,impulse="fedfunds",
               response="inflation",n.ahead=40,
               ortho=TRUE,cumulative=FALSE,
               boot=TRUE,seed=12345)
plot(oirf.fed1_2)
plot(oirf.fed1)

#------------------------------------------------------------------
#   FEVD of VAR-process
#------------------------------------------------------------------

# Forecast error variance decomposition
fevd.var<-fevd(var_us_macro,n.ahead=10)

# The argument n.ahead sets the number of forecasting steps

fevd.var 
# Notice how the results are partly due to ordering

# Table 1. FEV in fedfunds is predominantly due to uncertainty in fedfunds equation,
# but it is declining fast with additional quarters 
# Actually, already at step 8 ahead approximately 25% of FEV is due to uncertainty in unrate equation

# Table 2. FEV in inflation is consistently mainly due to uncertainty in inflation equation

# Table 3. FEV in unrate due to uncertainty in unrate equation grow with successive quarters.
# So, error variance in forecasting unemployment is almost fully due to unrate equation

plot(fevd.var,addbars=2)


# ======

# Keep in mind that in this example there are !n=6 posible orderings 
# Hopefully you can discard some on theoretical grounds, otherwise you will analyze 54 OIRFs, 54 FEVDs...

###################################################################
#
#                     Error-correction model
#
###################################################################

# spurious correlation due to stochastic trends in both series happen to be coincident 

# (run it run it couple of time excluding seed)

set.seed(10)

x <- rnorm(100); y <- rnorm(100)

for(i in 2:100) {
  x[i] <- x[i-1] + rnorm(1)
  y[i] <- y[i-1] + rnorm(1)}

plot(x, y)
ts.plot(y)


summary(lm(x~ y))
cor(x,y)

# Granger and Newbold (1974):
# if your R^2 is larger than your Durbin Watson statistics you have a spurious regression

library(lmtest)
dwtest(lm(x~ y))

#------------------------------------------------------------------
#   Engle-Granger procedure with generated data
#------------------------------------------------------------------

set.seed(123456)

e1<-rnorm(100)
e2<-rnorm(100)
y1<-cumsum(e1) # random walk
y2<-0.6*y1+e2

ts.plot(as.ts(y1),as.ts(y2), col=c("red","darkgreen"))

# y2, has been set to 0.6*y1 + e2, where e2 is a white noise process
# Hence, the cointegrating vector is (1,-0.6).

# Although we know there is cointegration let's test it with Phillips and Ouliaris test.
summary(ca.po(cbind(y1,y2),demean='none', type='Pu')) # variance ratio test
summary(ca.po(cbind(y1,y2),demean='none', type='Pz')) # multivariate trace statistic
# The null hypothesis is that no cointegration relationship exists.
# The inclusion of deterministic regressors can be set via the argument demean, 
# Lag length for estimating the long-run variance-covariance matrix can be set with the argument lag.

lr.reg<-lm(y2~y1) # estimate long-run equation
lr.reg

lines(predict(lr.reg)) # plot cointegrating relation

error<-residuals(lr.reg) # equilibrium error is stored

plot(error, type="l") # plot cointegrating relation

# If this was a regular (not simulated) analysis you would do augmented 
# Dickey-Fuller (ADF)-type test and Jarque-Bera test of normality

#-----------------------------------------------------------------

error.lagged <-error[-c(1,100)] 
# you have to take two (not one) values out, because of differencing and lagging of differences

dy1<-diff(y1)
dy2<-diff(y2)

diff.dat<-data.frame(embed(cbind(dy1,dy2),2))

colnames(diff.dat)<-c('dy1','dy2','dy1.1','dy2.1')
ecm.reg<-lm(dy2 ~ error.lagged + dy1 +  dy1.1 + dy2.1, data=diff.dat) # estimate lamda

summary(ecm.reg)

#-----------------------------------------------------------------
# Basic form, which is not adequate here would look like this:

ecm.reg<-lm(dy2 ~ error.lagged + dy1 , data=diff.dat)
summary(ecm.reg)


#----

lr.reg<-lm(y1~y2) 
lr.reg
error<-residuals(lr.reg) 

ts.plot(as.ts(y1),as.ts(y2), col=c("red","darkgreen"))
lines(predict(lr.reg)) 

error.lagged <-error[-c(1,100)] 
ecm.reg<-lm(dy1 ~ error.lagged + dy1.1 + dy2.1 + dy2 , data=diff.dat) # estimate lamda

summary(ecm.reg)



#-----------------------------------------------------------------


# You can also repeat procedure by estimating lm(y1~y2) - this  has more sense when there are more than two variables
# If two series are cointegrated, then there should be Granger-causation in at least one direction:
# at least one coefficient of the error term should be significant and with the correct sign (i.e., negative)

###################################################################
# 
#                Vector error-correction model
#
###################################################################

#------------------------------------------------------------------
# Johansen method with generated data
#------------------------------------------------------------------

# simulate the data
library(urca)
set.seed(12345)
e1<-rnorm(250,0,0.5)
e2<-rnorm(250,0,0.5)
e3<-rnorm(250,0,0.5)

u1.ar1<-arima.sim(model=list(ar=0.75),
                  innov=e1,n=250)
u2.ar1<-arima.sim(model=list(ar=0.3),
                  innov=e2,n=250)

y3<-cumsum(e3)
y1<-0.8*y3+u1.ar1
y2<--0.3*y3+u2.ar1
y.mat<-data.frame(y1,y2,y3)


ts.plot((y.mat), col =c("red","darkgreen","black"))

# ------
# Determining the cointegration rank

summary(ca.jo(y.mat,type='trace',K=2)) # The default value for the test statistic is the maximum eigenvalue test.
summary(ca.jo(y.mat,type='eigen',K=2))

# ------
# You should also select the appropriate number of lags, which can be done by using function VARselect from package vars.
# You should also choose whether constant, or trend terms are estimated besides cointegration relationship or not. 

# ca.jo(y.mat,type='eigen', ecdet = c("none", "const", "trend"),K=2)

# ------
# The hypothesis of two cointegrating vectors cannot be rejected for all significance levels.
# (Test statistics should be lower than critical values of the test)
# if deterministic term, constant, or trend should be included in cointegration is set by ecdet =
# model exogenous regressors can be provided by setting dumvar

model <- ca.jo(y.mat,type='eigen',K=2)

model@V   # - this is Beta matrix
# Eigenvectors, normalised to first column: (These are the cointegration relations) 

model@W   # - this is Alpha matrix 
# The loading or adjustment matrix - used to determine speed of adjustment to the long-run equilibrium 

# cointegrating vectors are all normalized to the first variable, 
# and the loading matrix alphA is adjusted accordingly

# =====
# Obtaining coefficients:

# Johansen [1995] proposed restricting beta by introduction of identity matrix in beta.
# These identifying restrictions and the estimation of the 
# VECM coefficients are implemented in the function cajorls()

vecm.r2<-cajorls(model,r=2)
vecm.r2

vecm.r2$beta   # normalized cointegration vectors 
# Notice I (identity) matrix at the top. 
# There are r=cointegrating relatioships, so there are r^2 restrictions

vecm.r2$rlm   
# First two rows are error-correction terms
# Second row are constants
# Third row are lagged differences

# =====
# Diagnostic tests:

# VECM as level-VAR representation 
library(vars)
vecm.level<-vec2var(model,r=2)

# Diagnostic tests
arch.test(vecm.level)
normality.test(vecm.level)
serial.test(vecm.level)

Forecasting 

# =====
# Forecast:

predict(vecm.level)

# IRF and FEVD
irf(vecm.level,boot=FALSE)
fevd(vecm.level)

