##  Conditional Variance, Residuals and VaR
# this file contains the R code to estimate the Conditional Variance
# run some diagnostic tests on the residuals and obtain the VaR

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(tseries)
source("Financial Econometrics week 2 R code llik_fun_GARCH.R")

########################
############ 1. Load data
########################

x <- scan("Financial Econometrics week 2 R code stock_returns.txt")

########################
############ 2. Estimate model (See the file Estimate_ML_GARCH.m for explanations)
########################

a <- 0.2 
b <- 0.6  
omega <- var(x)*(1-a-b) 

par_ini <- c(log(omega),log(a/(1-a)),log(b/(1-b)))

est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x), method = "BFGS")

omega_hat <- exp(est$par[1])
alpha_hat <- exp(est$par[2])/(1+exp(est$par[2]))
beta_hat <- exp(est$par[3])/(1+exp(est$par[3]))

########################
############ 3. Obtain and plot filtered volatility
########################

## 3a. Obtain filtered volatility

# In the following we obtain an estimate of the conditional variance. It is
# very simple we just need to run the GARCH recursion for the conditional 
# varaince and plug in the estimated parameters

n <- length(x)
sigma2 <- rep(0,n)
sigma2[1] <- var(x)

for(t in 2:n){
  sigma2[t] = omega_hat + alpha_hat*x[t-1]^2 + beta_hat*sigma2[t-1]
}

## 3b. Plot filtered volatility

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(x,type="l", main="Time series",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(sigma2,type="l",col=2, main="Estimated conditional variance",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

########################
############ 4. Diagnostic Analysis Residuals
########################

## 4a. Calculate residuals 

# Obtain the residuals
u=x/sqrt(sigma2)

## 4b. Normality test 

## Perform the Jarque Bera test with R function jarque.bera.test()

jarque.bera.test(u)

# 4c. Plot squared residuals and ACF squared residulas 

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1)) 
acf(u^2,main="")
title("ACF squared residuals", line = 0.3)
pacf(u^2, main="")
title("PACF squared residuals", line = 0.3)

########################
############ 5. Value-at-Risk
########################

## 5a. Calculate VaR at 10%, 5% and 1% level

VaR10 <- qnorm(0.1,0,sqrt(sigma2))
VaR05 <- qnorm(0.05,0,sqrt(sigma2))
VaR01 <- qnorm(0.01,0,sqrt(sigma2))

## 5b. Plot VaR 

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(x,type="l", main="Time series",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(VaR10,type="l",col="gray80", main="VaR at 10%, 5% and 1% levels",
     ylim=c(min(VaR01),max(VaR10)),ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
lines(VaR05,type="l",col="gray60")
lines(VaR01,type="l",col="gray40")