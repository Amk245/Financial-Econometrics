###  SIMULATE DATA FROM SV MODEL
#
#  Description: 
#  This code shows how to simulate data from a stochastic volatility model given by
#
#  x(t) = exp(f(t)/2) * eps(t)
#
#  f(t+1) = omega + beta * f(t) + eta(t)
#
#  where eps ~ NID(0,1) and eta ~ NID(0,sig2f)
#

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 1. Choose sample size and parameter values
########################

n <- 2500          # set the sample size n

omega <- 0     # intercept parameter in updating equation
beta <-  0.95  # autoregressive parameter in transition equation
sig2f <- 0.3  # error variance in transition equation

########################
############ 2. Simulate from SV
########################

## 2a. Generate Innovations
# The function rnorm(n) generates a vector of length n of independent N(0,1)

epsilon <- rnorm(n) # generate a vector of standard normal  
eta <- sqrt(sig2f)*rnorm(n) # generate a vector of normal with mean 0 and variance sig2f

## 2b. Define Time Series Vectors

f <- rep(0,n) # define vector of zeros of length n for ft
x <- rep(0,n) # define vector of zeros of length n for y_t

## 2c. Initialize the process ft

# generate ft at t=1 from its unconditional distribution, which is Normal
# with mean equal to omega/(1-beta) and variance sig2f/(1-beta^2)

f[1] <- omega/(1-beta) + sqrt(sig2f/(1-beta^2)) * rnorm(1)

## 2d. Genrate first observation

x[1] <- exp(f[1]/2) * epsilon[1] 

## 2e. Generate Time Series

for(t in 2:n){ # start recursion from t=2 to t=n
  
  f[t] <- omega + beta * f[t-1] + eta[t] # transition equation
  x[t] <- exp(f[t]/2) * epsilon[t]  # observation equation
  
} # end recursion

########################
############ 3. Plot the simulated data
########################

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1)) 
plot(x, type="l", main="simulated yt",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(exp(f/2), type="l", main="simulated sigmat",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

########################
############ 4. Plot autocorrelation functions
########################

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1)) 
acf(x^2,main="")
title("ACF squared yt", line = 0.3)
pacf(x, main="")
title("ACF yt", line = 0.3)