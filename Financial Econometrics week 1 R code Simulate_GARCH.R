### SIMULATE DATA FROM GARCH(1,1) MODEL
#
#  Description: 
#  This code shows how to simulate data from a generalized 
#  autoregressive conditional hetereoeskedasticity (GARCH) model given by:
#  x(t) = sqrt(sig(t)) * eps(t)
#
#  sig(t+1) = omega + alpha * x(t)^2 + beta * sig(t)
#
#  where eps ~ NID(0,1)

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 1. Choose sample size and parameter values
########################

n <- 1000          # set the sample size n

omega <- 0.1     # intercept parameter in update equation
alpha <- 0.2     # OD parameter in update equation
beta <-  0.75    # autoregressive parameter in update equation

########################
############ 2. Simulate from GARCH(1,1)
########################

## 2a. Generate Innovations
# The function rnorm(n) generates a vector of length n of independent N(0,1)

epsilon <- rnorm(n)   

## 2b. Define Time Series Vector

sig2 <- rep(0,n) # define vector of zeros of length n for sigma_t^2
x <- rep(0,n) # define vector of zeros of length n for y_t

## 2c. Define Initialization for conditional varaince

sig2[1] <- omega/(1-alpha-beta) 

## 2d. Genrate first observation

x[1] <- sqrt(sig2[1]) * epsilon[1] 

## 2e. Generate Time Series

for(t in 2:n){ # start recursion from t=2 to t=n
  
  sig2[t] <- omega + alpha * x[t-1]^2 + beta * sig2[t-1] # updating equation
  
  x[t] <- sqrt(sig2[t]) * epsilon[t] # observation equation
  
} # end recursion

########################
############ 3. Plot the simulated data
########################

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1)) 
plot(x, type="l", main="simulated yt",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(sqrt(sig2), type="l", main="simulated sigmat",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
