###  SIMULATE DATA from a bivariate SV MODEL
#
#  Description: 
#  This code shows how to simulate data from a bivariate SV model
#  For a more detailed description of the code see the Lecture Notes

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(MASS) 

########################
############ 1. Choose sample size and parameter values
########################

n <- 2500          # set the sample size n

# Parameter Values:

omega1 <- 0      # intercept parameter in transition equation
omega2 <- 0      

beta1 <- 0.95     # autoregressive parameter in transition equation
beta2 <- 0.95     

sig2f1 <-  0.10  # covariances transition equation error
sig2f2 <-  0.10  
sigf12 <-  0.05  

rho <- 0.5        # correlation observation equation error

## Define the correlation matrix R of epsilon_t

R <- cbind(c(1,rho),c(rho,1))

# Define the covariance matrix of eta_t

Sf <- cbind(c(sig2f1,sigf12),c(sigf12,sig2f2))

########################
############ 2. Simulate from SV
########################

## 2a. Generate Innovations

epsilon <- mvrnorm(n,rep(0,2),R)  #generate errors observation equation
eta <- mvrnorm(n,rep(0,2),Sf)     #generate errors transition equation

## 2b. Define Time Series Matrices

x <- matrix(0,nrow=n,ncol=2) # Define a matrix of dimension nx2 for yt
f <- matrix(0,nrow=n,ncol=2) # Define a matrix of dimension nx2 for ft

## 2c. Initialize the process ft

# generate ft at t=1 from its unconditional distribution

umf <- c(omega1/(1-beta1), omega2/(1-beta2))  # unocnditional mean ft
uSf <- matrix(0,nrow=2,ncol=2)                # unocnditional variance ft
uSf[1,1] <- sig2f1/(1-beta1^2)
uSf[2,2] <- sig2f2/(1-beta2^2)
uSf[2,1] <- sigf12/(1-beta1*beta2)
uSf[1,2] <- sigf12/(1-beta1*beta2)

f[1,] <- mvrnorm(1,umf,uSf)

## 2d. Genrate first observation yt

x[1,1] <- exp(f[1,1]/2) * epsilon[1,1] 
x[1,2] <- exp(f[1,2]/2) * epsilon[1,2] 

## 2e. Generate Time Series

#generate from a multivariate SV using observation equation and
#transition equation

for(t in 2:n){ # start recursion from t=2 to t=n
  
  f[t,1] <- omega1 + beta1*f[t-1,1] + eta[t,1] # transition equation
  f[t,2] <- omega2 + beta2*f[t-1,2] + eta[t,2]
  
  x[t,1] <- exp(f[t,1]/2) * epsilon[t,1] # observation equation
  x[t,2] <- exp(f[t,2]/2) * epsilon[t,2]

}

#obtain the dynamic covariance

s12 <- rho*exp(f[,1]/2)*exp(f[,2]/2)


########################
############ 3. Plot the simulated yt
########################

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1)) 
plot(x[,1], type="l", main="simulated y1t",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(x[,2], type="l", main="simulated y2t",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

########################
############ 4. Plot the simulated SIGMAt
########################

par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1)) 
plot(exp(f[,1]), type="l", main="simulated sigma2_1t",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(exp(f[,2]), type="l", main="simulated sigma2_2t",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s12, type="l", main="simulated sigma12t",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s12/(exp(f[,1]/2)*exp(f[,2]/2)), type="l", main="simulated rho12",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
