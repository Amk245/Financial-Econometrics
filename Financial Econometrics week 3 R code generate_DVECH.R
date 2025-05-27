## The following code can be used to simulate from a bivariate DVECH model with p=q=1

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(MASS) 

########################
############ 1. Set sample size and parameter values
########################

n <- 1000

w11 <- 0.1
w22 <- 0.2
w12 <- 0.02

b11 <- 0.7
b22 <- 0.7
b12 <- 0.7

a11 <- 0.2
a22 <- 0.2
a12 <- 0.15

########################
############ 2. Set initial conditions for the conditional variance matrix
########################

# define a nx2 matrix "x" that will contain the generated series and a nx3
# matrix "VECHt" that will contain the conditional covariance matrix in
# VECH form. n is the time series length.

x <- matrix(0,nrow = n, ncol = 2)
VECHt <- matrix(0,nrow=n,ncol=3)

# Set the unconditional covariance matrix as initialization

VECHt[1,1] <- w11/(1-b11-a11)
VECHt[1,3] <- w22/(1-b22-a22)
VECHt[1,2] <- w12/(1-b12-a12)

# Generate first observation from multivariate normal

SIGMAt <- cbind(c(VECHt[1,1],VECHt[1,2]),c(VECHt[1,2],VECHt[1,3]))
x[1,] <- mvrnorm(1,rep(0,2),SIGMAt)

########################
############ 3. Genrate from a DVECH model 
########################

#the following loop generates recursively from the DVECH model

for(t in 2:n){
  VECHt[t,1] <- w11 + b11*VECHt[t-1,1] + a11*x[t-1,1]^2 # updating equation variance series 1
  VECHt[t,3] <- w22 + b22*VECHt[t-1,3] + a22*x[t-1,2]^2 # updating equation variance series 2
  VECHt[t,2] <- w12 + b12*VECHt[t-1,2] + a12*x[t-1,1]*x[t-1,2] # updating equation covariance between series 1 and series 2

  SIGMAt <- cbind(c(VECHt[t,1],VECHt[t,2]),c(VECHt[t,2],VECHt[t,3]))
  x[t,] <- mvrnorm(1,rep(0,2),SIGMAt) # generate from the observation equation
}

########################
############ 4. Plot series 
########################

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(x[,1],type="l",main = "Series y1t",ylab="",xlab="")
plot(x[,2],type="l",main = "Series y2t",ylab="",xlab="")

########################
############ 5. Plot conditional covariance
########################

par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(VECHt[,1],type="l",main = "Conditional variance sigma1t",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(VECHt[,3],type="l",main = "Conditional variance sigma2t",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(VECHt[,2],type="l",main = "Conditional covariance sigma12t",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(VECHt[,2]/(sqrt(VECHt[,1])*sqrt(VECHt[,3])),type="l",main = "Conditional correlation rho12t",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
