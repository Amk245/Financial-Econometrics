########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

source("llik_OD_regression.R") 

########################
############ 1. Load data
########################

# Load the variables xt and yt that are contained in the files "xt.txt" and "yt.txt"

xt <- scan("xt.txt")
yt <- scan("yt.txt")

########################
############ 2. Set initialization for numerical optimization
########################

## 2a. choose the initail parameter values for the Newton Raphson optimization

a <- 0.2/sd(xt*yt) # initial value for alpha
phi <- 0.9  # initial value for beta
omega <- (cov(xt,yt)/var(xt))*(1-phi) # initial value for omega
sig2 <- var(yt)

## 2b. Transform intitial values using the inverse of the link funtions

par_ini <- c(omega,log(phi/(1-phi)),log(a),log(sig2))

########################
############ 3. Optimize Log Likelihood function
########################

est <- optim(par=par_ini,fn=function(par)-llik_OD_regression(yt,xt,par), method = "BFGS")

########################
############ 4. Obtain parameter estimates
########################

omega_hat <- est$par[1]
phi_hat <- exp(est$par[2])/(1+exp(est$par[2]))
alpha_hat <- exp(est$par[3])
sigma2_hat <- exp(est$par[4])

theta_hat <- c(omega_hat,phi_hat,alpha_hat,sigma2_hat)
cat("The parameter estimates are:")
round(theta_hat,4)

########################
############ 5. Obtain estimated betat
########################

n <- length(xt)
beta <- rep(0,n)
beta[1] <- omega_hat/(1-phi_hat) # initialize betat at unconditional expectation

for(t in 2:n){
  beta[t] <- omega_hat+phi_hat*beta[t-1]+alpha_hat*(yt[t-1]-beta[t-1]*xt[t-1])*xt[t-1];
}

########################
############ 6. Plot results
########################

par(mfrow=c(3,1),mar=c(4.1,4.1,1.1,2.1))
plot(yt,type="l", main="Series yt",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(xt,type="l", main="Series xt",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(beta,type="l", col=2,main="betat",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
