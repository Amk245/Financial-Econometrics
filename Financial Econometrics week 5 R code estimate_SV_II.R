### Estimate an SV model by II with R
# For a more detailed description of the code see the Lecture Notes

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(moments)
source("sim_m_SV.R")
source("filter_SV.R")

########################
############ 1. Obtain S&P500 weekly log-returns
########################

#The file ^GSPC.txt contains weekly prices of the S&P500 from 2005 to 2015

sep_prices <- scan("^GSPC.txt")
x <- diff(log(sep_prices))

########################
############ 2. Obtain sample moments
########################

n <- length(x)
xa <- abs(x)
sample_m <- c(var(x), kurtosis(x), cor(xa[2:n],xa[1:(n-1)]))

########################
############ 3. Generate errors for the simulations from the SV model
########################

set.seed(123)
H <- 50*n
epsilon <- rnorm(H)  
eta <- rnorm(H)  
e <- cbind(epsilon,eta)

########################
############ 4. Choose initial parameter values for the optimization
########################

#choose the initail parameter values for the numerical optimization

b <- 0.90
sig2f <- 0.1
omega <- log(var(x))*(1-b)

par_ini <- c(omega,log(b/(1-b)),log(sig2f))

########################
############ 5. Obtain parameter estimates
########################

# 3a. Minimize criterion function

est <- optim(par=par_ini,fn=function(par) mean((sim_m_SV(e,par)-sample_m)^2), method = "BFGS")

# 3a. Obtain parameter extimate using the link functions

omega_hat <- est$par[1]
beta_hat <- exp(est$par[2])/(1+exp(est$par[2]))
sig2f_hat <- exp(est$par[3])

theta_hat <- c(omega_hat,beta_hat,sig2f_hat)
cat("The parameter estimates are:")
round(theta_hat,4)

########################
############ 6. Obtain filtered volatility
########################

f <- rep(0,n)
f[1] <- log(var(x))

#Run a for loop minizing the function filter_SV() at each time period to
#obtain the filtered estimate of ft 

for(t in 2:n){ # start recursion from t=2 to t=T
  ft_ini <- f[t-1]
  f_est <- optim(par=ft_ini ,fn= function(ft) filter_SV(x[t],ft,f[t-1],theta_hat), method = "BFGS")
  f[t] <- f_est$par
}

########################
############ 7. Plot series and etimated volatility
########################

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(x,type="l", main="Time series",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(exp(f/2),type="l",col=2, main="Filtered sigmat",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

