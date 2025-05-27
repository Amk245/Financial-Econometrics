# This R file contains code for the estimation of the 
# bivariate scalar sDVECH(1,1) model
# Stock returns of Google and IBM are considered

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(yfR)
source("llik_fun_sDVECH.R") 

########################
############ 1. Obtain returns of Google and IBM
########################

# download Google stock prices from yahoo
gog_data <- yf_get(tickers = "GOOG", 
                   first_date = "2005-01-01",
                   last_date = "2018-01-01",freq_data = "daily")
# download IBM stock prices from yahoo
ibm_data <- yf_get(tickers = "IBM", 
                   first_date = "2005-01-01",
                   last_date = "2018-01-01",freq_data = "daily")

# Obtain Google stock prices  
p_gog <- gog_data$price_adjusted
# Obtain IBM stock prices  
p_ibm <- ibm_data$price_adjusted

# obtain log-returns for Google
r_gog <- diff(log(p_gog))*100
# obtain log-returns for IBM
r_ibm <- diff(log(p_ibm))*100

#combine the two series of returns in a single matrix
x <- cbind(r_gog,r_ibm)

#obtain the sample size
n <- length(r_gog)

########################
############ 2. plot prices and log-returns
########################

par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(p_gog,type="l",main = "Prices Google",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(r_gog,type="l",main = "log-returns Google",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(p_ibm,type="l",main = "Prices IBM",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(r_ibm,type="l",main = "Log-returns IBM",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

########################
############ 3. Initial paramter value for optimization
########################

a <- 0.2 
b <- 0.6  
w11 <- var(x[,1])*(1-a-b)
w22 <- var(x[,2])*(1-a-b)
w12 <- cov(x[,1],x[,2])*(1-a-b)

par_ini <- c(log(w11),w12,log(w22),log(a/(1-a)),log(b/(1-b)))

########################
############ 4. Optimize the log-likelihood function llik_fun_sDVECH()
########################

est <- optim(par=par_ini,fn=function(par)-llik_fun_sDVECH(par,x), method = "BFGS")

########################
############ 5. Display estimation results 
########################

# parameter values
(w11_hat <- exp(est$par[1]))
(w12_hat <- est$par[2])
(w22_hat <- exp(est$par[3]))
(a_hat <- exp(est$par[4])/(1+exp(est$par[4])))
(b_hat <- exp(est$par[5])/(1+exp(est$par[5])))

cat("log likelihood value:")
-est$value*n

########################
############ 6. Obtain the estimated conditional covariance matrix
########################

## Create nx3 matrix
VECHt <- matrix(0,nrow=n,ncol=3)

# set the initial value of the conditional variance equal to the sample
# covariance

C <- cov(x)
VECHt[1,] <- c(C[1,1],C[1,2],C[2,2])

for(t in 2:n){
  VECHt[t,1] <- w11_hat+b_hat*VECHt[t-1,1]+a_hat*x[t-1,1]^2
  VECHt[t,3] <- w22_hat+b_hat*VECHt[t-1,3]+a_hat*x[t-1,2]^2
  VECHt[t,2] <- w12_hat+b_hat*VECHt[t-1,2]+a_hat*x[t-1,1]*x[t-1,2]
}

########################
############ 7. Conditional variance matrix plots
########################

# we obtain the conditional standard deviations for IBM and Apple 
# log-returns and their conditional correlation

sd1 <- sqrt(VECHt[,1])
sd2 <- sqrt(VECHt[,3])
corr <- VECHt[,2]/(sd1*sd2)

# plot the results

par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(sd1,type="l",main = "Conditional sd google",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(sd2,type="l",main = "Conditional sd IBM",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(corr,type="l",main = "Conditional correlation",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
