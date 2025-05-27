## The following code shows how to obtain optimal portfolio weights 
## in terms of Sharpe Ratio. In this case we use a CCC model to estimate the
## conditional covariance matrix

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(yfR)
library(quadprog)

source("llik_fun_GARCH.R")
source("max_SR_portfolio.R")

########################
############ 1. Load data
########################

# download Microsoft stock prices from yahoo
mic_data <- yf_get(tickers = "MSFT", 
                   first_date = "1990-01-01",
                   last_date = "2017-01-01",freq_data = "monthly")
# download IBM stock prices from yahoo
ibm_data <- yf_get(tickers = "IBM", 
                   first_date = "1990-01-01",
                   last_date = "2017-01-01",freq_data = "monthly")

# Obtain Microsoft stock prices  
p_mic <- mic_data$price_adjusted
# Obtain IBM stock prices  
p_ibm <- ibm_data$price_adjusted

# obtain log-returns for Microsoft
r_mic <- diff(log(p_mic))*100
# obtain log-returns for IBM
r_ibm <- diff(log(p_ibm))*100

#combine the two series of returns in a single matrix
x <- cbind(r_mic,r_ibm)

########################
############ 2. Estimate bivariate CCC model
########################

## 2a. Estimate a univariate GARCH(1,1) for Microsoft log-returns x[,1]

alpha_ini <- 0.2 
beta_ini <- 0.6  
omega_ini <- var(x[,1])*(1-alpha_ini-beta_ini) 
par_ini <- c(log(omega_ini),log(alpha_ini/(1-alpha_ini)),log(beta_ini/(1-beta_ini)))

est1 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x[,1]), method = "BFGS")

(omega_hat1 <- exp(est1$par[1]))
(alpha_hat1 <- exp(est1$par[2])/(1+exp(est1$par[2])))
(beta_hat1 <- exp(est1$par[3])/(1+exp(est1$par[3])))

## 2b. Estimate a univariate GARCH(1,1) for IBM log-returns x[,2]

alpha_ini <- 0.2 
beta_ini <- 0.6  
omega_ini <- var(x[,2])*(1-alpha_ini-beta_ini) 
par_ini <- c(log(omega_ini),log(alpha_ini/(1-alpha_ini)),log(beta_ini/(1-beta_ini)))

est2 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x[,2]), method = "BFGS")

(omega_hat2 <- exp(est2$par[1]))
(alpha_hat2 <- exp(est2$par[2])/(1+exp(est2$par[2])))
(beta_hat2 <- exp(est2$par[3])/(1+exp(est2$par[3])))

## 2c. Obtain estimated conditional variances and covariance

n <- length(x[,1])
s1 <- rep(0,n)
s2 <- rep(0,n)

s1[1] <- var(x[,1])
s2[1] <- var(x[,2])

for(t in 2:n){
  s1[t] <- omega_hat1 + alpha_hat1*x[t-1,1]^2 + beta_hat1*s1[t-1]
  s2[t] <- omega_hat2 + alpha_hat2*x[t-1,2]^2 + beta_hat2*s2[t-1]
}

e1 <- x[,1]/sqrt(s1)
e2 <- x[,2]/sqrt(s2)
r <- cor(e1,e2)
s12=r*sqrt(s1)*sqrt(s2)

########################
############ 3. Optimal portfolio weights (max Sharpe Ratio)
########################

# Create a nx2 matrix that will contain the 

kt = matrix(0,nrow=n,ncol=2)

# set the conditional mean of the log-returns equal to the sample mean.
# Note that mu1 is the sample mean of Microsft log-returns instead mu2
# is the sample mean of IBM log returns.

mu1 <- mean(x[,1])
mu2 <- mean(x[,2])
mut <- cbind(mu1,mu2)

# Obtain optimal weights at each t (from 1 to T) using the function max_SR_portfolio()
# This function is contained in the file max_SR_portfolio.R
# max_SR_portfolio() requires as first argument the vector of expected mean
# and as second argument the covariance matrix.
# The function then returns the optimal weights that maximize the Sharpe Ratio.

for(t in 1:n){
  
  # Obtain the conditional covariance matrix. Note that the vector "s1" 
  # contains the conditional variance of Microsft, "s2" the conditional 
  # variance of IBM and "s12" the conditional covariance

  SIGMAt <- cbind(c(s1[t],s12[t]),c(s12[t],s2[t]))
  
  # Obtain the optimal weights at time t
  kt[t,] <- max_SR_portfolio(mut,SIGMAt)

}

########################
############ 3. Plot the portfolio weights
########################

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(kt[,1],type="l",main = "Portfolio weight on Microsoft",ylim=c(0,1),ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(kt[,2],type="l",main = "Portfolio weight on IBM",ylim=c(0,1),ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
