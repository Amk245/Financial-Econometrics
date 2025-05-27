# The following code shows how to estimate a bivariate CCC
# model using the equation by equation approach

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(yfR)
source("llik_fun_GARCH.R")

########################
############ 1. Load data
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

########################
############ 2. Estimate univariate GARCH for each of the series
########################

## 2a. Estimate a univariate GARCH(1,1) for Google log-returns x[,1]

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

########################
############ 3. Estimate the costant correlation from the residuals
########################

## 3a. Obtain estimated variances

n <- length(x[,1])

s1 <- rep(0,n)
s2 <- rep(0,n)

s1[1] <- var(x[,1])
s2[1] <- var(x[,2])

# we run the univariate GARCH recursions evaluated at estimated parameter theta_hat

for(t in 2:n){
  s1[t] <- omega_hat1 + alpha_hat1*x[t-1,1]^2 + beta_hat1*s1[t-1]
  s2[t] <- omega_hat2 + alpha_hat2*x[t-1,2]^2 + beta_hat2*s2[t-1]
}

## 3b. Obtain standardized series

e1 <- x[,1]/sqrt(s1)
e2 <- x[,2]/sqrt(s2)

## 3c. Estimate the correlation of the residuals

r <- cor(e1,e2)

## 3d. Obtain the conditional covariance between Microsoft and IBM

s12=r*sqrt(s1)*sqrt(s2)

########################
############ 4. Display estimation results
########################

# conditional correlation:
r

# w1, a1, b1 from variance equation 1st series:
(theta_hat1 <- c(omega_hat1,alpha_hat1,beta_hat1))

# w2, a2, b2 from variance equation 2nd series:
(theta_hat2 <- c(omega_hat2,alpha_hat2,beta_hat2))

########################
############ 5. Plot estimated conditional variances, covariance and correlation
########################

par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(s1,type="l",main = "Variance google",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s2,type="l",main = "Variance IBM",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s12,type="l",main = "Covariance",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s12/(sqrt(s1)*sqrt(s2)),type="l",main = "Correlation",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
