##  MAXIMUM LIKELIHOOD ESTIMATE PARAMETERS OF GARCH
#
#  Description: 
#  This code snippet shows how to optimize the log likelihood
#  and estimate the parameters of a GARCH(1,1) model by maximum likelihood.
#  Note that the log-likelihood is contained in the R file "llik_fun_GARCH.R"

########################
############ 00. Clean workspace
########################

# clear workspace
rm(list = ls())
graphics.off()

# Load package with custom functions for creating tables, 
# managing working directory and saving plots. (not required for running the notebook)

#devtools::load_all()
########################
############ 0. Load packages and/or functions
########################

# The average log-likelihood function for the GARCH(1,1) is llik_fun_GARCH(). 
# This function is contained in the R file "llik_fun_GARCH.R". 
# This file has to be saved in the working directory (getwd()).
# You can use: 

# setwd("C:\\Users\\Adnan\\Documents\\School; VU\\Bachelor; Econometrics and Operations Research\\Year 3.5; Financial Econometrics")

# C:/Users/R code will be replaced with
# the path of the folder in your PC
# We load the likelihood function with the command source(). 

source("Financial Econometrics week 2 R code llik_fun_GARCH.R") 

# The same as above applies for the function Hess_fun_GARCH()
# This function will be useful to obtain the covariance matrix of the ML estimates

source("Financial Econometrics week 2 R code Hess_fun_GARCH.R")

########################
############ 1. Load data
########################

# The command "scan" loads the txt file "stock_returns" which is stored in "x". 
# The file "stock_returns.txt" has to be saved in the working directory (getwd()).

x <- scan("Financial Econometrics week 2 R code stock_returns.txt")

########################
############ 2. Set initialization for numerical optimization
########################

## 2a. choose the initail parameter values for the Newton Raphson optimization

a <- 0.2 # initial value for alpha
b <- 0.6  # initial value for beta
omega <- var(x)*(1-a-b) # initial value for omega

## 2b. Transform intitial values using the inverse of the link funtions

# In the likelihood function llik_fun_GARCH() the parameter inputs are 
# transformed using link functions to ensure the parameter restrictions.
# This is done to avoid numerical problems in the optimization.
# The restrictions are the following (see llik_fun_GARCH): 
# - exp() to ensure omega>0
# - logistic()=exp()/(exp()) for 0<alpha<1
# - logistic()=exp()/(exp()) for 0<beta<1
# Therefore, we transform the parameter values with the inverse of these functions
# to give the desired initialization:

par_ini <- c(log(omega),log(a/(1-a)),log(b/(1-b)))

########################
############ 3. Optimize Log Likelihood function
########################

#Note that the numerical optimizer optim() minimizes a function. We want
#to maximize the log-likelihood. Therefore we need to provide a negative
#log-likelihood. Note also that for numerical reasons we optimize the
#average log-likelihood instead of the log-likelihood.

# Optim input:
# (1) Initial parameter (par): par_ini
# (2) Negative average log-likelihood function (fn): - llik_fun_GARCH()
# (3) Numerical algorithm used (method): BFGS

# fmincon output:
# (1) estimates: $par
# (2) log likelihood function value at theta_hat: $value
# (3) exit flag indicating convergence: $convergence

est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x), method = "BFGS")

########################
############ 4. Obtain parameter estimates
########################

# est$par gives the values that maximize the likelihood
# However, these are not the paramter estimates since we have used the link functions
# to guarantee omega>0, 0<alpha<1 and 0<beta<1.
# Therefore, we apply these link functions to obtain the parameter estimates

omega_hat <- exp(est$par[1])
alpha_hat <- exp(est$par[2])/(1+exp(est$par[2]))
beta_hat <- exp(est$par[3])/(1+exp(est$par[3]))

theta_hat <- c(omega_hat,alpha_hat,beta_hat)
cat("The parameter estimates are:")
round(theta_hat,4)

# display the log-likelihood. 
# Note that we multiply the average log-likelihood by the sample size and
# -1 to get the log-likelihood.

cat("The log-likelihood value is:")
-est$value*length(x)

# display the exit flag to see if convergence of the algorithm has been attained

cat("Exit flag:")
est$convergence # zero indicates succesfull optimization

########################
############ 5. Calculate covariance matrix of the ML estimator
########################

# An estimate of the covariance matrix of the ML estimator can be obtained 
# simply by inverting the hessian matrix 
# The hessian matrix is the second derivative of the negative average log-likelihood
# with respect to the parameter values
# The function Hess_fun_GARCH() gives the average log-likelihood but without
# transforming the parameter input with the link functions.
# This is needed because we want the Hessian with respect to the original parameters

# We compute the numerically the Hessian matrix using optimHess().

hessian <- optimHess(par=theta_hat, fn=function(par)-Hess_fun_GARCH(par,x))

# We obtain the covariance matrix of the ML estimator by inverting the Hessian matrix
SIGMA <- solve(hessian)

cat("Asymptotic covariance matrix of the ML estimator:")
SIGMA

########################
############ 6. Confidence intervals for the parameters
########################

# the follwing code shows how to obtain a 95% confidence interval for beta 

lb_beta <- theta_hat[3]-1.96*sqrt(SIGMA[3,3])/sqrt(length(x)) # lower bound of the interval
ub_beta <- theta_hat[3]+1.96*sqrt(SIGMA[3,3])/sqrt(length(x)) # upper bound of the interval

ci_beta <- c(lb_beta, ub_beta)

cat("95% confidence interval for beta:")
ci_beta

########################
############ 7. Obtain the AIC and BIC selection criteria
########################

# first we obtain the sample size n, the number of parameters p and the
# log-likelihood value

n <- length(x)   
p <- 3  
llik_val <- -est$value*n

# the AIC and BIC are caluclated as

aic <- 2*p-2*llik_val
bic <- log(n)*p-2*llik_val

cat("The AIC is:")
aic

cat("The BIC is:")
bic