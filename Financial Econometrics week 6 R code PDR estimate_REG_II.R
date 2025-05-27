### Estimate dynamic regression by indirect inference
# The following code can be used to estimate the dynamic regression model
# For a more detailed explantion of the code see the Lecture Notes

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

source("sim_m_REG.R") 

########################
############ 1. Load data
########################

# Load the variables xt and yt that are contained in the files "xt.txt" and "yt.txt"

xt <- scan("xt.txt")
yt <- scan("yt.txt")

########################
############ 2. Obtain sample moments
########################

hb <- cov(yt,xt)/var(xt)
yr <- yt-hb*xt

xy <- yr*xt
acvfxy <- acf(xy, lag.max=15, type ="covariance", plot=F)$acf[-1]

sample_m <- c(var(yr),hb,acvfxy)

########################
############ 3. Generate errors for the simulations from the regression model
########################

n <- length(xt)
M <- 20
H <- M*n

set.seed(123)
eta <- rnorm(H)
eps <- rnorm(H)
e <- cbind(eta,eps)
        
########################
############ 4. Choose initial parameter values for the optimization
########################    

a1 <- 0.95 
a0 <- cov(xt,yt)/var(xt)*(1-a1)
s_eta <- 0.2
s_eps <- var(yr)
            
par_ini <- c(a0, log(a1/(1-a1)), log(s_eta), log(s_eps))

########################
############ 5. Obtain parameter estimates
########################

# 3a. Minimize criterion function

est <- optim(par=par_ini,fn=function(par) mean((sim_m_REG(e,xt,par)-sample_m)^2), method = "BFGS")

# 3a. Obtain parameter extimate using the link functions

a0_hat <- est$par[1]
a1_hat <- exp(est$par[2])/(1+exp(est$par[2]))
s_eta <- exp(est$par[3])
s_eps <- exp(est$par[4])

theta_hat <- c(a0_hat,a1_hat,s_eta,s_eps)
cat("The parameter estimates are:")
round(theta_hat,4)