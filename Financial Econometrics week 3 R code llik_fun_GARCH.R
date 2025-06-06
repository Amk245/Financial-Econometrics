# This file contains the average log-likelihood function. 
# The function below takes the data, labeled "x", and the parameter vector, 
# labeled "par", as input and gives the average log-likelihood, labeled "llik", 
# as output. 


llik_fun_GARCH <- function(par,x){
  
  n <- length(x)
  
  #set paramter values from the input par using link functions for restrictions
  omega <- exp(par[1])                   #exp() to ensure omega>0
  alpha <- exp(par[2])/(1+exp(par[2]))   #logistic()=exp()/(exp()) for 0<alpha<1
  beta <- exp(par[3])/(1+exp(par[3]))    #logistic()=exp()/(exp()) for 0<beta<1
  
  ## Filter Volatility
  sig2 <- rep(0,n)
  sig2[1] <- var(x) #initialize volatility at unconditional variance
  
  for(t in 2:n){
    sig2[t] <- omega + alpha*x[t-1]^2 + beta*sig2[t-1]
  }
  
  ## Calculate Log Likelihood Values
  
  #construct sequence of log-likelihood contributions
  l <- -(1/2)*log(2*pi) - (1/2)*log(sig2) - (1/2)*x^2/sig2
  
  llik <- mean(l)  # obtain the average log-likelihood
  
  return(llik) # return the average log-likelihood as output
}