llik_fun_sDVECH <- function(par,x){
  
  w11 <- exp(par[1])
  w12 <- par[2]
  w22 <- exp(par[3])
  a <- exp(par[4])/(1+exp(par[4]))
  b <- exp(par[5])/(1+exp(par[5]))
  
  d <- dim(x)
  n <- d[1]
  
  VECHt <- matrix(0,nrow=n,ncol=3)
  llik <- 0
  
  C <- cov(x)
  VECHt[1,] <- c(C[1,1],C[1,2],C[2,2])
  
  for(t in 2:n){
    
    VECHt[t,1] <- w11+b*VECHt[t-1,1]+a*x[t-1,1]^2
    VECHt[t,3] <- w22+b*VECHt[t-1,3]+a*x[t-1,2]^2
    VECHt[t,2] <- w12+b*VECHt[t-1,2]+a*x[t-1,1]*x[t-1,2]
    
    SIGMAt <- cbind(c(VECHt[t,1],VECHt[t,2]),c(VECHt[t,2],VECHt[t,3]))
    llik <- llik-0.5*(log(det(SIGMAt))+x[t,]%*%solve(SIGMAt)%*%t(t(x[t,])))/n
  }
  
  return(llik)
}