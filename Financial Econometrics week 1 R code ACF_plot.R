########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(yfR)

########################
############ 1. Obtain log returns
########################

# 1a. As discussed in the R file "Load_Yahoo_Data.m", we download IBM prices and 
# obtain log-returns

ibm_data <- yf_get(tickers = "IBM", 
                   first_date = "2005-01-01",
                   last_date = "2018-01-01",freq_data = "daily")

ibm_prices <- ibm_data$price_adjusted
ibm_log_ret <- diff(log(ibm_prices))

# 1b. We obtain the squared log-returns simply by taking the square of the log-returns
# Note that R by default applies a function element by element 

ibm_log_ret2 <- ibm_log_ret^2;

########################
############ 2. Plot sample ACF and PACF of squared log-returns
########################

# We use the R functions "acf()" and "pacf()" to plot the ACF and PACF

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1)) 
acf(ibm_log_ret2,main="")
title("ACF squared log-returns IBM", line = 0.3)
pacf(ibm_log_ret2, main="")
title("PACF squared log-returns IBM", line = 0.3)

