# Download IBM stock data from YAHOO using the R function "yf_get".
# The function "yf_get" is part of the package "yfR"

# We need to set the following arguments to the function:
# tickers: IBM
# first_date: 2005-01-01 
# last_date: 1-Jan-2019
# freq_data: "daily"
# other frequencies can be chosen: "daily", "weekly", "monthly".

# Visist yahoo finance.com website to find out symbols of other stocks
# For instance, the symbol of the NASDAQ stock index is ^IXIC

# Keep in mind that:
# 1) the format of the dates is yy-mm-dd
# 2) "yfR" returns an R object that contains several
# information about IBM. We are interested in adjusted closing prices.

# The following lines of code downloads data on IBM stock from yahoo 
# and stores it in the R object "ibm_data"

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(yfR) #you may need to install the package first

########################
############ 1. Obtain log returns
########################

ibm_data <- yf_get(tickers = "IBM", 
                   first_date = "2005-01-01",
                   last_date = "2018-01-01",freq_data = "daily")

head(ibm_data)

# We now need to extract the closing prices from the object "ibm_data". 
# We store the closing prices in the R vector "ibm_prices".
# This can be done by selecting $price_adjusted.

ibm_prices <- ibm_data$price_adjusted

# We obtain log-returns taking the log of the IBM prices and then taking
# first differences with the command diff()

ibm_log_ret <- diff(log(ibm_prices))

########################
############ 2. Plot prices and log-returns
########################

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1)) 
# mfrow=c(2,1) sets 2 plots in 2 rows and 1 column 

plot(ibm_prices, type = "l", main="Prices IBM", ylab="", xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(ibm_log_ret, type = "l", main="Log-returns IBM", ylab="", xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
# "type" sets the type of line and "main" selects the title

