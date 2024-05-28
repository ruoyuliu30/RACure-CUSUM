
library(logitnorm)
library(penPHcure)
library(survival)
library(smoothHR)

### DESCRIPTION:
### Here are provided:
### 1. Load the dataset WHAS500 from R-package smoothHR;
### 2. Construct the risk adjustment for RACure-CUSUM 
###    by R-package penPHcure; 
### 3. Calculation of the upper control limit;
### 4. RACure-CUSUM for Phase II data of WHAS500.

### Version information:
### R v4.0.1
### penPHcure v1.0.1
### Survival v3.2-13
### smoothHR v1.0.3

### Note: 
### Firstly, we preprocess the WHAS500 dataset, 
### removing irrelevant information. Then, we convert 
### the entry time into the study from the format 
### "xx-xx-xxxx" to precise days and calculate the 
### total survival time. All other data remains the 
### same as in the original dataset. The processed 
### data has also been uploaded.


### Main code

### Load data
# original: data("whas500"); whas <- whas500
whas <- read.csv("RACure-CUSUM/data/WHAS500_prepro.csv")

### Transform patient data based on the expected 
### stopping time (i.e., ARL0).
# Specifically, considering the time point of entry 
# into the study, truncate patient observations based 
# on the expected stopping time. For example, a patient 
# who entered the study on the 100th day with a survival 
# time of 150 days should have their data truncated to 
# reflect survival for 100 days when observation ceases 
# on the 200th day, with no observed death.
data_tran <- function(tt, data)
{
  ### Parameter Description:
  ### tt is ARL0; data is the original data 
  data$st <- data$gt <- data$ct <- data$ds <- rep(0, dim(data)[1])
  for (i in 1:dim(data)[1]){
    if (data$zt[i] > tt){
      data$ct[i] <- tt; data$ds[i] <- 0
    } else {
      data$ct[i] <- data$zt[i]; data$ds[i] <- data$fstat[i]
    }
    data$gt[i] <- data$ct[i] - data$entry[i]
  }
  return(data)
}

# In our case study, let tt = 730, 
whas0 <- data_tran(tt, data = whas)

### Risk adjustment 
# training data (1997-2000) for risk adjustment
whas1 <- whas0[which(whas0$year != 3),]
# Phase I data (1999-2000) for UCL
whas2 <- whas0[which(whas0$year == 2),]
# Phase II data (1999-2000) to monitor
whas3 <- whas0[which(whas0$year == 3),]

ccm <- penPHcure(Surv(time = st, time2 = gt, 
                      event = ds) 
                 ~ age+gender+hr+bmi+chf, 
                 cureform = ~ hr+bmi+sysbp+diasbp+chf, 
                 data = whas1, standardize = FALSE, ties = "breslow")

### Obtain the risk-adjusted resluts for quickly computing 
# Phase I data
ccr1 <- predict(ccm, whas2)
whas2$p <- ccr1$CURE; whas2$s <- ccr1$SURV
# Phase II data
ccr2 <- predict(ccm, whas3)
whas3$p <- ccr2$CURE

### Calculate the approximate UCL
racure_u <- function(r1, r2, data)
{              
  ### Parameter Description:
  ### r1 and r2 are the values of RR and OR, respectively; 
  ### data is the Phase I data
  nd <- sum(data$ds); np <- sum(data$p)
  hh <- log(r2)*np + (log(r1)-r1+1)*nd
  return(hh)
}
# In our case study, let r1 = 1.5, r2 = 1.5, 
au <- racure_u(r1, r2, data = whas2)
# To facilitate analyzing, the UCL is scaled down. 
# The subsequent calculation of chart statistics will 
# be scaled down by the same factor.
# This approach does not alter the nature of RACure-CUSUM; 
# it is merely a numerical adjustment.
au <- au / 10

### Calculate the chart statistics of RACure-CUSUM
racure_st <- function(r1, r2, arl, ram, data)
{              
  ### Parameter Description:
  ### r1 and r2 are the values of RR and OR, respectively; 
  ### arl is ARL0; 
  ### ram is the risk-adjusted model
  ### data is the Phase II data
  rl <- 0; lrs <- c(0); sts <- c()
  while (rl >= 0){
    if (length(which(data$entry <= rl)) == 0){
      lrs <- c(lrs, 0); sts <- c(sts, 0)
      rl <- rl + 1
      next
    } else {
      sd <- data[which(data$entry <= rl),]
      sd$zt2 <- sd$zt; sd$fstat2 <- sd$fstat
      sd$H <- rep(0, dim(sd)[1]); sd$yy <- rep(0, dim(sd)[1])
      sd$yd <- rep(0, dim(sd)[1]); sd$yh <- rep(0, dim(sd)[1])
      break
    }
  }
  
  rl <- rl + 1
  while (rl > 0){
    for (j in 1:dim(sd)[1]){
      if (sd$zt[j] > rl){
        sd$zt[j] <- rl; sd$fstat[j] <- 0
      } else {
        sd$zt[j] <- sd$zt2[j]; sd$fstat[j] <- sd$fstat2[j]
      }
      sd$gt[j] <- sd$zt[j] - sd$entry[j]
      
      pres <- predict(ram, sd[j,]); ss <- pres$SURV
      # if (ss < 1e-40) {ss <- 1e-40}
      sd$H[j] <- -log(ss)
      sd$yy[j] <- sd$fstat[j] + 
        (1-sd$fstat[j]) * sd$p[j]*ss/(1-sd$p[j]+sd$p[j]*ss)
      sd$yd[j] <- sd$yy[j] * sd$fstat[j]
      sd$yh[j] <- sd$yy[j] * sd$H[j]
    }
    
    lr <- log(r2)*sum(sd$yy) + 
      log(r1)*sum(sd$yd) - (r1-1)*sum(sd$yh)
    lr <- lr / 10; lrs <- c(lrs, lr)
    lrz <- lr - min(lrs); sts <- c(sts, lrz)
    if (rl == arl){
      return(sts)
      break
    } else {
      sd <- data[which(data$entry <= rl),]
      sd$zt2 <- sd$zt; sd$fstat2 <- sd$fstat
      sd$H <- rep(0, dim(sd)[1]); sd$yy <- rep(0, dim(sd)[1])
      sd$yd <- rep(0, dim(sd)[1]); sd$yh <- rep(0, dim(sd)[1])
      rl <- rl + 1
    }
  }
}
# In our case study, let r1 = 1.5, r2 = 1.5, arl = 730, 
cst <- racure_st(r1, r2, arl, ram = ccm, data = whas3)
# Find out the time point of the first alarm
fat <- sort(which(cst >= au))[1]









