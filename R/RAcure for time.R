
library(logitnorm)
library(MASS)
library(foreach)
library(doParallel)


### DESCRIPTION:
### Here are provided:
### 1. Calculation of the upper control limit;
### 2. General form of RACure-CUSUM for 
###    Monte Carlo validation;
### 3. Code for parallel computation.

### Version information:
### R v4.0.1
### penPHcure v1.0.1
### Survival v3.2-13
### smoothHR v1.0.3

### Note: 
### We utilize parallel computing to accelerate 
### processing speed. Presented here is the general 
### form of the program for a single run of RACure-CUSUM.


### Main code
# Specific values for each parameter need to be 
# provided before use.
# For example in Scenario (1):
# a1 = 1, a2 = 0.5, b = 1, l1 = 1/30, l2 = 1/50, tr = 5
# r1 = 1.5, r2 = 1.5
# nt = 400

### The generation of data at each time point
data_generation <- function(a1, a0, b, sha, l1, l2, tr)
{
  ### Parameter Description:
  ### a1 and a0 are the regression coefficients and 
  ### intercept of the logit model, respectively;
  ### b is the regression parameter of the Cox PH model; 
  ### l2 is the parameter of the exponential 
  ### distribution to generate survival times; 
  ### l1 is the parameter of the exponential 
  ### distribution to generate censoring times; 
  ### tr is the number of new individuals entering 
  ### the study at each time point.
  
  x1 <- rnorm(tr, 0, 1); x2 <- rnorm(tr, 0, 1)
  p0 <- invlogit(x1*a1 + a0); e0 <- exp(x2*b)
  u0 <- runif(tr, 0, 1); u1 <- runif(tr, 0, 1)
  tc <- rexp(tr, l1); t <- rep(0, tr)
  ds <- rep(0, tr)
  for (i in 1:tr){
    if (u0[i] <= 1-p0[i]){
      t[i] <- 1000000000
    } else {
      t[i] <- -log(u1[i])/(l2*e0[i])
    }
    if (t[i] <= tc[i]){
      ds[i] <- 1
    } else {
      t[i] <- tc[i]
    }
  }
  sd0 <- data.frame(cbind(x1, x2, p0, e0, t, ds))
  return(sd0)
}

### Calculation of the upper control limit
ucl_time <- function(a1, a0, b, sha, l1, l2, tr, r1, r2, nt)
{ 
  ### Parameter Description:
  ### a1 and a0 are the regression coefficients and 
  ### intercept of the logit model, respectively;
  ### b is the regression parameter of the Cox PH model; 
  ### l2 is the parameter of the exponential 
  ### distribution to generate survival times; 
  ### l1 is the parameter of the exponential 
  ### distribution to generate censoring times; 
  ### tr is the number of new individuals entering 
  ### the study at each time point;
  ### r1 and r2 are the values of RR and OR, respectively; 
  ### nt is the excepted ARL0.
  tp <- nt * tr
  sd <- data_generation(a1 = a1, a0 = a0, b = b, 
                        l1 = l1, l2 = l2, tr = tp)
  ryt <- rep(c(0:(nt-1)), tr)
  sd$ht <- ryt; sd$zt <- sd$t + sd$ht
  sd$y1 <- rep(0, tp); sd$d1 <- rep(0, tp)
  for (j in 1:tp){
    if (sd$zt[j] > nt){
      sd$y1[j] <- nt - sd$ht[j]
      sd$d1[j] <- 0
    } else {
      sd$y1[j] <- sd$t[j]
      sd$d1[j] <- sd$ds[j]
    }
  }
  nd <- sum(sd$d1); np <- sum(sd$p0)
  hh <- log(r2)*np + (log(r1)-r1+1)*nd
  return(hh)
}

numCores <- detectCores()
registerDoParallel(numCores)
uts <- foreach(o = 1:10000, .combine = "c", 
               .packages = c("doParallel", "foreach", 
                             "MASS", "logitnorm")) %dopar% 
  {ucl_time(a1, a0, b, sha, l1, l2, tr, r1, r2, nt)
  }
stopImplicitCluster()
ut <- mean(uts)

### General form of RACure-CUSUM
racure_time <- function(a1, a0, b, sha, l1, l2, tr, r1, r2, h)
{ 
  ### Parameter Description:
  ### a1 and a0 are the regression coefficients and 
  ### intercept of the logit model, respectively;
  ### b is the regression parameter of the Cox PH model; 
  ### l2 is the parameter of the exponential 
  ### distribution to generate survival times; 
  ### l1 is the parameter of the exponential 
  ### distribution to generate censoring times; 
  ### tr is the number of new individuals entering 
  ### the study at each time point;
  ### r1 and r2 are the values of RR and OR, respectively; 
  ### h is the upper control limit.
  
  rl <- 0; lrs <- c()
  sd <- data_generation(a1 = a1, a0 = a0, b = b, 
                        l1 = l1, l2 = l2, tr = tr)
  sd$ht <- rep(rl, tr); sd$zt <- sd$t + sd$ht
  sd$nt <- rep(0, tr); sd$nd <- rep(0, tr)
  sd$H <- rep(0, tr); sd$yy <- rep(0, tr)
  sd$yd <- rep(0, tr); sd$yh <- rep(0, tr)
  dels <- c(0)
  
  rl <- 1
  while (rl > 0){
    tr1 <- tr*rl
    for (j in 1:tr1){
      if (j %in% dels){next}
      if (sd$zt[j] > rl){
        sd$nt[j] <- rl - sd$ht[j]
        sd$nd[j] <- 0
      } else {
        sd$nt[j] <- sd$t[j]
        sd$nd[j] <- sd$ds[j]
        dels <- c(dels, j); dels <- sort(unique(dels))
      }
      sd$H[j] <- l2*sd$nt[j]*sd$e0[j]
      ss <- exp(-sd$H[j])
      # if (ss < 1e-40) {ss <- 1e-40}
      sd$yy[j] <- sd$nd[j] + (1-sd$nd[j]) *
        sd$p0[j]*ss / (1-sd$p0[j]+sd$p0[j]*ss)
      sd$yd[j] <- sd$yy[j] * sd$nd[j]
      sd$yh[j] <- sd$yy[j] * sd$H[j]
    }
    
    lr <- log(r2) * sum(sd$yy) + 
      log(r1) * sum(sd$yd) - 
      (r1-1) * sum(sd$yh)
    lrs <- c(lrs, lr)
    lrz <- lr - min(lrs)
    if (lrz >= h){
      return(rl)
      break
    } else {
      sd0 <- data_generation(a1 = a1, a0 = a0, b = b, 
                             l1 = l1, l2 = l2, tr = tr)
      sd0$ht <- rep(rl, tr); sd0$zt <- sd0$t + sd0$ht
      sd0$nt <- rep(0, tr); sd0$nd <- rep(0, tr)
      sd0$H <- rep(0, tr); sd0$yy <- rep(0, tr)
      sd0$yd <- rep(0, tr); sd0$yh <- rep(0, tr)
      sd <- rbind(sd, sd0)
      rl <- rl + 1
    }
  }
}

numCores <- detectCores()
registerDoParallel(numCores)
rls <- foreach(o = 1:10, .combine = "c", 
               .packages = c("doParallel", "foreach", 
                             "MASS", "logitnorm")) %dopar% 
  {racure_time(a1, a0, b, sha, l1, l2, tr, r1, r2, h = ut)
  }
stopImplicitCluster()
arl <- mean(rls)

### Validation
round(arl - nt, 3) # Should be close enough


