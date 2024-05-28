
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
# a1 = 0.5, a2 = -0.5, b = 1, sha = 1.5, l1 = 1/30, l2 = 1/50, r1 = 1.5, r2 = 1.5
# nn = 400

### Calculation of the upper control limit
ucl_sample <- function(a1, a0, b, sha, l1, l2, r1, r2, nn)
{
  ### Parameter Description:
  ### a1 and a0 are the regression coefficients and 
  ### intercept of the logit model, respectively;
  ### b is the regression parameter of the Cox PH model; 
  ### l2 and sha are the scale and shape parameters 
  ### of the Weibull distribution, respectively; 
  ### l1 is the parameter of the exponential 
  ### distribution to generate censoring times; 
  ### r1 and r2 are the values of RR and OR, respectively; 
  ### nn is the excepted ARL0.
  
  j <- 1; ps <- 0; ds <- 0
  while (j <= nn){
    x1 <- rnorm(1, 0, 1); x2 <- rnorm(1, 0, 1)
    p0 <- invlogit(x1*a1 + a0); e0 <- exp(x2*b)
    u0 <- runif(1, 0, 1); u1 <- runif(1, 0, 1)
    tc <- rexp(1, l1)
    t <- ifelse(u0 <= 1-p0, 1000000000, 
                (-log(u1)/(l2*e0))^(1/sha))
    if (t <= tc){
      delta <- 1
    } else {
      t <- tc; delta <- 0
    }
    ps <- ps + p0
    ds <- ds + delta
    j <- j+1
  }
  hh <- log(r2)*ps + (log(r1)-r1+1)*ds
  return(hh)
}

numCores <- detectCores()
registerDoParallel(numCores)
uss <- foreach(o = 1:10000, .combine = "c", 
              .packages = c("doParallel", "foreach", 
                            "MASS", "logitnorm")) %dopar% 
  {ucl_sample(a1, a0, b, sha, l1, l2, r1, r2, nn)
  }
stopImplicitCluster()
us <- mean(uss)


### General form of RACure-CUSUM
racure_sample <- function(a1, a0, b, sha, l1, l2, r1, r2, h)
{ 
  ### Parameter Description:
  ### a1 and a0 are the regression coefficients and 
  ### intercept of the logit model, respectively;
  ### b is the regression parameter of the Cox PH model; 
  ### l2 and sha are the scale and shape parameters 
  ### of the Weibull distribution, respectively; 
  ### l1 is the parameter of the exponential 
  ### distribution to generate censoring times; 
  ### r1 and r2 are the values of RR and OR, respectively; 
  ### h is the upper control limit.

  rl <- 1; z <- 0
  while (rl > 0){
    x1 <- rnorm(1, 0, 1); x2 <- rnorm(1, 0, 1)
    p0 <- invlogit(x1*a1 + a0); e0 <- exp(x2*b)
    u0 <- runif(1, 0, 1); u1 <- runif(1, 0, 1)
    tc <- rexp(1, l1)
    t <- ifelse(u0 <= 1-p0, 1000000000, 
                (-log(u1)/(l2*e0))^(1/sha))
    if (t <= tc){
      delta <- 1
    } else {
      t <- tc; delta <- 0
    }

    H0 <- l2*(t^sha)*e0; ss <- exp(-H0)
    if (ss < 1e-40) {ss <- 1e-40}
    yy <- delta + (1-delta) * p0*ss/(1-p0+p0*ss)

    lr <- yy*log(r2) + delta*yy*log(r1) + 
      yy*(r1-1)*log(ss)
    z <- max(0, z + lr)
    if (z >= h){
      return(rl)
      break
    } else {
      rl <- rl + 1
    }
  }
}

numCores <- detectCores()
registerDoParallel(numCores)
rls <- foreach(o = 1:10000, .combine = "c", 
               .packages = c("doParallel", "foreach", 
                             "MASS", "logitnorm")) %dopar% 
  {racure_sample(a1, a0, b, sha, l1, l2, r1, r2, h = us)
  }
stopImplicitCluster()
arl <- mean(rls)

### Validation
round(arl - nn, 3) # Should be close enough

