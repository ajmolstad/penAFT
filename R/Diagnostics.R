# -------------------------------------------------
library(MASS)
library(Matrix)
genSurvData <- function(n, p, rho, scale=1, shape=1.5, cens.quant = 0.65){

  # --------------------------------
  # Generate predictors + beta
  # --------------------------------
  SigmaX <- matrix(0, p, p)
  for(j in 1:p){
    for(k in 1:p){
      SigmaX[j,k] <- rho^(abs(j-k))
    }
  }

  X <- mvrnorm(n = n, mu = rep(0, p), SigmaX, tol = 1e-06, empirical = FALSE)
  beta <- sample(c(0,1),p, replace = TRUE, prob=c(.95, .05))

  # ------------------------------------------------
  # Generate responses from Weibull AFT
  # ------------------------------------------------
  logY <- X%*%beta + rnorm(n, mean = 2, sd = 4)

  # -------------------------------------
  # Generate censoring times
  # -------------------------------------
  temp <- quantile(logY, cens.quant)
  C <- rexp(n=n, rate=1/temp)
  # follow-up times and event indicators
  time <- pmin(exp(logY), exp(C))
  status <- 1*(logY <= C)

  return(list(
    "beta" = beta,
    "time" = time,
    "status" = status,
    "logY" = logY,
    "X" = X
  ))
}


set.seed(1)
#set.seed(1)
temp <- genSurvData(n = 100, p = 500, rho = 0.25, scale=2.0, shape=1.5, cens.quant = .6)
X <- temp$X
logY <- log(temp$time)
#logY <- logY - rep(1, dim(logY)[1])%*%t(colMeans(logY))
X <- X - rep(1, dim(X)[1])%*%t(colMeans(X))
delta <- temp$status
betastar <- temp$beta
library(Rcpp)
source("/Users/aaron/Documents/GitHub/penAFT/R/penAFT.R")
source("/Users/aaron/Documents/GitHub/penAFT/R/ADMM.ENpath.R")
source("/Users/aaron/Documents/GitHub/penAFT/R/ADMM.ENrun.R")


ptm <- proc.time()
store <- penAFT(X, logY, delta, 
    nlambda = 100, 
    lambda.ratio.min = .5, lambda = NULL, 
    penalty = NULL,
    alpha = 1, weights = NULL, 
    groups = NULL, tol.abs = 1e-10, 
    tol.rel = 1e-4, 
    gamma = 0, center = TRUE, 
    standardize = FALSE)
proc.time() - ptm




plot(store$beta[1,], type="l", ylim=c(-.1,1))
for(j in 2:dim(X)[2]){
  if(any(store$beta[j,]!=0)){
    lines(store$beta[j,], col="black")
  }
}


plot(store$beta[1,], type="l", ylim=c(-.5,.5))
for(j in 2:dim(X)[2]){
  if(any(store$beta[j,]!=0)){
    lines(store$beta[j,],col=j)
  }
}

keep <- gehan.lasso(x = X,y = logY,delta = delta,lambda = 2*store$lambda[75],CENTER=TRUE) 

