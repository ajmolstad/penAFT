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
  beta <- sample(c(0,1),p, replace = TRUE, prob=c(.9, .1))/5

  # ------------------------------------------------
  # Generate responses from Weibull AFT
  # ------------------------------------------------
  logY <- X%*%beta + rnorm(n, mean = 2, sd = 1)

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


set.seed(10)
#set.seed(1)
temp <- genSurvData(n = 100, p = 1000, rho = 0.0, scale=2.0, shape=1.5, cens.quant = .6)
X <- temp$X
logY <- log(temp$time)
#logY <- logY - rep(1, dim(logY)[1])%*%t(colMeans(logY))
X <- X - rep(1, dim(X)[1])%*%t(colMeans(X))
delta <- temp$status
betastar <- temp$beta



store <- penAFT(X, logY, delta, 
                   nlambda = 50, 
                   lambda.ratio.min = NULL, lambda = NULL, 
                   penalty = NULL,
                   alpha = 1, weights = NULL, 
                   groups = NULL, tol.abs = 1e-10, 
                   tol.rel = 5e-4, 
                   gamma = .25, center = TRUE, 
                   standardize = FALSE)







