rm(list=ls())

library(irlba)

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
temp <- genSurvData(n = 200, p = 500, rho = 0.25, scale=2.0, shape=1.5, cens.quant = .6)
X <- temp$X
logY <- log(temp$time)
#logY <- logY - rep(1, dim(logY)[1])%*%t(colMeans(logY))
X <- X - rep(1, dim(X)[1])%*%t(colMeans(X))
delta <- temp$status
betastar <- temp$beta
library(Rcpp)
source("C:/Users/DELL/Desktop/To unpack 2/PenAFT_Loop/R/penAFT.R")
source("C:/Users/DELL/Desktop/To unpack 2/PenAFT_Loop/R/ADMM.ENpath.R")

source("C:/Users/DELL/Desktop/To unpack 2/PenAFT_Loop/R/ADMM.ENpath2.R")

source("C:/Users/DELL/Desktop/To unpack 2/PenAFT_Loop/R/ADMM.ENrun.R")

source("C:/Users/DELL/Desktop/To unpack 2/PenAFT_Loop/R/ADMM_ENrun.R")
source("C:/Users/DELL/Desktop/To unpack 2/PenAFT_Loop/R/ADMM.SGpath.R")
source("C:/Users/DELL/Desktop/To unpack 2/PenAFT_Loop/R/ADMM.SGrun.R")
sourceCpp("C:/Users/DELL/Desktop/To unpack 2/PenAFT_Loop/src/test.cpp")
sourceCpp("C:/Users/DELL/Desktop/To unpack 2/PenAFT_Loop/src/ADMM_ENrun.cpp")
sourceCpp("C:/Users/DELL/Desktop/To unpack 2/PenAFT_Loop/src/ADMM_ENrun2.cpp")

ptm <- proc.time()
store_cpp <- penAFT(
  X, logY, delta, 
    nlambda = 100,
    lambda.ratio.min = .5, lambda = NULL,
    penalty = "EN",
    alpha = .10, weights = NULL,
    groups = rep(1:10, each=dim(X)[2]/10), tol.abs = 1e-10, 
    tol.rel = 1e-5,
    gamma = 2, center = TRUE,
    standardize = FALSE, 
    admm.max.iter = 1e6,
    cpp = TRUE)
proc.time() - ptm

Beta1 <- store_cpp$beta

Beta2 <- store_cpp_2$beta

mean(abs(Beta1 - Beta2))

source("/Users/aaron/Documents/GitHub/PenalizedGehan/Code/Penalized Gehan/gehan_lasso.R")
 eval.obj <- function(logY, XB, beta, delta, lambda, alpha, w){
    out <- 0
    n <- dim(X)[1]
    E <- logY - XB
    for(i in which(delta==1)){
      for(j in 1:n){
        out <- out + max(E[j] - E[i], 0)
      }
    }
    return(out/n^2 + lambda*alpha*sum(w*abs(beta)) + lambda*0.5*(1-alpha)*sum(w*beta^2))
  }


for(j in 1:10){
  keep <- gehan.lasso(x = X,y = logY,delta = delta,lambda = 2*store$lambda[(j-1)*10 + 1],CENTER=FALSE)
  cat(eval.obj(logY, X%*%keep$b, keep$b, delta, store$lambda[(j-1)*10 + 1], alpha = 1, w = rep(1, dim(X)[2])), ";", eval.obj(logY, X%*%store$beta[,(j-1)*10 + 1], store$beta[,(j-1)*10 + 1], delta, store$lambda[(j-1)*10 + 1], alpha = 1, w = rep(1, dim(X)[2])), "\n")
}
