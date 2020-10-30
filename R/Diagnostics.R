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


set.seed(234)
#set.seed(1)
temp <- genSurvData(n = 200, p = 500, rho = 0.25, scale=2.0, shape=1.5, cens.quant = .6)
X <- temp$X
logY <- log(temp$time)
delta <- temp$status
betastar <- temp$beta
library(Rcpp)
source("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/R/penAFT.R")
source("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/R/ADMM.ENpath.R")
source("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/R/ADMM_ENrun.R")
source("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/R/ADMM.SGpath.R")
source("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/R/ADMM.SGrun.R")
sourceCpp("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/src/ADMM_ENrun.cpp")
sourceCpp("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/src/Test.cpp")

#sourceCpp("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/src/ADMM_ENrun_Eigen.cpp")

u <- 100
timing.vec_200 <- vector(length = u)
#for (i in 1:u) { 
ptm <- proc.time()
store_n_200 <- penAFT(
  X, logY, delta, 
    nlambda = 100,
    lambda.ratio.min = .1, lambda = NULL,
    penalty = "EN",
    alpha = 1, weights = NULL,
    tol.abs = 1e-10, 
    tol.rel = 5e-4,
    gamma = 2, center = TRUE,
    standardize = FALSE, 
    admm.max.iter = 1e4)
proc.time() - ptm

#}

# source("/Users/aaron/Documents/GitHub/PenalizedGehan/Code/Penalized Gehan/gehan_lasso.R")
#  eval.obj <- function(logY, XB, beta, delta, lambda, alpha, w){
#     out <- 0
#     n <- dim(X)[1]
#     E <- logY - XB
#     for(i in which(delta==1)){
#       for(j in 1:n){
#         out <- out + max(E[j] - E[i], 0)
#       }
#     }
#     return(out/n^2 + lambda*alpha*sum(w*abs(beta)) + lambda*0.5*(1-alpha)*sum(w*beta^2))
#   }


# for(j in 1:10){
#   keep <- gehan.lasso(x = X,y = logY,delta = delta,lambda = 2*store$lambda[(j-1)*10 + 1],CENTER=FALSE)
#   cat(eval.obj(logY, X%*%keep$b, keep$b, delta, store$lambda[(j-1)*10 + 1], alpha = 1, w = rep(1, dim(X)[2])), ";", eval.obj(logY, X%*%store$beta[,(j-1)*10 + 1], store$beta[,(j-1)*10 + 1], delta, store$lambda[(j-1)*10 + 1], alpha = 1, w = rep(1, dim(X)[2])), "\n")
# }



#store_1 - cpp
#store_2 - R

Beta_1 <- store_p_5000$beta
Beta_2 <- store_R$beta

obj.vec_1 <- vector(length = ncol(Beta_1))
obj.vec_2 <- vector(length = ncol(Beta_1))

alpha <- 1
p <- 500
w <- rep(1, p)

for (i in 1:ncol(Beta_1)) {
  obj.vec_1[i] <- eval.obj(logY, crossprod(t(X), Beta_1[,i]), Beta_1[,i], delta, lambda[i])
  obj.vec_2[i] <- eval.obj(logY, crossprod(t(X), Beta_2[,i]), Beta_2[,i], delta, lambda[i])
}

kk <- 6
Beta_1[,kk][abs(Beta_1[,kk]) > 1e-6]

Beta_2[,kk][Beta_2[,kk] != 0]

sum(abs(Beta_1 - Beta_2))

max(abs(Beta_1 - Beta_2))


#lambda <- store_1


Beta_175 <- store_175$beta

Beta_200 <- store_200$beta

length(Beta_175[Beta_175 != 0])
length(Beta_200[Beta_200 != 0])

lambda <- c(0.1318143,
0.102059,
0.07902061,
0.06118279,
0.04737161,
0.03667812,
0.02839854,
0.02198796,
0.01702448,
0.01318143)



###########################
###     SG       ##########
###########################

library(Rcpp)
source("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/R/penAFT.R")
source("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/R/ADMM_SGpath.R")
source("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/R/ADMM_SGrun.R")

source("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/R/ADMM.SGpath.R")
source("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/R/ADMM.SGrun.R")

sourceCpp("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/src/ADMM_SGrun.cpp")
sourceCpp("C:/Users/DELL/Desktop/Prof M penAFT/penAFT/src/Test.cpp")

set.seed(234)
#set.seed(1)
temp <- genSurvData(n = 175, p = 500, rho = 0.25, scale=2.0, shape=1.5, cens.quant = .6)
X <- temp$X
logY <- log(temp$time)
delta <- temp$status
betastar <- temp$beta

p <- 500
groups <- floor(runif(p, min = 1, max = 21))
groups <- sort(groups)
w <- runif(p, 1, 10)
v <- runif(length(unique(groups)), 1, 5)

weights <- data.frame(w, v)


#t0 <- Sys.time()
#for (i in 1:10) {
#out <- ADMM.ENpath(X.fit = X, logY = logY, delta = delta, nlambda = 1, max.iter = 20000, alpha = 0, w = w, tol.abs = 1e-20, tol.rel = 1e-4, gamma = 2)
#}
#t1 <- Sys.time()
#t1 - t0

#out_cpp <- ADMM.SGpath(X.fit = X, logY = logY, delta = delta, nlambda = 50, max.iter = 5000, alpha = 0.5, w = w, v = v, tol.abs = 1e-20, tol.rel = 1e-4, gamma = 2, groups = groups)

ptm <- proc.time()
store_R <- penAFT(
  X, logY, delta, 
  nlambda = 10,
  lambda.ratio.min = .1, lambda = NULL,
  penalty = "SG",
  alpha = 0.5, weights = weights, groups = groups,
  tol.abs = 1e-10, 
  tol.rel = 5e-4,
  gamma = 2, center = TRUE,
  standardize = FALSE, 
  admm.max.iter = 1e4)
proc.time() - ptm



Beta_1 <- store_cpp$beta
Beta_2 <- store_R$beta

sum(abs(Beta_1 - Beta_2))


Beta_1_1 <- Beta_1[,2]
Beta_2_1 <- Beta_2[,2]

sum(abs(Beta_1_1 - Beta_2_1))

Beta_1_1[Beta_1_1 != 0]
Beta_2_1[Beta_2_1 != 0]

