library(Rcpp)
source("C:/Users/DELL/Desktop/penAFT/R/ADMM_SGpath.R")
source("C:/Users/DELL/Desktop/penAFT/R/ADMM_SGrun.R")
sourceCpp("C:/Users/DELL/Desktop/penAFT/src/ADMM_SGrun.cpp")

#-----------------------------------------------

eval.obj <- function(logY, X, XB, beta, delta, lambda, alpha, w){
  out <- 0
  n <- dim(X)[1]
  E <- logY - XB
  for(i in which(delta==1)){
    for(j in 1:n){
      out <- out + max(E[j] - E[i], 0)
    }
  }
  penalty <- lambda * sum(w * (alpha * abs(beta) + (1/2) * (1 - alpha) * beta^2))
  return(out/n^2 + penalty)
}



# -------------------------------------------------
library(MASS)
library(Matrix)

set.seed(234)
#set.seed(1)
p <- 500
temp <- genSurvData(n = 175, p = p, rho = 0.0, scale=2.0, shape=1.5, cens.quant = .6)
X <- temp$X
logY <- log(temp$time)

X <- X - rep(1, dim(X)[1])%*%t(colMeans(X))
delta <- temp$status
betastar <- temp$beta
lambda.ratio <- .5
nlambda <- 100
#groups <- rep(1:20, each = dim(X)[2]/20)
#w <- rep(c(1,2,3), (p/3 + 3))[1:p]
#v <- rep(c(2,5), (groups[length(groups)] / 2) + 2)[1:groups[length(groups)]]


groups <- floor(runif(p, min = 1, max = 21))
w <- runif(p, 1, 10)
v <- runif(length(unique(groups)), 1, 5)




#t0 <- Sys.time()
#for (i in 1:10) {
#out <- ADMM.ENpath(X.fit = X, logY = logY, delta = delta, nlambda = 1, max.iter = 20000, alpha = 0, w = w, tol.abs = 1e-20, tol.rel = 1e-4, gamma = 2)
#}
#t1 <- Sys.time()
#t1 - t0

ptm <- proc.time()
out <- ADMM.SGpath(X.fit = X, logY = logY, delta = delta, nlambda = 50, max.iter = 5000, alpha = 0.5, w = w, v = v, tol.abs = 1e-20, tol.rel = 1e-4, gamma = 2, groups = groups)
proc.time() - ptm




Beta_out <- out$beta[,30]
Beta_out[Beta_out != 0]


Beta_out_mask <- abs(Beta_out) + groups

k <- 1
Beta_out_mask[Beta_out_mask > k - 1e-7 & Beta_out_mask < k + 1 - 1e-7]




#length(Beta_out[Beta_out != 0])

#eval.obj(logY, X, crossprod(t(X), Beta_out), Beta_out, delta, lambda = lambda, alpha = 1, w)

#lambda <- 0.1
#ft <- gehan.lasso(X, logY, delta, 2*lambda, CENTER=FALSE)

#eval.obj(logY, X, crossprod(t(X), ft$b), ft$b, delta, lambda = lambda, alpha = 1, w)


#ft$b[abs(ft$b) > 1e-7]

