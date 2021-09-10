genSurvData <- function(n, p, s, mag, cens.quant = 0.6){

  # --------------------------------
  # Generate predictors + beta
  # --------------------------------
  SigmaX <- matrix(0, nrow = p, ncol = p)
  for (j in 1:p) {
    for (k in 1:p) {
      SigmaX[j,k] <- 0.7^(abs(j-k))
    }
  }
  eo <- eigen(SigmaX)
  SigmaXSqrt <- eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)

  X <- tcrossprod(matrix(rnorm(n*p), nrow=n, ncol=p),SigmaXSqrt)
  beta <- sample(c(rep(0, p-s), rep(mag, s)))*sample(c(-1,1), p, replace = TRUE)

  # ------------------------------------------------
  # Generate responses from log-logistic AFT
  # ------------------------------------------------
  logtime <- X%*%beta + rlogis(n, location = 0, scale = 2)

  # -------------------------------------
  # Generate censoring times
  # -------------------------------------
  temp <- quantile(exp(logtime), cens.quant)
  C <- rexp(n=n, rate=1/temp)
  logY <- pmin(logtime, log(C))
  status <- 1*(logtime == logY)

  return(list(
    "beta" = beta,
    "logY" = logY,
    "status" = status,
    "X" = X
  ))
}
