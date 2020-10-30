# ----------------------------------------------------------------------------
# Function for computing the solution path using sparse group lasso penalty
# -----------------------------------------------------------------------------
ADMM.SGpath <- function(X.fit, logY, delta, max.iter = 5000, lambda, alpha, w, v, groups, tol.abs, tol.rel, gamma, quiet) {


  # -------------------------------------
  # Objective function evaluator
  # -------------------------------------

  eval.obj <- function(logY, XB, beta, delta, lambda, groups){
    out <- 0
    n <- dim(X)[1]
    E <- logY - XB
    for(i in which(delta==1)){
      for(j in 1:n){
        out <- out + max(E[j] - E[i], 0)
      }
    }
    pen <- 0
    for(g in 1:length(unique(groups))){
      pen <- pen + sqrt(sum(beta[which(groups == g)]^2))
    }
    return(out/n^2 + lambda*pen)
  }

  X <- X.fit; X.fit <- NULL

  # ---------------------------------
  # Preliminaries
  # ---------------------------------
  n <- length(logY)
  p <- dim(X)[2]



  gradient_g <- function(X, logY, Beta, delta) {
    n <- nrow(X)
    p <- ncol(X)
    grad <- rep(0, p)
    Xbeta <- X%*%Beta
    for (i in 1:n) {
      for (j in 1:n) {
        e_i <- logY[i] - Xbeta[i]
        e_j <- logY[j] - Xbeta[j]
        Indicator <- (e_i <= e_j)
        term <- delta[i]*(X[i,] - X[j,]) * Indicator
        grad <- grad + term
      }
    }
    grad/n^2
  }



  # --------------------------------------------------------------------------------
  # Sorting indexes by their groups and determining "border indexes" of groups
  # --------------------------------------------------------------------------------
  indexes.data <- data.frame(groups, w, 1:p)
  names(indexes.data) <- c("group", "w", "index")
  indexes.data <- indexes.data[order(indexes.data$group),]
  index.vec <- indexes.data$index
  groups <- indexes.data$group
  w <- indexes.data$w
  border.indexes <- vector(length = groups[length(groups)] + 1)
  current.group <- 0
  counter.indexes <- 1
  for (i in 1:length(groups)) {
    if (current.group != groups[i]) {
      border.indexes[counter.indexes] <- i
      counter.indexes <- counter.indexes + 1
      current.group <- groups[i]
    }
  }
  border.indexes[length(border.indexes)] <- p + 1

  X <- X[, index.vec]
  

  #inv.data <- data.frame(index.vec, 1:p)
  #inv.data <- inv.data[order(index.vec),]

  #inv.index.vec <- inv.data[,2]

  # --------------------------------
  # Get initial values
  # --------------------------------
  l <- 1
  for(j in 1:(n-1)){
    for(k in (j+1):n){
      if(delta[j]!=0 | delta[k]!=0){
        l <- l + 1
      }
    }
  }
  l <- l - 1

  Theta <- rep(0, l)
  counter <- 1
  for(j in 1:(n-1)){
    for(k in (j+1):n){
      if(delta[j]!=0 | delta[k]!=0){
        Theta[counter] <- logY[j] - logY[k]
        counter <- counter + 1
      }
    }
  }

  Gamma  <- -sign(Theta)
  Beta <- rep(0, p)
  D <- matrix(0, nrow=l, ncol=n)
  counter <- 1
  for(j in 1:(n-1)){
    for(k in (j+1):n){
      if(delta[j]!=0 | delta[k]!=0){
        D[counter, j] <- 1
        D[counter, k] <- -1
        counter <- counter + 1
      }
    }
  }

  tildelogY <- rep(0, l)
  counter <- 1
  for(j in 1:(n-1)){
    for(k in (j+1):n){
      if(delta[j]!=0 | delta[k]!=0){
        tildelogY[counter] <- logY[j] - logY[k]
        counter <- counter + 1
      }
    }
  }


  #tildedelta <- matrix(0, nrow = l, ncol = 2)
  #counter <- 1
  #for(j in 1:(n-1)){
  #  for(k in (j+1):n){
  #    if(delta[j]!=0 | delta[k]!=0){
  #      tildedelta[counter,] <- c(delta[j], delta[k])
  #      counter <- counter + 1
  #    }
  #  }
  #}
  
  Theta <- rep(0, l)
  D <- matrix(0, nrow=l, ncol=n)
  tildelogY <- rep(0, l)
  tildedelta <- matrix(0, nrow = l, ncol = 2)
  counter <- 1
  for(j in 1:(n-1)){
    for(k in (j+1):n){
      if(delta[j]!=0 | delta[k]!=0){
        Theta[counter] <- logY[j] - logY[k]
        D[counter, j] <- 1
        D[counter, k] <- -1
        tildelogY[counter] <- logY[j] - logY[k]
        tildedelta[counter,] <- c(delta[j], delta[k])
        counter <- counter + 1
      }
    }
  }

  if(n < 200){
    eta <- max(eigen(crossprod(crossprod(t(D), X)))$val)  
  } else {
    eta <- n*max(svd(X)$d)^2
  }
  
  
  
  Xbeta <- crossprod(t(X), Beta)
  tXB <-  crossprod(t(crossprod(t(D), X)), Beta)
  #eta <- eta/2
  rho <- 1
  BetaOut <- Matrix(0, nrow=p, ncol=length(lambda), sparse=TRUE)
  euc.tildelogY <- sqrt(sum(tildelogY^2))


  D <- Matrix(D, sparse=TRUE)
  
  for(kk in 1:length(lambda)){
    out <- ADMM.SGrun(tildelogY, X, D, tildedelta, rho = rho, eta = eta, tau = 1.5,
                      lambda = lambda[kk], alpha = alpha, w = w, v = v, border.indexes = border.indexes, Gamma = Gamma, Beta = Beta,
                      Theta = Theta,
                      max.iter = max.iter, tol.abs = tol.abs, tol.rel = tol.rel, gamma = gamma, euc.tildelogY = euc.tildelogY, max.iter.update = 2000, G = groups[length(groups)])

    Beta.data <- data.frame(out$Beta, index.vec)
    names(Beta.data) <- c("Beta", "indeces")
    Beta.data.unsorted <- Beta.data[order(Beta.data$indeces),]
    BetaOut[,kk] <- Beta.data.unsorted$Beta
    Beta <- out$Beta
    Gamma <- out$Gamma
    Theta <- out$Theta
    rho <- out$rho
    iter.counter <- out$iter.counter
    
    
    tildedelta_nrho <- out$tildedelta_nrho

    if (kk == 7) {
      BetaHist <- out$BetaHist

    }

    if(iter.counter == max.iter){
      warning("ADMM did not converge in max.iter iterations", "\n")
    }
    cat("Through tp ", kk, "\n")

  }


  result <- list("beta" = BetaOut, "lambda" = lambda, "tildedelta_nrho" = tildedelta_nrho)

}




# -------------------------------------------------
#ibrary(MASS)
#library(Matrix)
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
