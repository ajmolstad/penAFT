# ----------------------------------------------------------------------------
# Function for computing the solution path using sparse group lasso penalty
# -----------------------------------------------------------------------------
ADMM.SGpath <- function(X.fit, logY, delta, max.iter, lambda, alpha, w, v, groups, tol.abs, tol.rel, gamma, quiet) {
  
  X <- X.fit; X.fit <- NULL
  
  # ---------------------------------
  # Preliminaries
  # ---------------------------------
  n <- length(logY)
  p <- dim(X)[2]
  
  
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
  D <- matrix(0, nrow=l, ncol=n)
  
  D.pos <- matrix(0, nrow = l, ncol = 2)
  tildelogY <- rep(0, l)
  tildedelta <- matrix(0, nrow = l, ncol = 2)
  counter <- 1
  for(j in 1:(n-1)){
    for(k in (j+1):n){
      if(delta[j]!=0 | delta[k]!=0){
        Theta[counter] <- logY[j] - logY[k]
        D[counter, j] <- 1
        D[counter, k] <- -1
        
        D.pos[counter, 1] <- j
        D.pos[counter, 2] <- k
        
        tildelogY[counter] <- logY[j] - logY[k]
        tildedelta[counter,] <- c(delta[j], delta[k])
        counter <- counter + 1
      }
    }
  }
  
  D <- Matrix(D, sparse=TRUE)
  Gamma  <- -sign(Theta)
  Beta <- rep(0, p)
  if(n < 200){
    eta <- max(eigen(crossprod(crossprod(t(D), X)))$val)  
  } else {
    eta <- n*max(svd(X)$d)^2
  }
  Xbeta <- crossprod(t(X), Beta)
  tXB <-  crossprod(t(crossprod(t(D), X)), Beta)
  rho <- 1.5
  BetaOut <- Matrix(0, nrow=p, ncol=length(lambda), sparse=TRUE)
  euc.tildelogY <- sqrt(sum(tildelogY^2))
  
  for(kk in 1:length(lambda)){
    out <- ADMM.SGrun(tildelogY, X, D, D.pos, tildedelta, rho = rho, eta = eta, tau = 1.5, 
                      lambda = lambda[kk], alpha = alpha, w = w, v = v, border.indexes = border.indexes, Gamma = Gamma, Beta = Beta, 
                      Theta = Theta, 
                      max.iter = max.iter, tol.abs = tol.abs, tol.rel = tol.rel, gamma = gamma, euc.tildelogY = euc.tildelogY, G = groups[length(groups)])

    
    Beta.data <- data.frame(out$Beta, index.vec)
    names(Beta.data) <- c("Beta", "indeces")
    Beta.data.unsorted <- Beta.data[order(Beta.data$indeces),]
    BetaOut[,kk] <- Beta.data.unsorted$Beta

    Beta <- out$Beta
    Gamma <- out$Gamma
    Theta <- out$Theta
    rho <- out$rho
    
    
    #tildedelta_nrho <- out$tildedelta_nrho
    
    
    if (!quiet) {
      cat("Through ", kk,"th tuning parameter...", "\n")
    }
  }
  
  #BetaOutD <- Matrix(BetaOut, sparse = FALSE)
  #Beta.data <- data.frame(BetaOutD, index.vec)
  #Beta.data.unsorted <- Beta.data[order(Beta.data[,ncol(Beta.data)]),]
  #BetaOut <- Beta.data.unsorted[,1:(ncol(Beta.data)-1)]
  
  
  
  result <- list("beta" = BetaOut, "lambda" = lambda)
  
}



