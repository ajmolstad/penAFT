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
		out <- ADMM.SGrun(tildelogY, X, D, tildedelta, rho = rho, eta = eta, tau = 1.5, 
			lambda = lambda[kk], alpha = alpha, w = w, v = v, groups = groups, Gamma = Gamma, Beta = Beta, 
			Theta = Theta, 
			max.iter = max.iter, tol.abs = tol.abs, tol.rel = tol.rel, gamma = gamma, euc.tildelogY = euc.tildelogY)
		BetaOut[,kk] <- out$Beta
		Beta <- out$Beta
		Gamma <- out$Gamma
		Theta <- out$Theta
		rho <- out$rho
		
		
		tildedelta_nrho <- out$tildedelta_nrho
		
		
		if (!quiet) {
			cat("Through ", kk,"th tuning parameter...", "\n")
		}
	}
  	

	result <- list("beta" = BetaOut, "lambda" = lambda, "tildedelta_nrho" = tildedelta_nrho)
  
}













