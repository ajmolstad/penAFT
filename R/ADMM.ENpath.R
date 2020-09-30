ADMM.ENpath <- function(X.fit, logY, delta, lambda, alpha, w, tol.abs, tol.rel, gamma){
  
  	# -------------------------------------
	# Objective function evaluator 
	# -------------------------------------
	eval.obj <- function(logY, XB, beta, delta, lambda){
		out <- 0 
		n <- dim(X)[1]
		E <- logY - XB
		for(i in which(delta==1)){
			for(j in 1:n){
				out <- out + max(E[j] - E[i], 0)
			}
		}
		return(out/n^2 + lambda*alpha*w*sum(abs(beta)) + lambda*0.5*(1-alpha)*w*sum(beta^2))
	}

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
	if(n < 300){
		eta <- max(eigen(crossprod(crossprod(t(D), X)))$val)  
	} else {
		eta <- n*max(svd(X)$d)^2
	}
	Xbeta <- crossprod(t(X), Beta)
	tXB <-  crossprod(t(crossprod(t(D), X)), Beta)
	#eta <- eta/2
	rho <- 2
	BetaOut <- Matrix(0, nrow=p, ncol=length(lambda), sparse=TRUE)
	euc.tildelogY <- sqrt(sum(tildelogY^2))

	for(kk in 1:length(lambda)){
		out <- ADMM.ENrun(tildelogY, X, D, tildedelta, rho = rho, eta = eta, tau = 1.5, 
			lambda = lambda[kk], alpha = alpha, w = w, Gamma = Gamma, Beta = Beta, 
			Theta = Theta, 
			max.iter = 5000, tol.abs = tol.abs, tol.rel = tol.rel, gamma = gamma, euc.tildelogY = euc.tildelogY)
		BetaOut[,kk] <- out$Beta
		Beta <- out$Beta
		Gamma <- out$Gamma
		Theta <- out$Theta
		rho <- out$rho
		cat("Through tp ", kk, "\n")
	}


	result <- list("beta" = BetaOut, "lambda" = lambda)

}
