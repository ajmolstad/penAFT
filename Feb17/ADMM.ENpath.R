ADMM.ENpath <- function(X.fit, logY, delta, admm.max.iter, lambda, alpha, w, tol.abs, tol.rel, gamma, quiet){
  
  	# -------------------------------------
	# Objective function evaluator 
	# -------------------------------------
	eval.obj <- function(logY, XB, beta, delta, lambda){
		out <- 0 
		n <- length(logY)
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
	D.pos <- matrix(0, nrow = l, ncol = 2)
	tildedelta <- matrix(0, nrow = l, ncol = 2)
	counter <- 1
	for(j in 1:(n-1)){
		for(k in (j+1):n){
			if(delta[j]!=0 | delta[k]!=0){
				Theta[counter] <- logY[j] - logY[k]
				D.pos[counter, 1] <- j
				D.pos[counter, 2] <- k
				tildedelta[counter,] <- c(delta[j], delta[k])
				counter <- counter + 1
			}
		}
	}
	tildelogY <- Theta
	D.vert.1 <- matrix(0, nrow = n, ncol = n)
	D.vert.neg1 <- matrix(0, nrow = n, ncol = n)
	counter.1 <- rep(1, times = n)
	counter.neg1 <- rep(1, times = n)
	for (i in 1:l) {
	  dPos.1 <- D.pos[i,1]
	  dPos.2 <- D.pos[i,2]
	  D.vert.1[dPos.1, counter.1[dPos.1]] <- i
	  D.vert.neg1[dPos.2, counter.neg1[dPos.2]] <- i
	  counter.1[dPos.1] <- counter.1[dPos.1] + 1
	  counter.neg1[dPos.2] <- counter.neg1[dPos.2] + 1
	} 
	
	Gamma  <- -sign(Theta)
	Beta <- rep(0, p)
	eta <-  n*(irlba(X,1)$d^2)
	rho <- 0.1

	BetaOut <- Matrix(0, nrow=p, ncol=length(lambda), sparse=TRUE)
	euc.tildelogY <- sqrt(sum(tildelogY^2))

	for(kk in 1:length(lambda)){
		out <- ADMM_ENrun(tildelogY, X, D.pos, D.vert.1, D.vert.neg1, tildedelta, rho = rho, eta = eta, tau = 1.5, 
			lambda = lambda[kk], alpha = alpha, w = w, Gamma = Gamma, Beta = Beta, 
			Theta = Theta, 
			max_iter = 5000, tol_abs = tol.abs, tol_rel = tol.rel, 
			gamma = gamma, euc_tildelogY = euc.tildelogY)
		BetaOut[,kk] <- out$Beta
		Beta <- out$Beta
		Gamma <- out$Gamma
		Theta <- out$Theta
		rho <- out$rho
		if (!quiet) {
			cat("Through ", kk,"th tuning parameter...", "\n")
		}
	}


	result <- list("beta" = BetaOut, "lambda" = lambda)

}
