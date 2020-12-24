ADMM.ENpath <- function(X.fit, logY, delta, admm.max.iter, lambda, alpha, w, tol.abs, tol.rel, gamma, quiet, cpp){
  
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
	D <- matrix(0, nrow = l, ncol = n)
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
				tildedelta[counter,] <- c(delta[j], delta[k])
				counter <- counter + 1
			}
		}
	}
	tildelogY <- Theta
	D <- Matrix(D, sparse=TRUE)
	Gamma  <- -sign(Theta)
	Beta <- rep(0, p)
	#ptm <- proc.time()
	t1 <- as.matrix(tcrossprod(D, t(X)))
	eta <-  irlba(t1,1)$d^2
	#proc.time() - ptm
	t1 <- NULL
	rho <- 0.1

	BetaOut <- Matrix(0, nrow=p, ncol=length(lambda), sparse=TRUE)
	euc.tildelogY <- sqrt(sum(tildelogY^2))

	for(kk in 1:length(lambda)){
		if(!cpp){
			out <- ADMM.ENrun(tildelogY, X, D, D.pos, tildedelta, rho = rho, eta = eta, tau = 1.5, 
				lambda = lambda[kk], alpha = alpha, w = w, Gamma = Gamma, Beta = Beta, 
				Theta = Theta, 
				max.iter = 5000, tol.abs = tol.abs, tol.rel = tol.rel, 
				gamma = gamma, euc.tildelogY = euc.tildelogY)
		} 
		if(cpp){
			out <- ADMM_ENrun(tildelogY, X, D, D.pos, tildedelta, rho = rho, eta = eta, tau = 1.5, 
				lambda = lambda[kk], alpha = alpha, w = w, Gamma = Gamma, Beta = Beta, 
				Theta = Theta, 
				max_iter = 5000, tol_abs = tol.abs, tol_rel = tol.rel, 
				gamma = gamma, euc_tildelogY = euc.tildelogY)
		}
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
