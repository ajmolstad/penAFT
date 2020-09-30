ADMM.ENrun <- function(tildelogY, X, D, tildedelta, rho, eta, tau, 
			lambda, alpha, w, Gamma, Beta, 
			Theta,  
			max.iter, tol.abs, tol.rel, gamma, euc.tildelogY){

	p <- dim(X)[2]
	n <- dim(X)[1]
	l <- length(tildelogY)
	updateStep <- 1
	Xbeta <- crossprod(t(X), Beta)
	tXB <-  crossprod(t(D),Xbeta)

	for(lll in 1:max.iter){

		# ---------------------------------
		# Theta update
		# ---------------------------------
		tTheta <- Theta
		Theta <- rep(0, l)
		nrho <- (n^(2 - gamma))*rho
		temp <- tildelogY - tXB - Gamma/rho
		inds1 <- which(temp > tildedelta[,2]/(nrho))
		inds2 <- which(temp < -tildedelta[,1]/(nrho))
		Theta[inds1] <- temp[inds1] - tildedelta[inds1,2]/(nrho)
		Theta[inds2] <- temp[inds2] + tildedelta[inds2,1]/(nrho)

		# -------------------------------------
		# Beta update
		# -------------------------------------
		#for(jj in 1:5){
			W <- (1/(eta))*crossprod(X, crossprod(D, tildelogY - Theta - Gamma/rho - tXB)) + Beta
			Beta <- Matrix(pmax(abs(W[,1]) - ((n^gamma)*lambda*alpha*w)/(eta*rho), 0)*sign(W[,1]), sparse=TRUE)/(1 + ((n^gamma)*lambda*(1-alpha)*w)/(eta*rho))
			tXB <- crossprod(t(D), crossprod(t(X), Beta))
		#}
		# ----------------------------------
		# Gamma
		# ----------------------------------
		Gamma <- Gamma + 1.5*rho*(Theta - tildelogY + tXB)

		if(lll%%floor(2*updateStep) == 0){

	        s <- rho*sqrt(sum(crossprod(X,crossprod(D, Theta - tTheta))^2))#/length(tildelogY)
			r <- sqrt(sum((Theta - tildelogY + tXB)^2))#/length(tildelogY)

			eprim <- sqrt(l)*tol.abs + tol.rel*max(c(sqrt(sum(tXB^2)), sqrt(sum(Theta^2)), euc.tildelogY))
			edual <- sqrt(p)*tol.abs + tol.rel*sqrt(sum(crossprod(X, crossprod(D, Gamma))^2))
		
			if(r/eprim > 10*s/edual){
				rho <- rho*2
			} 
			if(s/edual > 10*r/eprim){
				rho <- rho/2
			}
		
			cat("r = ", r, ": s= ", s, "\n")

			if(r < eprim & s < edual){
				break
			}
			updateStep <- updateStep + 1

		}

	}

	if(lll == max.iter){
		warning("ADMM did not converge in max.iter iterations", "\n")
	}

	return(list(
		"Beta" = Beta,
		"Gamma" = Gamma,
		"Theta" = Theta,
		"rho" = rho
		))
}