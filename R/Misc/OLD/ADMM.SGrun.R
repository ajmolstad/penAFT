ADMM.SGrun <- function(tildelogY, X, D, tildedelta, rho, eta, tau, 
			lambda, alpha, w, v, groups, Gamma, Beta, 
			Theta,  
			max.iter, tol.abs, tol.rel, gamma, euc.tildelogY){

  #print(eta)
  #print("-----------------------------")
  #print(rho)
  #print("------------------------------------")
  #print(Theta[1:20])
  #print("------------------------------------")
  #print(Gamma[1:20])
  #print("------------------------------------")
  #print(Beta[Beta != 0])
  
  
  
	p <- dim(X)[2]
	n <- dim(X)[1]
	l <- length(tildelogY)
	updateStep <- 1
	tXB <-  crossprod(t(D),crossprod(t(X), Beta))
	lam <- ((n^gamma)*lambda*alpha)/(eta)

	inter.count <- 0
	for(lll in 1:max.iter){

		# ---------------------------------
		# Theta update
		# ---------------------------------
		tTheta <- Theta
		t0 <- tildelogY - tXB - Gamma/rho
		
		Theta <- ThetaUpdate(as.matrix(t0), tildedelta_nrho = tildedelta/((n^(2 - gamma))*rho))

		# -------------------------------------
		# Beta update
		# -------------------------------------
		W <- crossprod(X, crossprod(D, t0 - Theta))/eta + Beta
		
		for (j in 1:length(unique(groups))) {
			t0 <- pmax(abs(W[which(groups==j)]) - (lam*w[which(groups==j)])/rho, 0)*sign(W[which(groups==j)])
			t1 <- sqrt(sum(t0^2))
			
			if (t1 > 0) {
				Beta[which(groups==j)] <- max(1 - v[j]*lambda*(1-alpha)/(eta*rho*t1), 0)*t0
			}
		}
		#W <- Matrix(pmax(abs(W[,1]) - lam/rho, 0)*sign(W[,1]), sparse=TRUE)
		tXB <- tcrossprod(D, t(crossprod(t(X), Beta)))
		# ----------------------------------
		# Gamma
		# ----------------------------------
		Gamma <- Gamma + tau*rho*(Theta - tildelogY + tXB)

		if(lll%%floor(updateStep) == 0){

	        s <- rho*sqrt(sum(crossprod(X, SpMatMultDZ(D, as.matrix(Theta - tTheta)))^2))#/length(tildelogY)
			r <- sqrt(sum((Theta - tildelogY + tXB)^2))#/length(tildelogY)

			eprim <- sqrt(l)*tol.abs + tol.rel*max(c(sqrt(sum(tXB^2)), sqrt(sum(Theta^2)), euc.tildelogY))
			edual <- sqrt(p)*tol.abs + tol.rel*sqrt(sum(crossprod(X, crossprod(D, Gamma))^2))
		
			if(r/eprim > 10*s/edual){
				rho <- rho*2
			} 
			if(s/edual > 10*r/eprim){
				rho <- rho/2
			}
		
			#cat("r = ", r, ": s= ", s, "\n")
			if(lll > 10){
				if(r < eprim & s < edual){
					break
				}
			}
			updateStep <- (updateStep + 1)*1.1

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