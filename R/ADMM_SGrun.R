ADMM.SGrun <- function(tildelogY, X, D, tildedelta, rho = rho, eta = eta, tau = 1.5, 
           lambda, alpha, w, v, border.indexes, Gamma, Beta, Theta, max.iter = 5000, 
           tol.abs, tol.rel, gamma, euc.tildelogY) {
  
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
    A <- (1/(eta))*crossprod(X, crossprod(D, tildelogY - Theta - Gamma/rho - tXB)) + Beta
    #Beta <- Matrix(pmax(abs(W[,1]) - ((n^gamma)*lambda*alpha*w)/(eta*rho), 0)*sign(W[,1]), sparse=TRUE)/(1 + ((n^gamma)*lambda*(1-alpha)*w)/(eta*rho))
    #tXB <- crossprod(t(D), crossprod(t(X), Beta))
    
    for (g in 1:(length(border.indexes) - 1)) {
      i <- border.indexes[g]
      j <- border.indexes[g+1] - 1
      
      a.vec <- A[i:j]
      tau.vec <- w[i:j] * alpha * (n^gamma) * (lambda / (eta * rho))
      
      soft.vec <- pmax(abs(a.vec) - tau.vec, 0) * sign(a.vec)
      soft.denom <- sqrt(sum(soft.vec^2))
      r0 <- (v[g] * (n^gamma) * (lambda / (eta * rho)) * (1 - alpha)) / soft.denom
      r1 <- max(1 - r0, 0)
      
      #print(max(soft.vec))
      Beta[i:j] <- r1 * soft.vec
    }
    
    tXB <- crossprod(t(D), crossprod(t(X), Beta))
    
    # ----------------------------------
    # Gamma
    # ----------------------------------
    Gamma <- Gamma + tau*rho*(Theta - tildelogY + tXB)
    
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