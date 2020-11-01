ADMM.SGrun <- function(tildelogY, X, D, tildedelta, rho, eta, tau, lambda, alpha, w, v, border.indexes, Gamma, Beta, Theta, max.iter, tol.abs, tol.rel, gamma, euc.tildelogY, G) {
  p <- dim(X)[2]
  n <- dim(X)[1]
  l <- length(tildelogY)
  
  
  Xbeta <- crossprod(t(X), Beta)
  tXB <-  crossprod(t(D),Xbeta)
  tXB <- Matrix(tXB, sparse = FALSE)
  #D <- Matrix(D, sparse=TRUE)
  
  
  #print(eta)
  #print("-----------------------------")
  #print(rho)
  #print("------------------------------------")
  #print(Theta[1:20])
  #print("------------------------------------")
  #print(Gamma[1:20])
  #print("------------------------------------")
  #print(Beta[Beta != 0])
  
  
  out <- ADMM_SGrun(tildelogY, X, D, tildedelta, rho, eta, tau, lambda, alpha, w, v, borderIndexes = border.indexes, Gamma, Beta, Theta, max_iter = max.iter, tol_abs = tol.abs, tol_rel = tol.rel, gamma, euc_tildelogY = euc.tildelogY, n = n, l = l, p = p, G = G)
  out
}
