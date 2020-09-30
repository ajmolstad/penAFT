ADMM.ENrun <- function(tildelogY, X, D, tildedelta, rho, eta, tau, lambda, alpha, w, Gamma, Beta, Theta, max.iter, tol.abs, tol.rel, gamma, euc.tildelogY) {
  p <- dim(X)[2]
  n <- dim(X)[1]
  l <- length(tildelogY)
  Xbeta <- crossprod(t(X), Beta)
  tXB <-  crossprod(t(D),Xbeta)
  #D <- Matrix(D, sparse=TRUE)

  out <- ADMM_ENrun(tildelogY, X, D, tildedelta, rho, eta, tau, lambda, alpha, w, Gamma, Beta, Theta, max_iter = max.iter, tol_abs = tol.abs, tol_rel = tol.rel, gamma, euc_tildelogY = euc.tildelogY, Xbeta, tXB, n, l, p)
  out
}
