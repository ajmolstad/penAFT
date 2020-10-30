ADMM.SGrun <- function(tildelogY, X, D, tildedelta, rho, eta, tau, lambda, alpha, w, v, border.indexes, Gamma, Beta, Theta, max.iter, tol.abs, tol.rel, gamma, euc.tildelogY, max.iter.update, G) {
  p <- dim(X)[2]
  n <- dim(X)[1]
  l <- length(tildelogY)


  Xbeta <- crossprod(t(X), Beta)
  tXB <-  crossprod(t(D),Xbeta)
  tXB <- Matrix(tXB, sparse = FALSE)
  #D <- Matrix(D, sparse=TRUE)

  out <- ADMM_SGrun(tildelogY, X, D, tildedelta, rho, eta, tau, lambda, alpha, w, v, borderIndexes = border.indexes, Gamma, Beta, Theta, max_iter = max.iter, tol_abs = tol.abs, tol_rel = tol.rel, gamma, euc_tildelogY = euc.tildelogY, n, l, p, max_iter_update = max.iter.update, G)
  out
}
