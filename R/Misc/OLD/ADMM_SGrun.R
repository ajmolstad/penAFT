ADMM.SGrun2 <- function(tildelogY, X, D.pos, D.vert.1, D.vert.neg1, tildedelta, rho, eta, tau, lambda, alpha, w, v, border.indexes, Gamma, Beta, Theta, max.iter, tol.abs, tol.rel, gamma, euc.tildelogY, G) {
  p <- dim(X)[2]
  n <- dim(X)[1]
  l <- length(tildelogY)
  out <- ADMM_SGrun2(tildelogY, X, D.pos, D.vert.1, D.vert.neg1, tildedelta, rho, eta, tau, lambda, alpha, w, v, borderIndexes = border.indexes, Gamma, Beta, Theta, max_iter = max.iter, tol_abs = tol.abs, tol_rel = tol.rel, gamma, euc_tildelogY = euc.tildelogY, n = n, l = l, p = p, G = G)
  out
}
