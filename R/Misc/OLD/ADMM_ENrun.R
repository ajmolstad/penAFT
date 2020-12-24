ADMM.ENrun <- function(tildelogY, X, D, D.pos, tildedelta, rho, eta, tau, 
                       lambda, alpha, w, Gamma, Beta, 
                       Theta,  
                       max.iter, tol.abs, tol.rel, gamma, euc.tildelogY){


  out <- ADMM_ENrun(tildelogY, X, D, D.pos, tildedelta, rho, eta, tau, 
             lambda, alpha, w, Gamma, Beta, 
             Theta,  
             max.iter, tol.abs, tol.rel, gamma, euc.tildelogY)
  
  if(out$iter.counter == max.iter){
    warning("ADMM did not converge in max.iter iterations", "\n")
  }
  
  
  out
}