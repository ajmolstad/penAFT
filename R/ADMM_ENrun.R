ADMM.ENrun <- function(tildelogY, X, D, tildedelta, rho, eta, tau, 
                       lambda, alpha, w, Gamma, Beta, 
                       Theta,  
                       max.iter, tol.abs, tol.rel, gamma, euc.tildelogY){
  max.iter <- 1
  out <- ADMM_ENrun(tildelogY, X, D, tildedelta, rho, eta, tau, 
             lambda, alpha, w, Gamma, Beta, 
             Theta,  
             max.iter, tol.abs, tol.rel, gamma, euc.tildelogY)
  
  if(out$iter.counter == max.iter){
    warning("ADMM did not converge in max.iter iterations", "\n")
  }
  
  print(out$Theta[out$Theta != 0][1:30])
  
  out
}