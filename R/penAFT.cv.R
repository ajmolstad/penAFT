# ----------------------------------------------------------------------
# penAFT function for fitting solution path
# ----------------------------------------------------------------------
# Args:
#   X: n times p design matrix 
#   logY: n dimensional vector of log survival times or censoring times 
#   delta: n dimensional vector censoring indicator (1 = uncensored, 0 = censored)
#   nlambda: number of candidate tuning parameters to consider
#   lambda.ratio.min: ratio of largest to small cnadidate tuning parmaeter (e.g., 0.1)
#   lambda: a vector of candidat tuning parameters that can be used to override internal choice
#   penalty: which penalty should be used? EN is elastic net, SG is sparse group lasso
#   alpha: balance parameter between 0 and 1 -- has a different meaning depending on the penalty
#   weights: a list containing w and v: for pen = EN, only w is used, for pen = SG, both w and v are used. If either is not input, all weights are set to 1
#   groups: if pen = SG, a p-dimensional vector of integers indicating group membership
#   tol.abs: ADMM absolute convergence tolerance
#   tol.rel: ADMM relative convergence 
#   gamma: ADMM balance parameter
#   centered: should predictors be centered for model fitting?
#   standardized: should predictors be standardized for model fitting? 
# --------------------------------------------------------------------

penAFT.cv <- function(X, logY, delta, 
                   nlambda = 50, 
                   lambda.ratio.min = NULL, lambda = NULL, 
                   penalty = NULL,
                   alpha = 1, weights = NULL, 
                   groups = NULL, tol.abs = 1e-10, 
                   tol.rel = 5e-4, 
                   gamma = 0, centered = TRUE, 
                   standardized = FALSE,
                   nfolds = 5, fold.id = NULL) {
                     
  # ----------------------------------------------------------
  # Preliminary checks
  # ----------------------------------------------------------
  p <- dim(X)[2]
  n <- dim(X)[1]
  if (length(logY) != n) {
    stop("Dimensions of X and logY do not match: see documentation")
  }
  if (length(delta) != n) {
    stop("Dimension of X and delta do not match: see documentation")
  }
  if (!any(delta==0 | delta==1)) {
    stop("delta  must be a vector with elements in {0,1}: see documentation")
  }
  if (alpha > 1 | alpha < 0) {
    stop("alpha must be in [0,1]: see documentation.")
  }
  
  if(length(unique(logY))!=length(logY)){
    stop("logY contains duplicate survival or censoring times (i.e., ties).")
  }
  
  # -----------------------------------------------------------
  # Center and standardize
  # -----------------------------------------------------------
  if(center & !standardize){
    X.fit <- X - tcrossprod(rep(1, n), colMeans(X))
  } 
  if(standardize){
    X.fit <- (X - tcrossprod(rep(1, n), colMeans(X)))/(tcrossprod(rep(1, n), apply(X, 1, sd)))
  }
  if(!center& !standardize){
    X.fit <- X
  }
  
  # -----------------------------------------------------------
  # Get candidate tuning parameters
  # -----------------------------------------------------------
  gradient_g <- function(X.fit, logY, beta, delta) {
    n <- nrow(X.fit)
    p <- ncol(X.fit)
    grad <- rep(0, p)
    Xbeta <- crossprod(t(X.fit), beta)
    E <- logY - Xbeta
    for (i in 1:n) {
      for (j in 1:n) {
        term <- delta[i]*(X.fit[i,] - X.fit[j,])*(E[i] <= E[j])
        grad <- grad + term
      }
    }
    return(grad/n^2)
  }
  gradG <- gradient_g(X.fit, logY, rep(0, p), delta)
  
  if (penalty == "SG"){
    if (is.null(groups)) {
      stop("To use group-lasso penalty, must specify 'groups'!")
    } 
    G <- length(unique(groups))
    if(is.null(weights)){
      w <- rep(1, p)
      v <- rep(1, G)
    } else {
      if (is.null(weights$w)) {
        stop("Need to specify both w and v using weighted group lasso")
      } else {
        w <- weights$w
      }
      if (is.null(weights$v)) {
        stop("Need to specify both w and v using weighted group lasso")
      } else {
        v <- weights$v
      }
    }
    
    if (is.null(lambda)) {
     lambda.max <- max(abs(gradG/w))/alpha
     if (is.null(lambda.ratio.min)) {
       lambda.ratio.min <- 0.01
     }
     lambda.min <- lambda.ratio.min*lambda.max
     lambda <- 10^seq(log10(lambda.max), log10(lambda.min), length=nlambda)
      
    } else {
      warning("It is recommended to let tuning parameters be chosen automatically: see documentation.")
    }
    
    # ------------------------------------------------------
    # Fit the solution path for the sparse group lasso
    # ------------------------------------------------------
    getPath <- ADMM.SGpath(X.fit, logY, delta, lambda, w, v, groups, tol.abs, tol.rel, gamma)
    
  } 
  
  
  if (pen == "EN" | is.null(penalty)){
    if (is.null(penalty)) {
      warning("No penalty specified: using elastic net!")
    }
    if(is.null(weights)){
      w <- rep(1, p)
    } else {
      if (is.null(weights$w)) {
        stop("Need to specify both w to use weighted elastic net")
      } else {
        w <- weights$w
      }
    }
    
    if (is.null(lambda)) {
      lambda.max <- max(abs(gradG/w))/alpha
      if (is.null(lambda.ratio.min)) {
        lambda.ratio.min <- 0.1
      }
      lambda.min <- lambda.ratio.min*lambda.max
      lambda <- 10^seq(log10(lambda.max), log10(lambda.min), length=nlambda)
      
    } else {
      warning("It is recommended to let tuning parameters be chosen automatically: see documentation.")
    }
    
    # -----------------------------------------------
    # Fit the solution path for the elastic net
    # -----------------------------------------------
    getPath <- ADMM.ENpath(X.fit, logY, delta, lambda, w, v, groups, tol.abs, tol.rel, gamma)
    
  }     
  
  
  # --------------------------------------------
  # Append useful info to getpath
  # --------------------------------------------
  getPath$center
  getPath$standardize
  getPath$X.mean
  getPath$X.sd 
  
  return(getPath)
}