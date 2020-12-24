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
#   standardize: should predictors be standardize for model fitting? 
# --------------------------------------------------------------------

penAFT.cv <- function(X, logY, delta, 
                   nlambda = 50, 
                   lambda.ratio.min = NULL, 
                   lambda = NULL, 
                   penalty = NULL,
                   alpha = 1, 
                   weight.set = NULL, 
                   groups = NULL, 
                   tol.abs = 1e-10, 
                   tol.rel = 5e-4, 
                   center = TRUE, 
                   standardize = FALSE,
                   nfolds = 5, 
                   fold.id = NULL,
                   admm.max.iter = 1e4, 
                   quiet = TRUE, 
                   cpp = TRUE) {
                     
  # ----------------------------------------------------------
  # Preliminary checks 
  # ----------------------------------------------------------
  p <- dim(X)[2]
  n <- dim(X)[1]
  gamma <- 0

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
  if (center & !standardize) {
    X.fit <- X - tcrossprod(rep(1, n), colMeans(X))
  } 
  if (standardize) {
    X.fit <- (X - tcrossprod(rep(1, n), colMeans(X)))/(tcrossprod(rep(1, n), apply(X, 1, sd)))
  }
  if (!center & !standardize) {
    X.fit <- X
  }

  if (is.null(penalty)) {
     penalty <- "EN"
  }

# -------------------------------------------------
# Elastic net fit
# ------------------------------------------------
  if (penalty == "EN") {

    if (is.null(weight.set)){
      w <- rep(1, p)
    } else {
      if (is.null(weight.set$w)) {
        stop("Need to specify weight.set as a list with element w to use weighted elastic net")
      } else {
        w <- weight.set$w
      }
    }
    
    # -----------------------------------------------------------
    # Get candidate tuning parameters
    # -----------------------------------------------------------
    if (is.null(lambda)) {
      if (alpha != 0) {
        EN_TPcalc <- function(X.fit, logY, delta) {
          n <- nrow(logY)
          grad <- rep(0, dim(X.fit)[2])
          grad2 <- rep(0, dim(X.fit)[2])
          for (i in which(delta==1)) {
            for (j in 1:n) {
              if(logY[i] != logY[j]){
                grad <- grad + delta[i]*(X.fit[i,] - X.fit[j,])*(logY[i] < logY[j])
              } else {
                grad2 <- grad2 + delta[i]*abs(X.fit[i,] - X.fit[j,])
              } 
            }
          }
          return(abs(grad)/n^2 + grad2/n^2)
        }
        gradG <- EN_TPcalc(X.fit, logY, delta)
        wTemp <- w
        if(any(wTemp==0)){
          wTemp[which(wTemp == 0)] <- Inf
        }
        lambda.max <- max(gradG/wTemp)/alpha + 1e-6
        if (is.null(lambda.ratio.min)) {
          lambda.ratio.min <- 0.1
        }
        lambda.min <- lambda.ratio.min*lambda.max
        lambda <- 10^seq(log10(lambda.max), log10(lambda.min), length=nlambda) 
      } else {
        lambda <- 10^seq(-4, 4, length=nlambda)
        warning("Setting alpha = 0 corresponds to a ridge-penalty: may need to check candidate tuning parameter values.")
      }
      
    } else {
      warning("It is recommended to let tuning parameters be chosen automatically: see documentation.")
    }
    
    # -----------------------------------------------
    # Fit the solution path for the elastic net
    # -----------------------------------------------
    getPath <- ADMM.ENpath(X.fit, logY, delta, admm.max.iter, lambda, alpha, w, tol.abs, tol.rel, gamma, quiet)

    # --------------------------------------------------------------------
    # Objective function to evaluate CV metrics
    # --------------------------------------------------------------------
    eval.obj.inner <- function(logY, XB, delta){
      out <- 0 
      n <- length(logY)
      E <- logY - XB
      for(i in which(delta==1)){
        out <- out + sum(pmax(E - E[i], 0))
      }
      return(out/n^2)
    }

    # --------------------------------------------------------------------
    # Perform cross-validation
    # --------------------------------------------------------------------
    fold1 <- sample(rep(1:nfolds, length=length(which(delta==1))))
    fold0 <- sample(rep(1:nfolds, length=length(which(delta==0))))

    cv.index1 <- split(which(delta==1), fold1)
    cv.index0 <- split(which(delta==0), fold0)
    cv.index <- list()
    for (ll in 1:nfolds) {
      cv.index[[ll]] <- c(cv.index1[[ll]],cv.index0[[ll]])
    }

    preds <- matrix(Inf, nrow = n, ncol = length(lambda))
    cv.err.obj <- matrix(Inf, nrow=nfolds, ncol=length(lambda))

    for (k in 1:nfolds) {
      # -----------------------------------------------------------
      # Get training predictors and responses
      # -----------------------------------------------------------
      ntrain <- dim(X[-cv.index[[k]],])[1]
      ntest <- length(cv.index[[k]])
      if (center & !standardize) {
        X.train.fit <- X[-cv.index[[k]], ] - tcrossprod(rep(1, ntrain), colMeans(X[-cv.index[[k]], ]))
        X.test <- X[cv.index[[k]], ] - tcrossprod(rep(1, ntest), colMeans(X[-cv.index[[k]], ]))
      } 
      if (standardize) {
        X.train.fit <- (X[-cv.index[[k]], ] - tcrossprod(rep(1, ntrain), colMeans(X[-cv.index[[k]], ])))/(tcrossprod(rep(1, ntrain), apply(X[-cv.index[[k]], ], 1, sd)))
        X.test <- (X[cv.index[[k]], ] - tcrossprod(rep(1, ntest), colMeans(X[-cv.index[[k]], ])))/(tcrossprod(rep(1, ntest), apply(X[-cv.index[[k]], ], 1, sd)))
      }
      if (!center & !standardize) {
        X.train.fit <- X[-cv.index[[k]], ]
        X.test <- X[cv.index[[k]], ]
      }

      logY.train <- logY[-cv.index[[k]]]
      logY.test <- logY[cv.index[[k]]]
      delta.train <- delta[-cv.index[[k]]]
      delta.test <- delta[cv.index[[k]]]

      cv.getPath <- ADMM.ENpath(X.train.fit, logY.train, delta.train, admm.max.iter, lambda, alpha, w, tol.abs, tol.rel, gamma, quiet)
      
      # ----------------------------------------------------------
      # cross-validated linear predictors 
      # ----------------------------------------------------------
      for (ll in 1:length(lambda)) {
        preds[cv.index[[k]],ll] <- X.test%*%cv.getPath$beta[,ll]
        cv.err.obj[k, ll] <- eval.obj.inner(logY.test, preds[cv.index[[k]],ll], delta.test)
      }
      cat("CV through: ", rep("###", k), rep("   ", nfolds - k), (k/nfolds)*100, "%", "\n")
    }
      
    cv.err.linPred <- rep(Inf, length(lambda))
    for (ll in 1:length(lambda)) {
       cv.err.linPred[ll] <- eval.obj.inner(logY, preds[,ll], delta)
    }


  } else {

# ------------------------------------------------------
#
# ------------------------------------------------------

      if (penalty == "SG") {
        if (is.null(groups)) {
          stop("To use group-lasso penalty, must specify \"groups\"!")
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
      if (length(unique(logY)) == length(logY)) {

        getGrad <- function(X.fit, logY, delta) {
          n <- nrow(logY)
          grad <- rep(0, dim(X.fit)[2])
          for (i in which(delta==1)) {
            for (j in 1:n) {
                grad <- grad + delta[i]*(X.fit[i,] - X.fit[j,])*(logY[i] < logY[j])
              } 
          }
          return(grad/n^2)
        }

        grad <- getGrad(X.fit, logY, delta)

        lam.check <- 10^seq(-4, 4, length=500)
        check.array <- matrix(0, nrow=length(lam.check), ncol=G)
        for(j in 1:length(lam.check)){
          for(k in 1:G){
            t0 <- pmax(abs(grad[which(groups==k)]) - alpha*lam.check[j]*w[which(groups==k)], 0)*sign(grad[which(groups==k)])
            check.array[j,k] <- sqrt(sum(t0^2)) < v[k]*(1-alpha)*lam.check[j] 
          }
        }

      lambda.max <- lam.check[min(which(rowSums(check.array) == G))] + 1e-6
      if (is.null(lambda.ratio.min)) {
        lambda.ratio.min <- 0.1
      }
      lambda.min <- lambda.ratio.min*lambda.max
      lambda <- 10^seq(log10(lambda.max), log10(lambda.min), length=nlambda)

      } else {
        stop("Ties present in the input logY!")
      }
    } else {
      warning("It is recommended to let tuning parameters be chosen automatically: see documentation.")
    }


    # ------------------------------------------------------
    # Fit the solution path for the sparse group lasso
    # ------------------------------------------------------
    getPath <- ADMM.SGpath(X.fit, logY, delta, admm.max.iter, lambda, alpha, w, v, groups, tol.abs, tol.rel, gamma, quiet)

    # --------------------------------------------------------------------
    # Objective function to evaluate CV metrics
    # --------------------------------------------------------------------
    eval.obj.inner <- function(logY, XB, delta){
      out <- 0 
      n <- length(logY)
      E <- logY - XB
      for(i in which(delta==1)){
        out <- out + sum(pmax(E - E[i], 0))
      }
      return(out/n^2)
    }

    # --------------------------------------------------------------------
    # Perform cross-validation
    # --------------------------------------------------------------------
    fold1 <- sample(rep(1:nfolds, length=length(which(delta==1))))
    fold0 <- sample(rep(1:nfolds, length=length(which(delta==0))))

    cv.index1 <- split(which(delta==1), fold1)
    cv.index0 <- split(which(delta==0), fold0)
    cv.index <- list()
    for (ll in 1:nfolds) {
      cv.index[[ll]] <- c(cv.index1[[ll]],cv.index0[[ll]])
    }

    preds <- matrix(Inf, nrow = n, ncol = length(lambda))
    cv.err.obj <- matrix(Inf, nrow=nfolds, ncol=length(lambda))

    for (k in 1:nfolds) {
      # -----------------------------------------------------------
      # Get training predictors and responses
      # -----------------------------------------------------------
      ntrain <- dim(X[-cv.index[[k]],])[1]
      ntest <- length(cv.index[[k]])
      if (center & !standardize) {
        X.train.fit <- X[-cv.index[[k]], ] - tcrossprod(rep(1, ntrain), colMeans(X[-cv.index[[k]], ]))
        X.test <- X[cv.index[[k]], ] - tcrossprod(rep(1, ntest), colMeans(X[-cv.index[[k]], ]))
      } 
      if (standardize) {
        X.train.fit <- (X[-cv.index[[k]], ] - tcrossprod(rep(1, ntrain), colMeans(X[-cv.index[[k]], ])))/(tcrossprod(rep(1, ntrain), apply(X[-cv.index[[k]], ], 1, sd)))
        X.test <- (X[cv.index[[k]], ] - tcrossprod(rep(1, ntest), colMeans(X[-cv.index[[k]], ])))/(tcrossprod(rep(1, ntest), apply(X[-cv.index[[k]], ], 1, sd)))
      }
      if (!center & !standardize) {
        X.train.fit <- X[-cv.index[[k]], ]
        X.test <- X[cv.index[[k]], ]
      }

      logY.train <- logY[-cv.index[[k]]]
      logY.test <- logY[cv.index[[k]]]
      delta.train <- delta[-cv.index[[k]]]
      delta.test <- delta[cv.index[[k]]]

      cv.getPath <- ADMM.ENpath(X.train.fit, logY.train, delta.train, admm.max.iter, lambda, alpha, w, v, groups, tol.abs, tol.rel, gamma, quiet)
      
      # ----------------------------------------------------------
      # cross-validated linear predictors 
      # ----------------------------------------------------------
      for (ll in 1:length(lambda)) {
        preds[cv.index[[k]],ll] <- X.test%*%cv.getPath$beta[,ll]
        cv.err.obj[k, ll] <- eval.obj.inner(logY.test, preds[cv.index[[k]],ll], delta.test)
      }
      cat("CV through: ", rep("###", k), rep("   ", nfolds - k), (k/nfolds)*100, "%", "\n")
    }
      
    cv.err.linPred <- rep(Inf, length(lambda))
    for (ll in 1:length(lambda)) {
       cv.err.linPred[ll] <- eval.obj.inner(logY, preds[,ll], delta)
    }




    }
  } 

 
  return(list("full.fit" = getPath, 
    "cv.err.linPred" = cv.err.linPred,
    "cv.err.obj" = cv.err.obj))
  
}