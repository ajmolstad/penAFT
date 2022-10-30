\name{penAFT.predict}
\alias{penAFT.predict}
\title{Obtain linear predictor for new subjects using fitted model from \code{penAFT} of \code{penAFT.cv}}
\description{A function for prediction along the solution path of the regularized semiparametric accelerated failure time model estimator.}
\usage{
penAFT.predict(fit, Xnew, lambda = NULL)
}
\arguments{
  \item{fit}{A fitted model from \code{penAFT} or \code{penAFT.cv}.}
  \item{Xnew}{A matrix of dimension \eqn{n_{\rm new} \times p}. Must be a matrix, even if \eqn{n_{\rm new}= 1}.}
  \item{lambda}{The value of \eqn{\lambda} used to estimate \eqn{\beta}. If \code{NULL} and \code{fit} was obtained using \code{nfolds} non-\code{NULL}, the function will use the tuning parameter which minimized cross-validation linear predictor scores.  }
}
\value{
  \item{preds}{The matrix of linear predictors: rows correspond to rows of \code{Xnew}.  }
}


\details{
	It is important to note that the output of this function should not be treated as an estimate of the log-survival time. Because the Gehan loss function is location invariant, the intercept is absored into the error. If predictors were standardized for model fitting, this function returns \eqn{\tilde{X}_{\rm new}\hat{\beta}} where \eqn{\tilde{X}_{\rm new}} is the version of input \code{Xnew} which has been centered and standardized according to the design matrix used to fit the \code{penAFT} or \code{penAFT.cv} object. If predictors were not standardized, this function returns \eqn{X_{\rm new}\hat{\beta}}. 
	
 We recommend input \code{Xnew} as a matrix, although if a \eqn{p}-dimensional vector is input, the function will detect this. 
}

\examples{
# --------------------------------------
# Generate data  
# --------------------------------------
set.seed(1)
genData <- genSurvData(n = 50, p = 50, s = 10, mag = 2, cens.quant = 0.6)
X <- genData$X
logY <- genData$logY
delta <- genData$status

# --- generate data for two new subjects
p <- dim(X)[2]
Xnew <- rbind(rnorm(p), rnorm(p))

# -----------------------------------------------
# Fit elastic net penalized estimator without CV
# -----------------------------------------------
fit <- penAFT(X = X, logY = logY, delta = delta,
                   nlambda = 10, lambda.ratio.min = 0.1,
                   penalty = "EN",
                   alpha = 1)

# predict at 10th candidate tuning parameter
linPred.10 <- penAFT.predict(fit, Xnew = Xnew, lambda = fit$lambda[10])

\donttest{
  # ------------------------------------------
  # Fit elastic net penalized estimator with CV
  # -------------------------------------------
  fit.cv <- penAFT.cv(X = X, logY = logY, delta = delta,
                    nlambda = 50,
                    penalty = "EN",
                     alpha = 1, nfolds = 5)
  
  # --- return linear predictor at lambda minimizing cross-validation error 
  linPred.cv <- penAFT.predict(fit.cv, Xnew = Xnew) 
  
  # --- predict at 10th candidate tuning parameter
  linPred.cv10 <- penAFT.predict(fit.cv, Xnew = Xnew, lambda = fit.cv$full.fit$lambda[10]) 
  
  
  # ------------------------------------------
  # Fit penAFT with cross-validation
  # -------------------------------------------
  groups <- rep(1:5, each = 10)
  fit.sg.cv <- penAFT.cv(X = X, logY = logY, delta = delta,
                    nlambda = 50, groups = groups,
                    penalty = "SG",
                    alpha = 0.5, nfolds = 5)
  
  # ---- return linear predictor at lambda minimizing cross-validation error 
  linPred.sg.cv <- penAFT.predict(fit.sg.cv, Xnew = Xnew) 
  
  # --- predict at 10th candidate tuning parameter
  linPred.sg.cv10 <- penAFT.predict(fit.sg.cv, Xnew = Xnew, lambda = fit.sg.cv$full.fit$lambda[10]) 

}

}