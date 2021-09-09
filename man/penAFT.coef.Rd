\name{penAFT.coef}
\alias{penAFT.coef}
\title{Extract regression coefficients from fitted model object}
\description{A function to extract coefficients along the solution path for the regularized Gehan estimator.}
\usage{
penAFT.coef(fit, lambda = NULL)
}

\arguments{
\item{fit}{A fitted model from \code{penAFT} or \code{penAFT.cv}.}
\item{lambda}{The tuning parameter value at which to extract coefficients. If \code{NULL} and \code{fit} is a \code{penAFT.cv} object, will use the tuning parameter value with minimum cross-validation error.}
}

\value{
 \item{beta}{The coefficient estimates}
}



\details{
  The regression coefficients stored in the fitted model objects coming from \code{penAFT} or \code{penAFT.cv} will (i) be on the scale of standardized predictors if standardization was used (which is the default) and (ii) are stored as a specific sparse matrix so that coefficient extraction is cumbersome. This function returns the regression coefficient estimates on the original scale of the predictors for a particular tuning parmaeter value. It is important to note that this method does not return an estimate of the intercept: the intercept is absored into the error term as the Gehan loss function is invariant to location change. 
}

\examples{
# --------------------------------------
# Generate data  
# --------------------------------------
set.seed(1)
genData <- genSurvData(n = 50, p = 100, s = 25, cens.quant = 0.6)
X <- genData$X
logY <- genData$logY
delta <- genData$status


# --------------------------------------
# Fit elastic net penalized estimator without CV
# --------------------------------------
fit <- penAFT(X = X, logY = logY, delta = delta,
                   nlambda = 50,
                   penalty = "EN",
                   alpha = 1)

coef.10 <- penAFT.coef(fit, lambda = fit$lambda[10])
coef.20 <- penAFT.coef(fit, lambda = fit$lambda[20])

# Cannot obtain fit at lambda not in fit$lambda
# coef.error <- penAFT.coef(fit, lambda = 10) # throws error


# ------------------------------------------
# Fit elastic net penalized estimator with CV
# -------------------------------------------
fit.cv <- penAFT.cv(X = X, logY = logY, delta = delta,
                   nlambda = 50,
                   penalty = "EN",
                   alpha = 1, nfolds = 5)

# coefficients at lambda minimizing cross-validation error
coef.cv <- penAFT.coef(fit.cv) 

# coefficients at 10th considered lambda 
coef.cv10 <- penAFT.coef(fit.cv, lambda = fit.cv$full.fit$lambda[10]) 

# -------------------------------------------
# Repeat with sparse group lasso without CV
# -------------------------------------------
groups <- rep(1:10, each = 10)
fit.sg <- penAFT(X = X, logY = logY, delta = delta,
                   nlambda = 50, groups = groups,
                   penalty = "SG",
                   alpha = 0.5)

coef.sg.10 <- penAFT.coef(fit.sg, lambda = fit.sg$lambda[10])
coef.sg.20 <- penAFT.coef(fit.sg, lambda = fit.sg$lambda[20])


# -------------------------------------------
# Finally, fit sparse group lasso with CV
# -------------------------------------------
groups <- rep(1:10, each = 10)
fit.sg.cv <- penAFT.cv(X = X, logY = logY, delta = delta,
                   nlambda = 50, groups = groups,
                   penalty = "SG",
                   alpha = 0.5, nfolds = 5)

coef.sg.cv <- penAFT.coef(fit.sg.cv)
coef.sg.cv10 <- penAFT.coef(fit.sg.cv, lambda = fit.sg$full.fit$lambda[20])

}