\name{penAFT.plot}
\alias{penAFT.plot}
\title{Plot cross-validation curves}
\description{A function for plotting the cross-validation curves for the regularized Gehan estimator.}
\usage{
penAFT.plot(fit)
}

\arguments{
\item{fit}{A fitted model from \code{penAFT.cv}.}
}


\details{
	This function returns a plot with the cross-validation curves for the regularized Gehan estimator.  The vertical blue line indicates the tuning parameter which minimized cross-validated linear predictor error and the vertical black line indicates the tuning parameter minimized the cross-validated deviance. The vertical axis (and blue line) denotes the cross-validated linear predictor error whereas the right vertical axis (and black line) denotes cross-validated deviance. To make matters simple, we do not allow for customization of the plot: please refer to the source code if extensive customization is desired. 
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


# ------------------------------------------
# Fit elastic net penalized estimator with CV
# -------------------------------------------
fit.cv <- penAFT.cv(X = X, logY = logY, delta = delta,
                   nlambda = 50,
                   penalty = "EN", tol.rel= 1e-5,
                   alpha = 1, nfolds = 5)

# penAFT.plot(fit.cv)

}