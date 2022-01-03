\name{penAFT.trace}
\alias{penAFT.trace}
\title{Print trace plot for the regularized Gehan estimator fit using \code{penAFT} or \code{penAFT.cv}}
\description{Print the trace plot for the regularized Gehan estimator.}
\usage{
penAFT.trace(fit, groupNames=NULL)
}

\arguments{
\item{fit}{A fitted model from \code{penAFT} or \code{penAFT.cv}.}
\item{groupNames}{A list of groupnames to be printed when \code{fit} used penalty "SG".}
}


\details{
	The function \code{penAFT.trace} returns a trace plot for a fitted model obtained from either \code{penAFT} or \code{penAFT.cv}. If the model is fit using the sparse group-lasso penalty, you may provide names for the groups (in order of the integer values specifying the groups). This feature may not be desired if there are a large number of groups, however. The vertical blue line indicates the tuning parameter which minimized cross-validated linear predictor scores and the vertical black line indicates the tuning parameter minimizing the cross-validated Gehan loss according to the one-standard error rule. 
}

\examples{
# --------------------------------------
# Generate data  
# --------------------------------------
set.seed(1)
genData <- genSurvData(n = 50, p = 50, s = 25, mag = 0.5, cens.quant = 0.6)
X <- genData$X
logY <- genData$logY
delta <- genData$status


# ------------------------------------------
# Fit elastic net penalized estimator with CV
# -------------------------------------------
fit.cv <- penAFT.cv(X = X, logY = logY, delta = delta,
                   nlambda = 10, tol.rel= 1e-5,
                   penalty = "EN",
                   alpha = 1, nfolds = 5)

# penAFT.trace(fit.cv)


# ------------------------------------------
# Fit sparse group-lasso estimator with CV
# -------------------------------------------
# groups <- rep(1:5, length=100)
# fit.sg.cv <- penAFT.cv(X = X, logY = logY, delta = delta,
#                   nlambda = 50, groups = groups,
#                   penalty = "SG", tol.rel= 1e-5,
#                   alpha = 0, nfolds = 5)

# penAFT.trace(fit.sg.cv, groupNames = paste("Group", 1:5, sep="-"))

}
