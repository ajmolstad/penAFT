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
 \item{mean.adjustment}{The mean of predictors used for adjustment in the case that \code{center=TRUE}.}

}



\details{
	See \url{ajmolstad.github.io/docs/penAFT_Example.html} for usage examples.
}
