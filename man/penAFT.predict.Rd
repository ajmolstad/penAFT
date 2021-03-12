\name{penAFT.predict}
\alias{penAFT.predict}
\title{Predict for new subjects using fitted model from \code{penAFT} of \code{penAFT.cv}}
\description{A function for prediction along the solution path of the regularized Gehan estimator.}
\usage{
penAFT.predict(fit, Xnew, lambda = NULL)
}
\arguments{
  \item{fit}{A fitted model from \code{penAFT} or \code{penAFT.cv}.}
  \item{Xnew}{A matrix of dimension \eqn{n_{\rm new} \times p}. Must be a matrix, even if \eqn{n_{\rm new}= 1}.}
  \item{lambda}{The value of \eqn{\lambda} used to estimate \eqn{\beta}. If \code{NULL} and \code{fit} was obtained using \code{nfolds} non-\code{NULL}, will use tuning parameters which minimized prediction error in cross-validation.  }
}
\value{
  \item{preds}{The matrix of predictors: rows correspond to rows of \code{Xnew}.  }
}



\details{
	See \url{ajmolstad.github.io/docs/penAFT_Example.html} for usage examples.
}

