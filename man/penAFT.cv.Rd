\name{penAFT.cv}
\alias{penAFT.cv}
\title{Cross-validation function for fitting regularized Gehan estimator}
\description{A function to perform cross-validation and compute the solution path for the regularized Gehan estimator.}
\usage{
penAFT.cv(X, logY, delta, 
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
                   cv.index = NULL,
                   admm.max.iter = 1e4, 
                   quiet = TRUE)
}

\arguments{
  \item{X}{An \eqn{n \times p} matrix of predictors. Observations should be organized by row.}
  \item{logY}{An \eqn{n}-dimensional vector of log-survival or log-censoring times.}
  \item{delta}{An \eqn{n}-dimensional binary vector indicating whether the \eqn{j}th component of \code{logY} is an observed log-survival time (\eqn{\delta_j = 1}) or a log-censoring time (\eqn{\delta_j = 0}) for \eqn{j=1, \dots, n}.}
  \item{nlambda}{The number of candidate tuning parameters to consider.}
  \item{lambda.ratio.min}{The ratio of maximum to minimum candidate tuning parameter value. As a default, we suggest 0.1, but standard model selection procedures should be applied to select \eqn{\lambda}.}
  \item{lambda}{An optional (not recommended) prespecified vector of candidate tuning parameters. Should be in descending order. }
  \item{penalty}{Either "EN" or "SG" for elastic net or sparse group lasso penalties.}
  \item{alpha}{The tuning parameter \eqn{\alpha}. See documentation.}
  \item{weight.set}{A list of weights. For both penalties, \eqn{w} is an \eqn{n}-dimensional vector of nonnegative weights. For "SG" penalty, can also include \eqn{v} -- a non-negative vector the length of the number of groups. See documentation for usage example.}
  \item{groups}{When using penalty "SG", a \eqn{p}-dimensional vector of integers corresponding the to group assignment of each predictor (i.e., column of \code{X}).}
  \item{tol.abs}{Absolute convergence tolerance.}
  \item{tol.rel}{Relative convergence tolerance.}
  \item{center}{Should predictors be centered for model fitting?}
  \item{standardize}{Should predictors be standardized (i.e., scaled to have unit variance) for model fitting?}
  \item{nfolds}{The number of folds to be used for cross-validation. Default is five.}
  \item{cv.index}{A list of length \code{nfolds} of indices to be used for cross-validation. This is to be used if trying to perform cross-validation for both \eqn{\alpha} and \eqn{\lambda}.}
  \item{admm.max.iter}{Maximum number of ADMM iterations.}
  \item{quiet}{\code{TRUE} or \code{FALSE} variable indicating whether progress should be printed.}
}

\value{
\item{full.fit}{A model fit with the same output as a model fit using \code{penAFT}. See documentation for \code{penAFT} for more.}
\item{cv.err.linPred}{A \code{nlambda}-dimensional vector of linear-predictor-based cross-valdiation errors.}
\item{cv.err.obj}{A \code{nfolds }\eqn{\times}\code{ nlambda} matrix of within-fold cross-valdiation errors.}
\item{cv.index}{A list of length \code{nfolds}. Each element contains the indices for subjects belonging to that particular fold.}
}

 
\details{
  See \url{ajmolstad.github.io/docs/penAFT_Example.html} for usage examples.
}
  