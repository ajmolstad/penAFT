\name{penAFT}
\alias{penAFT}
\title{Fit the solution path for the regularized Gehan estimator.}
\description{A function to fit the solution path for the regularized Gehan estimator.}
\usage{

penAFT(X, logY, delta,
                   nlambda = 50,
                   lambda.ratio.min = NULL, lambda = NULL,
                   penalty = NULL,
                   alpha = 1, weight.set = NULL,
                   groups = NULL, tol.abs = 1e-10,
                   tol.rel = 5e-4,
                   gamma = 0, center = TRUE,
                   standardize = FALSE,
                   admm.max.iter = 1e4, quiet=TRUE)
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
  \item{gamma}{A non-negative optimization parameter which can improve convergence speed in certain settings. It is recommended to set equal to zero.}
  \item{center}{Should predictors be centered for model fitting?}
  \item{standardize}{Should predictors be standardized (i.e., scaled to have unit variance) for model fitting?}
  \item{admm.max.iter}{Maximum number of ADMM iterations.}
  \item{quiet}{\code{TRUE} or \code{FALSE} variable indicating whether progress should be printed.}
}

\value{
  \item{beta}{A \eqn{p \times} \code{nlambda} sparse matrix consisting of the estimates of \eqn{\beta} for the candidate values of \eqn{\lambda}. It is recommended to use \code{penAFT.coef} to extract coefficients.}
  \item{lambda}{The candidate tuning parameter values. }
  \item{center}{Were predictors centered for model fitting? }
  \item{standardize}{Were predictors standardized to have unit variance for model fitting? }
  \item{X.mean}{The mean of the predictors.}
  \item{X.sd}{The standard deviation of the predictors.}
  \item{alpha}{The tuning parameter \eqn{\alpha}. See documentation.}
}


\details{
  See \url{ajmolstad.github.io/docs/penAFT_Example.html} for usage examples.
}
