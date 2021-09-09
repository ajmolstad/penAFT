\name{genSurvData}
\alias{genSurvData}
\title{Generate a survival dataset from the log-logistic accelerated failure time model.}
\description{This is a function for generating synthetic datasets from the log-logistic accelerated failure time model. The purpose of this function is to provide structured data for the examples of the other functions' usage. }
\usage{
genSurvData(n, p, s, cens.quant = 0.6)
}

\arguments{
\item{n}{The numer of subjects to be included in the dataset.}
\item{p}{Dimension of the predictor. Note that the function computes the square-root of a \eqn{p \times p} covariance matrix, so setting \code{p} large may be time-consuming.}
\item{s}{The number of coefficients of \eqn{\beta} to be set equal to one. All others will be zero.}
\item{cens.quant}{The quantile of true survival times used to set the mean of the exponential distribution from which censoring times are drawn. Default is 0.6. }
}


\details{
 This function generates predictors to follow a \eqn{p}-dimensional multivariate normal distribution whose covariance has an AR(1) structure with lag 0.7. Then, log survival times are generated as
 \deqn{\log(T) = X \beta + \epsilon}
 where \eqn{\epsilon} has independent components drawn from a logistic distribution with location parmeter zero and scale parameter two. Then censoring times are drawn from an exponential distribution with mean equal to the quantile \code{cens.quant} of \eqn{T}.
 
 For example usage, see the example usage of \code{penAFT}. 
}
