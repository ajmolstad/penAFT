\name{penAFT}
\alias{penAFT}
\title{Fit the solution path for the regularized Gehan estimator.}
\description{A function to fit the solution path for the regularized Gehan estimator.}
\usage{
penAFT(X, logY, delta, nlambda = 50, 
  lambda.ratio.min = 0.1, lambda = NULL, 
  penalty = NULL, alpha = 1, weight.set = NULL, 
  groups = NULL, tol.abs = 1e-8, tol.rel = 2.5e-4, 
  gamma = 0,  standardize = TRUE, 
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
  \item{gamma}{A non-negative optimization parameter which can improve convergence speed in certain settings. It is highly recommended to set equal to zero.}
  \item{standardize}{Should predictors be standardized (i.e., column-wise average zero and scaled to have unit variance) for model fitting?}
  \item{admm.max.iter}{Maximum number of ADMM iterations.}
  \item{quiet}{\code{TRUE} or \code{FALSE} variable indicating whether progress should be printed.}
}

\value{
  \item{beta}{A \eqn{p \times} \code{nlambda} sparse matrix consisting of the estimates of \eqn{\beta} for the candidate values of \eqn{\lambda}. It is recommended to use \code{penAFT.coef} to extract coefficients.}
  \item{lambda}{The candidate tuning parameter values. }
  \item{standardize}{Were predictors standardized to have unit variance for model fitting? }
  \item{X.mean}{The mean of the predictors.}
  \item{X.sd}{The standard deviation of the predictors.}
  \item{alpha}{The tuning parameter \eqn{\alpha}. See documentation.}
}


\details{
Given \eqn{(\log(y_1), x_1, \delta_1),\dots,(\log(y_n), x_n, \delta_n)} where \eqn{y_i} is the minimum of the survival time and censoring time, \eqn{x_i} is a \eqn{p}-dimensional predictor, and \eqn{\delta_i} is the indicator of censoring,  \code{penAFT} fits the solution path for the argument minimizing
\deqn{\frac{1}{n^2}\sum_{i=1}^n \sum_{j=1}^n \delta_i \{ \log(y_i) - \log(y_j) - (x_i - x_j)'\beta \}^{-} + \lambda g(\beta)}
where \eqn{\{a \}^{-} := \max(-a, 0) }, \eqn{\lambda > 0}, and \eqn{g} is either the weighted elastic net penalty (\code{penalty = "EN"}) or weighted sparse group lasso penalty (\code{penalty = "SG"}).
The weighted elastic net penalty is defined as 
  \deqn{\alpha \| w \circ \beta\|_1 + \frac{(1-\alpha)}{2}\|\beta\|_2^2}
  where \eqn{w} is a set of non-negative weights (which can be specified in the \code{weight.set} argument). The weighted sparse group-lasso penalty we consider is
   \deqn{\alpha \| w \circ \beta\|_1 + (1-\alpha)\sum_{l=1}^G v_l\|\beta_{\mathcal{G}_l}\|_2}
   where again, \eqn{w} is a set of non-negative weights and \eqn{v_l} are weights applied to each of the \eqn{G} groups. 
}

\examples{
# --------------------------------------
# Generate data  
# --------------------------------------
set.seed(1)
genData <- genSurvData(n = 50, p = 100, s = 25, mag = 0.5, cens.quant = 0.6)
X <- genData$X
logY <- genData$logY
delta <- genData$status


# -----------------------------------------------
# Fit elastic net penalized estimator
# -----------------------------------------------
fit.en <- penAFT(X = X, logY = logY, delta = delta,
                   nlambda = 50, lambda.ratio.min = 0.01,
                   penalty = "EN",
                   alpha = 1)
                   

coef.en.10 <- penAFT.coef(fit.en, lambda = fit.en$lambda[10])

# ------------------------------------------------
# Fit weighted elastic net penalized estimator
# ------------------------------------------------
weight.set <- list("w" = c(0, 0, rep(1, 98)))
fit.weighted.en <- penAFT(X = X, logY = logY, delta = delta,
                   nlambda = 50, weight.set = weight.set,
                   penalty = "EN",
                   alpha = 1)
coef.wighted.en.10 <- penAFT.coef(fit.weighted.en, lambda = fit.weighted.en$lambda[10])
                   
                   
# ------------------------------------------------
# Fit ridge penalized estimator with user-specified lambda
# ------------------------------------------------
fit.ridge <- penAFT(X = X, logY = logY, delta = delta,
                   lambda = 10^seq(-4, 4, length=50), 
                   penalty = "EN",
                   alpha = 0)
                   
                   
# -----------------------------------------------
# Fit sparse group penalized estimator
# -----------------------------------------------
groups <- rep(1:5, each = 20)
fit.sg <- penAFT(X = X, logY = logY, delta = delta,
                   nlambda = 50, lambda.ratio.min = 0.01,
                   penalty = "SG", groups = groups, 
                   alpha = 0.5)
                   
# -----------------------------------------------
# Fit weighted sparse group penalized estimator
# -----------------------------------------------
groups <- rep(1:5, each = 20)
weight.set <- list("w" = c(0, 0, rep(1, 98)), 
      "v" = 1:5)
fit.weighted.sg <- penAFT(X = X, logY = logY, delta = delta,
                   nlambda = 100, 
                   weight.set = weight.set,
                   penalty = "SG", groups = groups, 
                   alpha = 0.5)

coef.weighted.sg.20 <- penAFT.coef(fit.weighted.sg, lambda = fit.weighted.sg$lambda[20])

}