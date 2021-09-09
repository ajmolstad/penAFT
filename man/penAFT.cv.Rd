\name{penAFT.cv}
\alias{penAFT.cv}
\title{Cross-validation function for fitting regularized Gehan estimator}
\description{A function to perform cross-validation and compute the solution path for the regularized Gehan estimator.}
\usage{
penAFT.cv(X, logY, delta, nlambda = 50, 
  lambda.ratio.min = 0.1, lambda = NULL, 
  penalty = NULL, alpha = 1,weight.set = NULL, 
  groups = NULL, tol.abs = 1e-10, tol.rel = 5e-4, 
  standardize = TRUE, nfolds = 5, cv.index = NULL, 
  admm.max.iter = 1e4,quiet = TRUE)
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
  \item{standardize}{Should predictors be standardized (i.e., scaled to have unit variance) for model fitting?}
  \item{nfolds}{The number of folds to be used for cross-validation. Default is five.}
  \item{cv.index}{A list of length \code{nfolds} of indices to be used for cross-validation. This is to be used if trying to perform cross-validation for both \eqn{\alpha} and \eqn{\lambda}. Use with extreme caution: this overwrites \code{nfolds}. }
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
Given \eqn{(\log(y_1), x_1, \delta_1),\dots,(\log(y_n), x_n, \delta_n)} where for subject \eqn{i} (\eqn{i = 1, \dots, n}), \eqn{y_i} is the minimum of the survival time and censoring time, \eqn{x_i} is a \eqn{p}-dimensional predictor, and \eqn{\delta_i} is the indicator of censoring,  \code{penAFT.cv} performs \code{nfolds} cross-validation for selecting the tuning parameter to be used in the argument minimizing
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
genData <- genSurvData(n = 50, p = 100, s = 25, cens.quant = 0.6)
X <- genData$X
logY <- genData$logY
delta <- genData$status


# -----------------------------------------------
# Fit elastic net penalized estimator
# -----------------------------------------------
fit.en <- penAFT.cv(X = X, logY = logY, delta = delta,
                   nlambda = 50, lambda.ratio.min = 0.1,
                   penalty = "EN", nfolds = 10,
                   alpha = 1)
# coefficients at tuning parameter minimizing cross-valdiation error
coef.en <- penAFT.coef(fit.en)

# predict at 30th tuning parameter from full fit
Xnew <- matrix(rnorm(10*100), nrow=10)
predict.en <- penAFT.predict(fit.en, Xnew = Xnew, lambda = fit.en$full.fit$lambda[30])


# -----------------------------------------------
# Fit sparse group penalized estimator
# -----------------------------------------------
groups <- rep(1:5, each = 20)
fit.sg <- penAFT.cv(X = X, logY = logY, delta = delta,
                   nlambda = 50, lambda.ratio.min = 0.01,
                   penalty = "SG", groups = groups, nfolds = 5,
                   alpha = 0.5)
                   
# -----------------------------------------------
# Pass fold indices
# -----------------------------------------------
groups <- rep(1:5, each = 20)
cv.index <- list()
for(k in 1:5){
  cv.index[[k]] <- which(rep(1:5, length=50) == k)
}
fit.sg.cvIndex <- penAFT.cv(X = X, logY = logY, delta = delta,
                   nlambda = 50, lambda.ratio.min = 0.01,
                   penalty = "SG", groups = groups, 
                   cv.index = cv.index,
                   alpha = 0.5)

# fit.sg.cvIndex$cv.index  == cv.index
}
  