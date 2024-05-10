\name{penAFT-package}
\alias{penAFT-package}
\docType{package}
\title{Fit and tune the semiparameteric accelerated failure time model by minimizing a rank-based estimation criterion plus weighted elastic net or weighted sparse group-lasso penalty.}
\description{
This package contains numerous functions related t a penalized rank-based estimator of the semiparameteric AFT, the so-called penalized Gehan estimator. In particular, the main functions are for solution path computation, cross-validation, prediction, and coefficient extraction. 
}
\details{The primary functions are \code{penAFT} and \code{penAFT.cv}, the latter of which performs cross-validation. In general, both functions fit the penalized Gehan estimator. Given \eqn{(\log(y_1), x_1, \delta_1),\dots,(\log(y_n), x_n, \delta_n)} where \eqn{y_i} is the minimum of the survival time and censoring time, \eqn{x_i} is a \eqn{p}-dimensional predictor, and \eqn{\delta_i} is the indicator of censoring,  \code{penAFT} fits the solution path for the argument minimizing
\deqn{\frac{1}{n^2}\sum_{i=1}^n \sum_{j=1}^n \delta_i \{ \log(y_i) - \log(y_j) - (x_i - x_j)'\beta \}^{-} + \lambda g(\beta)}
where \eqn{\{a \}^{-} := \max(-a, 0) }, \eqn{\lambda > 0}, and \eqn{g} is either the weighted elastic net penalty or weighted sparse group lasso penalty. The weighted elastic net penalty is defined as 
  \deqn{\alpha \| w \circ \beta\|_1 + \frac{(1-\alpha)}{2}\|\beta\|_2^2}
  where \eqn{w} is a set of non-negative weights (which can be specified in the \code{weight.set} argument). The weighted sparse group-lasso penalty we consider is
   \deqn{\alpha \| w \circ \beta\|_1 + (1-\alpha)\sum_{l=1}^G v_l\|\beta_{\mathcal{G}_l}\|_2}
   where again, \eqn{w} is a set of non-negative weights and \eqn{v_l} are weights applied to each of the \eqn{G} (user-specified) groups. 
   
   For a comprehensive description of the algorithm, and more details about rank-based estimation in general, please refer to the referenced manuscript. 
}
\author{Aaron J. Molstad and Piotr M. Suder
Maintainer: Aaron J. Molstad <amolstad@ufl.edu>
}


\keyword{package}