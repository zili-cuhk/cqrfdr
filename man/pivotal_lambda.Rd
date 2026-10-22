\name{pivotal_lambda}
\alias{pivotal_lambda}
\title{Pivotal Penalty Level for Censored Quantile Regression}
\description{
Monte Carlo pivotal penalty
\eqn{\lambda_0 = 2 C\, Q_{\mathrm{level}}(\max_j |\sum_i Z_{ij}\xi_i|)}, with
\eqn{\xi_i = 1\{U_i \le \tau\} - H(U_i \wedge \tau)}, \eqn{U_i \sim} Uniform(0,1)
and \eqn{H(u) = -\log(1-u)}: the score of the Peng--Huang objective at the truth
without censoring. Under censoring it is conservative.
}
\usage{
pivotal_lambda(Zs, tau, nsim = 500L, level = 0.9, C = 1.1)
}
\arguments{
  \item{Zs}{standardized covariate matrix without the intercept column.}
  \item{tau}{target quantile level.}
  \item{nsim}{number of Monte Carlo draws.}
  \item{level}{quantile of the maximum score.}
  \item{C}{inflation constant.}
}
\value{The penalty level at \code{tau}. On a grid the penalties are
\eqn{\lambda_0(\tau_k/\tau)^{1/2}}.}
\examples{
set.seed(5)
pivotal_lambda(scale(matrix(rnorm(200 * 50), 200)), tau = 0.5)
}
