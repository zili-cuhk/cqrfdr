\name{hdcqr_path}
\alias{hdcqr_path}
\title{Sequential Peng--Huang Path with a Given Penalty}
\description{
Low-level solver: computes the (penalized) sequential Peng--Huang estimator on a
grid of quantile levels for a given design and penalty sequence.
}
\usage{
hdcqr_path(y, delta, Z, taus, lambda = 0, R0 = 1e4, R_growth = 10,
           R_tries = 5L, eps = 1e-6, working_set = ncol(Z) > 100L)
}
\arguments{
  \item{y}{observed response (e.g. log time).}
  \item{delta}{event indicator.}
  \item{Z}{design matrix whose \strong{first column is the intercept} (not
    penalized). Covariates should be standardized.}
  \item{taus}{increasing grid of quantile levels; the first element is
    \eqn{\nu}.}
  \item{lambda}{penalty at each grid point (recycled); 0 gives the unpenalized
    Peng--Huang estimator.}
  \item{R0, R_growth, R_tries}{the two linear terms of the objective are written
    as absolute-value rows with a large pseudo-response \eqn{R}; \eqn{R} starts at
    \code{R0} and is multiplied by \code{R_growth} (at most \code{R_tries} times)
    until the artificial residuals are positive.}
  \item{eps}{duality-gap tolerance of the interior-point solver.}
  \item{working_set}{use an exact working-set (active-set with KKT check) scheme
    for wide designs.}
}
\details{
Each grid point is one LAD-lasso linear program
\eqn{\min_b \sum |y - A b| + \sum_j \lambda_j |b_j|}, solved by an R port of the
Frisch--Newton interior-point method of \pkg{quantreg} (\code{rq.fit.lasso}).
If the estimating equation has no root at some \eqn{\tau_k} (\eqn{\tau_k} beyond
the identifiable range, or too many columns for the number of events), the path
stops with status \code{"no_root"}.
}
\value{
A list with \code{B}, the \code{ncol(Z)} \eqn{\times}{x} \code{length(taus)}
coefficient matrix (columns after a failure are \code{NA}), \code{status}
(\code{"ok"}, \code{"no_root"} or \code{"solver_error"}), the index \code{k} of
the last grid point attempted, and a \code{message} on failure.
}
\seealso{\code{\link{hdcqr}}, \code{\link{pivotal_lambda}}}
\examples{
set.seed(4)
d <- sim_cqr(n = 300, p = 5, s = 2, cens = 0.2, calib_n = 20000)
Z <- cbind(1, scale(d$x))
f <- hdcqr_path(d$y, d$delta, Z, taus = seq(0.1, 0.5, length.out = 20))
f$status
round(f$B[, 20], 3)
}
\keyword{survival}
