\name{hdcqr}
\alias{hdcqr}
\alias{coef.hdcqr}
\alias{predict.hdcqr}
\alias{print.hdcqr}
\title{Penalized Sequential Censored Quantile Regression (L-HDCQR)}
\description{
Fits the \eqn{L_1}-penalized sequential Peng--Huang estimator (L-HDCQR) of
Zheng, Peng and He (2018) over a grid of quantile levels on the full sample.
With \code{lambda = "zheng"} the penalty follows their Section 5.1: a pivotal
level at the first grid point, increased at every grid point by an increment
chosen by cross-validation. With \code{"cv"} or \code{"pivotal"} one penalty
level, scaled along the grid, is chosen by cross-validation of its multiplier or
by the pivotal rule. The nonzero coefficients at the target quantile are the
L-HDCQR selection (no FDR control).
}
\usage{
hdcqr(x, y, delta = NULL, tau = 0.5, nu = 0.1, n_grid = 40L, taus = NULL,
      log_time = NULL, lambda = c("cv", "pivotal", "zheng"), mult_grid = NULL,
      fixed_mult = 0.5, nfolds = 5L, foldid = NULL,
      loss = c("ph", "martingale", "deviance"), std = c("sd", "robust"),
      coef_tol = 1e-6, nsim = 500L, cores = 1L,
      increments = c(0, exp(seq(log(0.02), log(20), length.out = 12))))

\method{coef}{hdcqr}(object, tau = object$tau,
     scale = c("original", "standardized"), ...)

\method{predict}{hdcqr}(object, newx, tau = object$tau, ...)

\method{print}{hdcqr}(x, ...)
}
\arguments{
  \item{x}{covariate matrix (for \code{hdcqr}) or a fitted object (for \code{print}).}
  \item{y, delta, nu, n_grid, taus, log_time, mult_grid, fixed_mult, nfolds,
    foldid, loss, std, coef_tol, nsim, cores}{as in \code{\link{ssms}}.}
  \item{lambda}{\code{"cv"} (default) or \code{"pivotal"}: one penalty level,
    scaled along the grid, with the multiplier chosen by \code{nfolds}-fold CV or
    fixed at \code{fixed_mult}; \code{"zheng"}: the sequential penalty of Zheng,
    Peng and He (2018, Section 5.1), used for L-HDCQR in the paper; or one
    positive number, the penalty at \code{tau}.}
  \item{increments}{with \code{lambda = "zheng"}: candidate increments, in units
    of \eqn{\sqrt{n\log(p \vee n)}} times the grid spacing, from which the
    increment at each grid point is chosen by \code{nfolds}-fold CV of the
    held-out Peng--Huang objective. Zheng et al. do not report their constants;
    the default spans them on the log scale and includes 0.}
  \item{tau}{target quantile level; in \code{coef} and \code{predict}, one or more
    levels, each mapped to the nearest grid point (levels above the target are not
    estimated).}
  \item{object}{an \code{hdcqr} fit.}
  \item{scale}{coefficients on the original covariate scale (default) or on the
    standardized scale used in the fit.}
  \item{newx}{matrix of new covariate values.}
  \item{...}{unused.}
}
\value{
\code{hdcqr} returns an object of class \code{"hdcqr"} with the standardized
coefficient path \code{coef_std} (\eqn{(p+1) \times}{(p+1) x} grid), the
selection \code{selected}, the penalties \code{lambda} on the grid, the chosen
multiplier \code{mult} (or, with \code{lambda = "zheng"}, the chosen
\code{increment} at each grid point), the pivotal level, the CV curve and
\code{tuning}. \code{coef} returns
a vector (one \code{tau}) or matrix; \code{predict} returns predicted
conditional quantiles of the response (on the scale of \code{y}).
}
\references{
Zheng, Q., Peng, L. and He, X. (2018). High dimensional censored quantile
regression. \emph{The Annals of Statistics}, 46, 308--343.
}
\seealso{\code{\link{hdcqr_path}}, \code{\link{ssms}}}
\examples{
set.seed(3)
d <- sim_cqr(n = 200, p = 40, beta = 0.8, cens = 0.3, calib_n = 20000)
fit <- hdcqr(d$x, d$y, d$delta, tau = 0.5, n_grid = 10, nfolds = 3)
fit
round(coef(fit)[1:12], 3)
head(predict(fit, d$x[1:5, ], tau = c(0.25, 0.5)))
\donttest{
## the penalty sequence of Zheng, Peng and He (2018)
fz <- hdcqr(d$x, d$y, d$delta, tau = 0.5, n_grid = 10, nfolds = 3, lambda = "zheng")
fz
}
}
\keyword{survival}
\keyword{regression}
