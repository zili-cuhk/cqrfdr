\name{sim_cqr}
\alias{sim_cqr}
\title{Simulate Censored Data from a Linear Quantile Model}
\description{
Generates data from the designs of the simulation study:
\itemize{
  \item \code{model = "location"}:
    \eqn{\log T = \alpha + Z^\top b^* + \sigma_e\{E - F^{-1}(\tau)\}}, where the
    first \code{s} coefficients are
    \eqn{b^*_j = \beta\,(\pm 1) + \beta_{\mathrm{slope}}(\tau - 0.5)(\pm 1)}
    with alternating signs, so \eqn{b^*} is the coefficient of the \eqn{\tau}-th
    conditional quantile.
  \item \code{model = "hetero"}:
    \eqn{\log T = \alpha + Z^\top b + \sigma(1 + \gamma Z_h) E} with
    \eqn{Z_h = \Phi(W_h)}, \eqn{h = s + 1}; the coefficient of \eqn{Z_h} is
    \eqn{\gamma\sigma F^{-1}(u)}, so the active set depends on the quantile level.
}
Covariates have Toeplitz correlation \eqn{\rho^{|j-k|}} (Gaussian or
multivariate \eqn{t}); \eqn{\log C \sim N(c_0, \mathrm{censor\_sd}^2)}
independent of \eqn{(T, Z)}, with \eqn{c_0} calibrated by Monte Carlo to the
target censoring rate.
}
\usage{
sim_cqr(n = 400L, p = 200L, s = 10L, beta = 0.8, beta_slope = beta,
        tau = 0.5, rho = 0.5, cens = 0.3,
        xdist = c("gaussian", "t"), edist = c("normal", "t", "exp"),
        model = c("location", "hetero"), xdf = 5, edf = 5,
        err_sd = 0.5, alpha = 0.5, gamma = 3, sigma = 0.5,
        censor_sd = 1, c0 = NULL, calib_n = 200000L)
}
\arguments{
  \item{n, p, s}{sample size, number of covariates, number of signals.}
  \item{beta, beta_slope}{signal level and its slope in \eqn{\tau} (location model).}
  \item{tau}{target quantile level.}
  \item{model}{\code{"location"} or \code{"hetero"}; see Description.}
  \item{rho}{AR(1) correlation of the covariates.}
  \item{cens}{target censoring rate (0 for no censoring).}
  \item{xdist, xdf}{covariate distribution and degrees of freedom of the
    multivariate \eqn{t} (rescaled to unit variance when \code{xdf > 2}).}
  \item{edist, edf}{error distribution (unit variance) and \eqn{t} degrees of freedom.}
  \item{err_sd}{error scale \eqn{\sigma_e} in the location model.}
  \item{alpha, gamma, sigma}{intercept and heteroscedasticity parameters.}
  \item{censor_sd}{standard deviation of \eqn{\log C}.}
  \item{c0}{censoring location; calibrated if \code{NULL}. Pass the value from a
    previous call to generate replicates without recalibration.}
  \item{calib_n}{Monte Carlo sample size of the calibration.}
}
\value{
A list with \code{x} (\eqn{n \times p}{n x p}, columns \code{X1}, ...), \code{y}
(observed log time), \code{delta}, \code{time} (\code{exp(y)}), \code{active}
(indices of the nonzero coefficients at \code{tau}), \code{beta} (true
coefficients at \code{tau}, intercept first), the realised \code{cens_rate},
\code{c0}, the \code{calibration} diagnostics (including the probability of
remaining at risk at the \eqn{\tau}-th quantile) and the full \code{spec}.
}
\examples{
set.seed(6)
d <- sim_cqr(n = 300, p = 100, cens = 0.3, calib_n = 50000)
d$cens_rate
d$active
## replicate with the same censoring distribution
d2 <- sim_cqr(n = 300, p = 100, c0 = d$c0)
}
