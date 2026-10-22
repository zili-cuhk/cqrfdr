\name{mds}
\alias{mds}
\title{Multiple Data Splitting (MDS) for Censored Quantile Regression}
\description{
Runs the single-split procedure of \code{\link{ssms}} on \code{B} independent
random splits and aggregates the selections through inclusion rates, following
Dai, Lin, Xing and Liu (2023). For split \eqn{b} with selection
\eqn{\hat S_b}, covariate \eqn{j} receives
\eqn{I_j = B^{-1}\sum_b 1\{j \in \hat S_b\}/(|\hat S_b| \vee 1)}; the rates are
sorted increasingly, the largest \eqn{l} with cumulative sum \eqn{\le q} is
found, and the covariates outside the \eqn{l} smallest rates (and with
\eqn{I_j > 0}) are selected.
}
\usage{
mds(x, y, delta = NULL, tau = 0.5, q = c(0.1, 0.2), B = 50L,
    nu = 0.1, n_grid = 40L, taus = NULL, log_time = NULL,
    tune = c("cv_each", "cv_once", "pivotal"), ties = c("dai", "random"),
    fixed_mult = 0.5, mult_grid = NULL, nfolds = 10L, foldid = NULL,
    loss = c("ph", "martingale", "deviance"), screen_cap = "none",
    offset = 0, std = c("sd", "robust"), coef_tol = 1e-6, nsim = 500L,
    cores = 1L, verbose = FALSE)
}
\arguments{
  \item{x, y, delta, tau, q, nu, n_grid, taus, log_time, mult_grid, nfolds, foldid,
    loss, screen_cap, std, coef_tol, nsim, cores}{as in \code{\link{ssms}};
    \code{nfolds} defaults to 10 and \code{screen_cap} to \code{"none"}, as in
    Dai et al. (2023).}
  \item{B}{number of random splits.}
  \item{tune}{\code{"cv_each"} (default): every split is a complete DS fit with its
    own \code{nfolds}-fold CV on its first half, as MDS repeats DS in Dai et al.
    (2023). \code{"cv_once"}: the multiplier is chosen once by CV on one random
    half and reused in all splits. \code{"pivotal"}: \code{fixed_mult} without CV.}
  \item{ties}{\code{"dai"} (default) is the rule of Dai et al. (2023): every
    covariate with \eqn{I_j > I_{(l)}} is selected; when all signals have tied
    rates this can remove the whole block of signals at once. \code{"random"}
    removes exactly the \eqn{l} smallest rates, breaking ties at random.}
  \item{fixed_mult}{multiplier of the pivotal penalty when \code{tune = "pivotal"}.}
  \item{offset}{threshold offset within each split; 0 as in Dai et al. (2023).}
  \item{verbose}{print the status of every split.}
}
\details{
The defaults reproduce the design of Dai et al. (2023): each split is DS
(\code{\link{ds}}), i.e. the CV support (10-fold CV) without a cap as the screen,
the unpenalized refit on the second half and the cutoff \eqn{U}. Splits that fail numerically
are dropped and counted. The computing time is roughly \code{B} times that of a
single SS-MS fit without CV.
}
\value{
An object of class \code{c("cqrfdr_mds", "cqrfdr")}: a list with
\code{selected} (named list, one element per \code{q}), \code{inclusion}
(\eqn{p \times}{p x} \code{length(q)} matrix of inclusion rates),
\code{n_failed_splits}, \code{split_status}, the multiplier \code{mult}, and the
settings.
}
\references{
Dai, C., Lin, B., Xing, X. and Liu, J. S. (2023). False discovery rate control
via data splitting. \emph{Journal of the American Statistical Association}.
}
\seealso{\code{\link{ds}}, \code{\link{ssms}}, \code{\link{msms}}}
\examples{
\donttest{
set.seed(2)
d <- sim_cqr(n = 200, p = 60, beta = 0.8, cens = 0.3, calib_n = 20000)
fit <- mds(d$x, d$y, d$delta, tau = 0.5, B = 5, n_grid = 10, nfolds = 3)
fit
head(summary(fit))
}
}
\keyword{survival}
\keyword{regression}
