\name{msms}
\alias{msms}
\title{Multiple Splits of SS-MS (MS-MS)}
\description{
The multiple-splitting version of \code{\link{ssms}}, constructed as Dai, Lin,
Xing and Liu (2023) construct MDS from their single data-splitting procedure
DS. SS-MS (path-ordered screen with the cap \eqn{\lfloor d_2/\log d_2\rfloor}
and the cutoff \eqn{U_+}) is run on \code{B} independent random splits, and the
selections are aggregated through inclusion rates: for split \eqn{b} with
selection \eqn{\hat S_b}, covariate \eqn{j} receives
\eqn{I_j = B^{-1}\sum_b 1\{j \in \hat S_b\}/(|\hat S_b| \vee 1)}; the rates are
sorted increasingly, the largest \eqn{l} with cumulative sum \eqn{\le q} is
found, and the covariates outside the \eqn{l} smallest rates (and with
\eqn{I_j > 0}) are selected. The aggregation removes the dependence of the
selection on a single random split.
}
\usage{
msms(x, y, delta = NULL, tau = 0.5, q = c(0.1, 0.2), B = 50L,
     nu = 0.1, n_grid = 40L, taus = NULL, log_time = NULL,
     tune = c("cv_once", "cv_each", "pivotal"), ties = c("random", "dai"),
     fixed_mult = 0.5, mult_grid = NULL, nfolds = 5L, foldid = NULL,
     loss = c("ph", "martingale", "deviance"),
     screen = c("path", "cv", "size"), screen_cap = "d2",
     offset = 1, std = c("sd", "robust"), coef_tol = 1e-6, nsim = 500L,
     cores = 1L, verbose = FALSE)
}
\arguments{
  \item{x, y, delta, tau, q, nu, n_grid, taus, log_time, mult_grid, nfolds, foldid,
    loss, screen, screen_cap, std, coef_tol, nsim, cores}{as in \code{\link{ssms}};
    \code{screen} and \code{screen_cap} apply within every split.}
  \item{B}{number of random splits.}
  \item{tune, ties, fixed_mult, verbose}{as in \code{\link{mds}}. With
    \code{tune = "cv_once"} the penalty multiplier is chosen once by CV on one
    random half and reused in all splits; the path-ordered screen then moves
    down the multiplier grid from that value within each split. With
    \code{tune = "cv_each"} every split runs its own CV.}
  \item{offset}{threshold offset within each split; 1 gives the cutoff
    \eqn{U_+} of SS-MS.}
}
\details{
\code{msms()} differs from \code{\link{mds}} only in what is run within a
split: \code{mds()} uses the CV support and the cutoff \eqn{U}, as DS in Dai et
al. (2023); \code{msms()} uses SS-MS. The finite-sample bound of SS-MS holds for
each split; for the aggregated selection, as for MDS, the guarantee is
asymptotic (Dai et al., 2023). The result also records the share of splits with
an empty selection (\code{share_empty}) and the split-wise selection sizes
(\code{split_sizes}), which show how unstable a single split would have been.
The computing time is roughly \code{B} times that of a single SS-MS fit without
CV.
}
\value{
An object of class \code{c("cqrfdr_msms", "cqrfdr_mds", "cqrfdr")} with the
components of \code{\link{mds}}, plus \code{split_sizes}, \code{share_empty},
\code{screen_rule} and \code{offset}.
}
\references{
Dai, C., Lin, B., Xing, X. and Liu, J. S. (2023). False discovery rate control
via data splitting. \emph{Journal of the American Statistical Association}.
}
\seealso{\code{\link{ssms}}, \code{\link{mds}}}
\examples{
\donttest{
set.seed(2)
d <- sim_cqr(n = 200, p = 60, beta = 0.8, cens = 0.3, calib_n = 20000)
fit <- msms(d$x, d$y, d$delta, tau = 0.5, B = 5, n_grid = 10, nfolds = 3)
fit
head(summary(fit))
}
}
\keyword{survival}
\keyword{regression}
