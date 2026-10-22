\name{ds}
\alias{ds}
\title{Single Data Splitting (DS) of Dai et al. for Censored Quantile Regression}
\description{
The single data-splitting procedure of Dai, Lin, Xing and Liu (2023)
transferred to censored quantile regression: a half-half split; on the first
half, the L1-penalized censored quantile regression path with the penalty
chosen by K-fold cross-validation (10 folds, as in their simulations); the
cross-validated support, without a cap, as the candidate set; on the second
half, the unpenalized Peng--Huang estimator on the candidate set (in place of
their least-squares refit); mirror statistics; and the cutoff
\eqn{\#\{M_j < -t\}/(\#\{M_j > t\} \vee 1) \le q} without the \eqn{+1}.
}
\usage{
ds(x, y, delta = NULL, tau = 0.5, q = c(0.1, 0.2), nu = 0.1, n_grid = 40L,
   taus = NULL, log_time = NULL, nfolds = 10L, foldid = NULL,
   loss = c("ph", "martingale", "deviance"), split = NULL,
   std = c("sd", "robust"), coef_tol = 1e-6, nsim = 500L, cores = 1L)
}
\arguments{
  \item{x, y, delta, tau, q, nu, n_grid, taus, log_time, nfolds, foldid, loss,
    split, std, coef_tol, nsim, cores}{as in \code{\link{ssms}}.}
}
\details{
\code{ds()} is \code{ssms(..., screen = "cv", screen_cap = "none", offset = 0,
nfolds = 10)}. It is provided as the comparator of SS-MS; its FDR guarantee is
asymptotic (Dai et al., 2023), and without a cap the refit fails when the CV
support is large relative to the events in the second half.
}
\value{An object of class \code{c("cqrfdr_ds", "cqrfdr_ssms", "cqrfdr")} with the
components of \code{\link{ssms}} and \code{method = "DS"}.}
\references{
Dai, C., Lin, B., Xing, X. and Liu, J. S. (2023). False discovery rate control
via data splitting. \emph{Journal of the American Statistical Association}.
}
\seealso{\code{\link{ssms}}, \code{\link{mds}}}
\examples{
\donttest{
set.seed(2)
d <- sim_cqr(n = 200, p = 60, beta = 0.8, cens = 0.3, calib_n = 20000)
fit <- ds(d$x, d$y, d$delta, tau = 0.5, n_grid = 10, nfolds = 3)
fit
}
}
\keyword{survival}
\keyword{regression}
