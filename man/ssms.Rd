\name{ssms}
\alias{ssms}
\title{Single-Split Mirror-Statistic Selection (SS-MS) for Censored Quantile Regression}
\description{
Selects covariates of the conditional \eqn{\tau}-th quantile of a
right-censored response with the false discovery rate controlled at level
\code{q}. The data are split at random into two halves. On the first half the
sequential \eqn{L_1}-penalized Peng--Huang estimator (L-HDCQR) is fitted over
a grid of quantile levels and a screened set of candidates is formed; on the
second half the screened model is refitted without penalty. The two estimates
at \eqn{\tau} are combined into mirror statistics
\eqn{M_j = \mathrm{sign}(\hat b^{(1)}_j \hat b^{(2)}_j)(|\hat b^{(1)}_j| + |\hat b^{(2)}_j|)},
and covariates with \eqn{M_j \ge \hat t} are selected, where
\deqn{\hat t = \min\Bigl\{t > 0:\ \frac{\mathrm{offset} + \#\{j: M_j \le -t\}}{\#\{j: M_j \ge t\} \vee 1} \le q\Bigr\}.}
}
\usage{
ssms(x, y, delta = NULL, tau = 0.5, q = c(0.1, 0.2),
     nu = 0.1, n_grid = 40L, taus = NULL, log_time = NULL,
     lambda = c("cv", "pivotal"), mult_grid = NULL, fixed_mult = 0.5,
     nfolds = 5L, foldid = NULL, loss = c("ph", "martingale", "deviance"),
     screen = c("path", "cv", "size"), screen_cap = "d2",
     offset = 1, split = NULL, std = c("sd", "robust"),
     coef_tol = 1e-6, nsim = 500L, cores = 1L)
}
\arguments{
  \item{x}{numeric \eqn{n \times p}{n x p} covariate matrix (no intercept column).
    Column names, if present, are used to label the selection.}
  \item{y}{observed response \eqn{Y = \min(T, C)} on the scale of the quantile
    model (usually log survival time), or a right-censored
    \code{survival::Surv} object, in which case the times are log-transformed
    unless \code{log_time = FALSE}.}
  \item{delta}{event indicator (1 = event, 0 = censored); ignored when \code{y}
    is a \code{Surv} object.}
  \item{tau}{target quantile level.}
  \item{q}{one or more nominal FDR levels. Only the threshold depends on
    \code{q}, so all levels are obtained from one fit.}
  \item{nu, n_grid}{the quantile grid is \code{seq(nu, tau, length.out = n_grid)};
    \code{nu} is the lower end \eqn{\nu} of the grid.}
  \item{taus}{optional increasing grid of quantile levels ending at the target;
    overrides \code{tau}, \code{nu} and \code{n_grid}.}
  \item{log_time}{if \code{TRUE}, \code{y} is log-transformed. Defaults to
    \code{TRUE} for a \code{Surv} object and \code{FALSE} otherwise.}
  \item{lambda}{penalty rule for the first-stage fit: \code{"cv"} (default)
    chooses a multiplier of the pivotal penalty by \code{nfolds}-fold
    cross-validation; \code{"pivotal"} uses \code{fixed_mult} times the pivotal
    penalty; a positive number is used directly as the penalty at \code{tau}
    (standardized covariates). See Details.}
  \item{mult_grid}{multipliers of the pivotal penalty searched by CV. Default
    \code{c(0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1, 1.25)}.}
  \item{fixed_mult}{multiplier used when \code{lambda = "pivotal"}.}
  \item{nfolds}{number of CV folds (stratified by the event indicator).}
  \item{foldid}{optional fold labels for the full sample (length \eqn{n});
    SS-MS uses them restricted to the first half.}
  \item{loss}{CV criterion: \code{"ph"} (held-out Peng--Huang objective,
    default), \code{"martingale"} or \code{"deviance"}.}
  \item{screen}{how the screened set is formed; see Details.}
  \item{screen_cap}{maximum number of screened candidates: \code{"d2"} (default)
    \eqn{\lfloor d_2/\log d_2\rfloor} with \eqn{d_2} the number of events in the
    second half; \code{"d1"} the same with the events of the first half (so the
    screen depends on the first half only); a positive integer; or
    \code{"none"}.}
  \item{offset}{1 gives the threshold \eqn{U_+} with the finite-sample FDR bound;
    0 gives the threshold of Dai et al. (2023).}
  \item{split}{optional indices of the first half; random by default.}
  \item{std}{standardization of the covariates, computed on the first half:
    mean and standard deviation (\code{"sd"}) or median and MAD
    (\code{"robust"}, for heavy-tailed covariates).}
  \item{coef_tol}{coefficients with absolute value below this are treated as zero.}
  \item{nsim}{Monte Carlo draws for the pivotal penalty.}
  \item{cores}{number of cores for the CV fits (\code{parallel::mclapply};
    ignored on Windows).}
}
\details{
\strong{Estimator.} At grid point \eqn{\tau_k} the sequential estimator of
Zheng, Peng and He (2018) minimises
\deqn{\sum_i \delta_i |Y_i - Z_i^\top b| + b^\top \sum_i \delta_i Z_i - 2 b^\top \sum_i a_{ik} Z_i + \lambda_k \sum_{j \ge 1} |b_j|,}
with \eqn{a_{ik} = H(\nu) + \sum_{r<k} 1\{Y_i \ge Z_i^\top \hat b(\tau_r)\}\{H(\tau_{r+1}) - H(\tau_r)\}}
and \eqn{H(u) = -\log(1-u)}. Each step is one LAD-lasso linear program, solved by
an interior-point method (see \code{\link{hdcqr_path}}).

\strong{Penalty.} The pivotal level is
\eqn{\lambda_0 = 2 \cdot 1.1 \cdot Q_{0.9}(\max_j |\sum_i Z_{ij}\xi_i|)} with
\eqn{\xi_i = 1\{U_i \le \tau\} - H(U_i \wedge \tau)}, \eqn{U_i} uniform (see
\code{\link{pivotal_lambda}}); on the grid
\eqn{\lambda_k = c\,\lambda_0 (\tau_k/\tau)^{1/2}}, and \eqn{c} is chosen by CV.

\strong{Screen.} With \code{screen = "path"} (default, the procedure studied in
the paper) the screened set is the support of the CV fit; if it holds fewer than
the cap, the multiplier is lowered along \code{mult_grid} and candidates are
added in the order in which they enter the lasso path (ties broken by larger
\eqn{|b|}), each keeping its coefficient from the fit at which it entered, until
the cap is reached. \code{screen = "cv"} keeps the CV support only (capped by
largest \eqn{|b|}); \code{screen = "size"} uses the coefficients of the first
fit along the grid with at least \code{cap} nonzero coefficients. All rules use
the first half only, so the mirror-statistic argument is unchanged.

A nonempty selection with \code{offset = 1} contains at least
\eqn{\lceil 1/q \rceil} covariates; when fewer signals are expected, use a
larger \code{q} or \code{offset = 0}.

The procedure is random through the split, the CV folds and the Monte Carlo
pivotal penalty; call \code{set.seed} for reproducibility.
}
\value{
An object of class \code{c("cqrfdr_ssms", "cqrfdr")}, a list with components
  \item{selected}{named list (one element per \code{q}) of selected column indices,
    named by the column names of \code{x}.}
  \item{threshold}{selection threshold for each \code{q} (\code{Inf} if no threshold
    satisfies the bound).}
  \item{M, b1, b2}{mirror statistics and the first- and second-half coefficients
    at \code{tau} (standardized scale; zero outside the screen).}
  \item{screen, cap, screen_rule}{the screened candidates and the cap.}
  \item{split}{indices of the first half.}
  \item{tuning}{list with the CV multiplier (\code{mult}), the multiplier at which
    the screen stopped (\code{mult_screen}), the pivotal level and the CV curve.}
  \item{status}{\code{"ok"}, \code{"empty_screen"}, \code{"no_threshold"}, or a
    failure code (\code{"stage1_no_root"}, \code{"stage2_no_root"},
    \code{"tuning_failed"}, ...).}
plus \code{call}, \code{q}, \code{tau}, \code{taus}, \code{n}, \code{p},
\code{events}. Use \code{\link{selected}} to extract a selection and
\code{\link[=summary.cqrfdr]{summary}} for the screened covariates.
}
\references{
Peng, L. and Huang, Y. (2008). Survival analysis with quantile regression models.
\emph{Journal of the American Statistical Association}, 103, 637--649.

Zheng, Q., Peng, L. and He, X. (2018). High dimensional censored quantile
regression. \emph{The Annals of Statistics}, 46, 308--343.

Dai, C., Lin, B., Xing, X. and Liu, J. S. (2023). False discovery rate control
via data splitting. \emph{Journal of the American Statistical Association}.
}
\seealso{\code{\link{mds}}, \code{\link{hdcqr}}, \code{\link{sim_cqr}}}
\examples{
set.seed(1)
d <- sim_cqr(n = 200, p = 60, s = 10, beta = 0.8, cens = 0.3, calib_n = 20000)
fit <- ssms(d$x, d$y, d$delta, tau = 0.5, q = c(0.1, 0.2),
            n_grid = 10, nfolds = 3)
fit
selected(fit, q = 0.2)
fdp_power(selected(fit, q = 0.2), paste0("X", d$active))
head(summary(fit))
plot(fit)
}
\keyword{survival}
\keyword{regression}
