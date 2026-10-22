\name{mirror_stat}
\alias{mirror_stat}
\alias{mirror_threshold}
\title{Mirror Statistics and the Data-Driven Threshold}
\description{
\code{mirror_stat} combines two independent estimates of the same coefficients,
\eqn{M_j = \mathrm{sign}(b_{1j} b_{2j})(|b_{1j}| + |b_{2j}|)}.
\code{mirror_threshold} returns the smallest \eqn{t > 0} among the observed
\eqn{|M_j|} with
\eqn{(\mathrm{offset} + \#\{M_j \le -t\}) / (\#\{M_j \ge t\} \vee 1) \le q}.
}
\usage{
mirror_stat(b1, b2)
mirror_threshold(M, q, offset = 1)
}
\arguments{
  \item{b1, b2}{numeric vectors of the same length.}
  \item{M}{mirror statistics (zeros are ignored).}
  \item{q}{nominal FDR level, one number in (0, 1).}
  \item{offset}{1 for \eqn{U_+} (finite-sample bound), 0 for the threshold of
    Dai et al. (2023).}
}
\value{\code{mirror_stat}: vector of mirror statistics.
\code{mirror_threshold}: the threshold, or \code{Inf} if no \eqn{t} satisfies
the bound (nothing is selected).}
\examples{
M <- mirror_stat(c(2, 1.5, -1, 0.3, 0.2), c(1.8, 1.2, 0.4, -0.1, 0.25))
M
mirror_threshold(M, q = 0.5, offset = 0)
}
