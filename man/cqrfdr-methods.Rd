\name{cqrfdr-methods}
\alias{selected}
\alias{selected.cqrfdr}
\alias{selected.hdcqr}
\alias{print.cqrfdr}
\alias{summary.cqrfdr}
\alias{plot.cqrfdr}
\title{Methods for SS-MS and MDS Fits}
\description{
Extract the selection, print, summarize and plot fits from \code{\link{ssms}}
and \code{\link{mds}}.
}
\usage{
selected(object, ...)
\method{selected}{cqrfdr}(object, q = object$q[1L], names = TRUE, ...)
\method{selected}{hdcqr}(object, names = TRUE, ...)
\method{print}{cqrfdr}(x, ...)
\method{summary}{cqrfdr}(object, ...)
\method{plot}{cqrfdr}(x, q = x$q, label = TRUE, ...)
}
\arguments{
  \item{object, x}{a fit from \code{ssms}, \code{mds} or \code{hdcqr}.}
  \item{q}{FDR level(s), among those used in the fit. For the MDS plot only the
    first is used.}
  \item{names}{return column names (default) rather than indices.}
  \item{label}{label the ten largest statistics in the plot.}
  \item{...}{passed to \code{plot} (for \code{plot}); otherwise unused.}
}
\value{
\code{selected}: character vector of names (or integer indices).
\code{summary}: for SS-MS, a data frame of the screened covariates with
\code{b1}, \code{b2}, \code{M} and one selection indicator per \code{q}, sorted
by \code{M}; for MDS, the covariates with positive inclusion rate, their rates
and selection indicators. \code{plot}: the ordered mirror statistics with the
thresholds (SS-MS) or the ordered inclusion rates (MDS), drawn with base
graphics.
}
\seealso{\code{\link{ssms}}, \code{\link{mds}}}
