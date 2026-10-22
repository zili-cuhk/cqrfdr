\name{fdp_power}
\alias{fdp_power}
\title{False Discovery Proportion and Power of a Selection}
\description{Compares a selected set with the true active set.}
\usage{fdp_power(selected, active)}
\arguments{
  \item{selected}{selected covariates (indices or names).}
  \item{active}{true active covariates, of the same type.}
}
\value{Named vector: number \code{selected}, number of \code{false} and
\code{true} discoveries, \code{FDP} (0 for an empty selection) and \code{power}
(true positive proportion).}
\examples{
fdp_power(c(1, 2, 3, 15), active = 1:10)
}
