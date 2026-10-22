\name{cqrfdr-package}
\alias{cqrfdr-package}
\alias{cqrfdr}
\docType{package}
\title{FDR-Controlled Variable Selection for High-Dimensional Censored Quantile Regression}
\description{
The package implements the single-split mirror-statistic procedure (SS-MS) for
selecting covariates of a conditional quantile of a right-censored response with
the false discovery rate controlled, together with the building blocks it needs.
}
\details{
Main functions:
\itemize{
  \item \code{\link{ssms}}: SS-MS selection (the proposed method).
  \item \code{\link{ds}}: single data splitting (Dai et al., 2023).
  \item \code{\link{mds}}: multiple data splitting aggregation (Dai et al., 2023).
  \item \code{\link{msms}}: SS-MS on many random splits, aggregated by inclusion
    rates (the multiple-splitting version of SS-MS).
  \item \code{\link{hdcqr}}: the \eqn{L_1}-penalized sequential Peng--Huang
    estimator (L-HDCQR), with the penalty sequence of Zheng, Peng and He (2018)
    (\code{lambda = "zheng"}) or one CV or pivotal penalty level;
    \code{\link{hdcqr_path}} is the low-level solver.
  \item \code{\link{mirror_stat}}, \code{\link{mirror_threshold}},
    \code{\link{pivotal_lambda}}: components of the procedure.
  \item \code{\link{sim_cqr}}, \code{\link{fdp_power}}: simulation designs and
    evaluation.
}
All computations are in base R (no compiled code). See
\code{vignette("cqrfdr")} for a worked example.
}
\author{Jinzhao Yu, Zili Liu (maintainer, \email{ZLiu20@mdanderson.org})}
\keyword{package}
