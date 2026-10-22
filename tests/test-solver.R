## Solver checks (base R; quantreg comparisons run only if it is installed)
library(cqrfdr)
lad_lasso <- cqrfdr:::lad_lasso
lad_lasso_ws <- cqrfdr:::lad_lasso_ws
set.seed(1)
obj <- function(A, y, lam, b) sum(abs(y - A %*% b)) + sum(lam * abs(b))

## 1. working-set solver = full solver (same LP)
n <- 60; p <- 150
A <- cbind(1, matrix(rnorm(n * (p - 1)), n)); y <- A[, 2] - A[, 3] + rnorm(n)
lam <- c(0, rep(4, p - 1))
o1 <- obj(A, y, lam, lad_lasso(A, y, lam)$coef)
o2 <- obj(A, y, lam, lad_lasso_ws(A, y, lam)$coef)
stopifnot(abs(o1 - o2) / o1 < 1e-6)

## 2. agreement with quantreg::rq.fit.lasso and rq.fit.fnb
if (requireNamespace("quantreg", quietly = TRUE)) {
  for (cfg in list(c(60, 20), c(50, 120))) {
    n <- cfg[1]; p <- cfg[2]
    A <- cbind(1, matrix(rnorm(n * (p - 1)), n)); y <- A[, 2] - A[, 3] + rnorm(n)
    lam <- c(0, rep(2, p - 1))
    b1 <- lad_lasso(A, y, lam)$coef
    b2 <- quantreg::rq.fit.lasso(A, y, tau = 0.5, lambda = lam)$coefficients
    stopifnot(abs(obj(A, y, lam, b1) - obj(A, y, lam, b2)) / obj(A, y, lam, b2) < 1e-6)
  }
  b0 <- lad_lasso(A[, 1:10], y, 0)$coef
  stopifnot(max(abs(b0 - quantreg::rq.fit.fnb(A[, 1:10], y, 0.5)$coefficients)) < 1e-4)

  ## 3. without censoring the Peng-Huang path recovers quantile regression
  n <- 1500; X <- matrix(rnorm(n * 3), n)
  y <- drop(1 + X %*% c(1, -1, 0.5) + (1 + 0.5 * pnorm(X[, 1])) * rnorm(n))
  taus <- seq(0.02, 0.6, length.out = 59)
  f <- hdcqr_path(y, rep(1, n), cbind(1, X), taus, 0)
  stopifnot(f$status == "ok",
            max(abs(f$B[, 59] - quantreg::rq.fit.fnb(cbind(1, X), y, 0.6)$coefficients)) < 0.1)
}

## 4. too many columns for the number of events is reported, not hidden
set.seed(2)
d <- sim_cqr(n = 200, p = 150, cens = 0.5, calib_n = 20000)
f <- hdcqr_path(d$y, d$delta, cbind(1, scale(d$x[, 1:140])), seq(0.1, 0.5, length.out = 10), 0)
stopifnot(f$status != "ok")

## 5. pivotal penalty scales like sqrt(n log p)
set.seed(3)
l1 <- pivotal_lambda(scale(matrix(rnorm(200 * 100), 200)), 0.5)
l2 <- pivotal_lambda(scale(matrix(rnorm(800 * 100), 800)), 0.5)
stopifnot(l2 / l1 > 1.6, l2 / l1 < 2.4)
cat("solver tests passed\n")
