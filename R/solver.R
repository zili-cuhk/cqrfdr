## ---------------------------------------------------------------------------
## L-HDCQR / Peng-Huang sequential estimator
## ---------------------------------------------------------------------------
##
## At grid point tau_k the estimator minimises (Zheng, Peng & He, 2018)
##   L_k(b) = sum_i d_i |y_i - Z_i'b| + b' sum_i d_i Z_i - 2 b' sum_i a_ik Z_i
##            + sum_j lambda_kj |b_j|,
##   a_ik  = H(tau_0) + sum_{r<k} 1{y_i >= Z_i' bhat(tau_r)} {H(tau_{r+1}) - H(tau_r)},
## where H(u) = -log(1-u).  The two linear terms are written as two artificial
## absolute-value rows with a large response R, so every step is one
## LAD-lasso problem.  At the solution the artificial residuals must be
## positive; if they are not even for very large R the estimating equation
## has no root at tau_k: either tau_k is beyond the identifiable range tau_U,
## or the design has too many columns for the number of events (status
## "no_root"); a singular normal matrix gives status "solver_error".
##
## lad_lasso() is an R port of quantreg's Frisch-Newton solver (rqfnb, the
## routine behind rq.fit.lasso) that builds the normal matrix with a level-3
## BLAS crossprod and adds the penalty rows as a diagonal instead of p dense
## rank-one updates.  It solves exactly the same LP as
## quantreg::rq.fit.lasso(x, y, tau = 0.5, lambda); tests/test-solver.R
## checks the agreement.

## cumulative hazard of the standard uniform: H(u) = -log(1 - u)
H <- function(u) -log1p(-u)

.ip_step <- function(a, da, b, db) {
  r1 <- da < 0; r2 <- db < 0
  min(c(Inf, -a[r1] / da[r1], -b[r2] / db[r2]))
}

## min_b sum|y - A b| + sum_j lam_j |b_j|   (lam_j = 0: unpenalised column)
lad_lasso <- function(A, y, lam, beta = 0.99995, eps = 1e-6, maxit = 500L) {
  A <- as.matrix(A); P <- ncol(A); nA <- nrow(A)
  if (length(lam) == 1L) lam <- rep(lam, P)
  pen <- which(lam > 0); lp <- lam[pen]; npen <- length(pen); N <- nA + npen
  iA <- seq_len(nA); iP <- nA + seq_len(npen)
  Xmul  <- function(b) c(drop(A %*% b), lp * b[pen])
  Xtmul <- function(v) {
    r <- drop(crossprod(A, v[iA]))
    if (npen) r[pen] <- r[pen] + lp * v[iP]
    r
  }
  fact <- function(d) {
    M <- crossprod(A * sqrt(d[iA]))
    if (npen) { ii <- cbind(pen, pen); M[ii] <- M[ii] + lp^2 * d[iP] }
    chol(M)
  }
  solveF <- function(R, b) backsolve(R, backsolve(R, b, transpose = TRUE))

  cvec <- -c(y, numeric(npen))
  rhs <- 0.5 * Xtmul(rep(1, N))
  x <- rep(0.5, N); u <- rep(1, N)
  R <- fact(rep(1, N))
  yv <- solveF(R, Xtmul(cvec))
  s <- cvec - Xmul(yv)
  z <- pmax(s, 0); w <- pmax(-s, 0)
  sm <- abs(s) < eps; z[sm] <- z[sm] + eps; w[sm] <- w[sm] + eps
  s <- u - x
  gap <- sum(z * x) + sum(w * s); it <- 0L
  while (gap > eps && it < maxit) {
    it <- it + 1L
    d <- 1 / (z / x + w / s)
    ds <- z - w
    dz <- d * ds
    rhs0 <- rhs - Xtmul(x) + Xtmul(dz)
    R <- fact(d)
    dy <- solveF(R, rhs0)
    ds <- Xmul(dy) - ds
    dx <- d * ds; ds <- -dx
    dz <- -z * (dx / x + 1); dw <- -w * (ds / s + 1)
    dp <- min(beta * .ip_step(x, dx, s, ds), 1)
    dd <- min(beta * .ip_step(z, dz, w, dw), 1)
    if (min(dp, dd) < 1) {
      mu <- sum(x * z) + sum(s * w)
      g <- mu + dp * sum(dx * z) + dd * sum(dz * x) + dp * dd * sum(dz * dx) +
        dp * sum(ds * w) + dd * sum(dw * s) + dp * dd * sum(ds * dw)
      mu <- mu * (g / mu)^3 / (2 * N)
      dr <- d * (mu * (1 / s - 1 / x) + dx * dz / x - ds * dw / s)
      dy <- solveF(R, rhs0 + Xtmul(dr))
      uu <- Xmul(dy)
      dxdz <- dx * dz; dsdw <- ds * dw
      dx <- d * (uu - z + w) - dr; ds <- -dx
      dz <- -z + (mu - z * dx - dxdz) / x
      dw <- -w + (mu - w * ds - dsdw) / s
      dp <- min(beta * .ip_step(x, dx, s, ds), 1)
      dd <- min(beta * .ip_step(z, dz, w, dw), 1)
    }
    x <- x + dp * dx; s <- s + dp * ds; yv <- yv + dd * dy
    z <- z + dd * dz; w <- w + dd * dw
    gap <- sum(z * x) + sum(w * s)
  }
  if (!all(is.finite(yv))) stop("non-finite LP solution")
  ## x is the dual vector a in [0,1]^N of the quantile-regression LP:
  ## X'(a - 1/2) = 0 at the optimum (used for the KKT check below)
  list(coef = -yv, dual = x[iA], iterations = it, converged = gap <= eps)
}

## Exact working-set version for sparse solutions.  Solve on a set W of
## columns (b_j = 0 elsewhere), then check optimality of the excluded columns
## with the dual vector: column j may stay at zero iff
##     | sum_i A_ij (a_i - 1/2) | <= lam_j / 2.
## Violators are added and the problem re-solved; when no column violates,
## the restricted solution satisfies the KKT conditions of the full problem
## and is therefore a global minimiser.  Columns with lam_j = 0 are always in W.
lad_lasso_ws <- function(A, y, lam, init = integer(0), max_add = 100L,
                         kkt_tol = 1e-6, max_rounds = 50L, ...) {
  P <- ncol(A)
  W <- sort(unique(c(which(lam == 0), init)))
  for (round in seq_len(max_rounds)) {
    fit <- lad_lasso(A[, W, drop = FALSE], y, lam[W], ...)
    out <- setdiff(seq_len(P), W)
    if (!length(out)) break
    g <- abs(drop(crossprod(A[, out, drop = FALSE], fit$dual - 0.5)))
    viol <- g - lam[out] / 2 * (1 + kkt_tol) - kkt_tol
    bad <- which(viol > 0)
    if (!length(bad)) break
    add <- out[bad[order(-viol[bad])][seq_len(min(length(bad), max_add))]]
    W <- sort(c(W, add))
  }
  b <- numeric(P); b[W] <- fit$coef
  list(coef = b, dual = fit$dual, working_set = W, rounds = round)
}

## LAD fits interpolate about P observations (residual exactly 0 in exact
## arithmetic, +-1e-10 from an interior-point solver).  The indicators
## 1{y >= Z'b} and 1{y <= Z'b} are evaluated with a small tolerance so that
## interpolated points are treated identically by every solver/platform.
IND_TOL <- 1e-5
at_risk <- function(y, fit) y >= fit - IND_TOL * (1 + abs(fit))
had_event_by <- function(y, fit) y <= fit + IND_TOL * (1 + abs(fit))

## Sequential path ------------------------------------------------------------
## y: log observed time; delta: event indicator; Z: n x P design whose first
## column is the (unpenalised) intercept; taus: increasing grid ending at the
## target; lambda: one penalty per grid point (0 = unpenalised Peng-Huang).
## Returns B (P x m).  On failure returns status != "ok" and the grid index.
hdcqr_path <- function(y, delta, Z, taus, lambda = 0, R0 = 1e4,
                       R_growth = 10, R_tries = 5L, eps = 1e-6,
                       working_set = ncol(Z) > 100L) {
  n <- length(y); P <- ncol(Z); m <- length(taus)
  lambda <- rep_len(lambda, m)
  ev <- delta == 1
  Ze <- Z[ev, , drop = FALSE]; ye <- y[ev]
  zminus <- -colSums(Ze)
  dH <- diff(H(taus))
  a <- rep(H(taus[1L]), n)
  B <- matrix(NA_real_, P, m, dimnames = list(colnames(Z), NULL))
  for (k in seq_len(m)) {
    if (k > 1L) a <- a + at_risk(y, drop(Z %*% B[, k - 1L])) * dH[k - 1L]
    zplus <- 2 * colSums(Z * a)
    pen <- c(0, rep(lambda[k], P - 1L))
    R <- R0; ok <- FALSE
    init <- if (k > 1L) which(B[, k - 1L] != 0) else integer(0)
    for (t in seq_len(R_tries)) {
      fit <- tryCatch(
        if (working_set && lambda[k] > 0)
          lad_lasso_ws(rbind(Ze, zminus, zplus), c(ye, R, R), pen, init = init, eps = eps)
        else lad_lasso(rbind(Ze, zminus, zplus), c(ye, R, R), pen, eps = eps),
        error = function(e) e)
      if (inherits(fit, "error"))
        return(list(B = B, status = "solver_error", k = k,
                    message = conditionMessage(fit)))
      b <- fit$coef
      art <- R - c(sum(zminus * b), sum(zplus * b))
      if (all(art > 0)) { ok <- TRUE; break }
      R <- R * R_growth
    }
    if (!ok) return(list(B = B, status = "no_root", k = k,
                         message = sprintf("no root of the estimating equation at tau = %.3f", taus[k])))
    B[, k] <- b
  }
  list(B = B, status = "ok", k = m)
}

## Pivotal penalty -------------------------------------------------------------
## At the truth the subgradient of L_k is 2 sum_i Z_i M_i(tau_k), where M_i is
## the counting-process martingale.  Without censoring M_i(tau) has the
## distribution of 1{U <= tau} - H(min(U, tau)), U ~ Unif(0,1), with variance
## tau, so
##   lambda(tau) = C * q_{level}( max_j |2 sum_i Z_ij xi_i| ),
## is a Belloni-Chernozhukov-type pivotal choice (conservative under
## censoring).  The grid-wise penalty is lambda(tau_k) = lambda(tau) sqrt(tau_k/tau).
pivotal_lambda <- function(Zs, tau, nsim = 500L, level = 0.9, C = 1.1) {
  n <- nrow(Zs)
  U <- matrix(stats::runif(n * nsim), n, nsim)
  xi <- (U <= tau) - H(pmin(U, tau))
  mx <- apply(abs(crossprod(Zs, xi)), 2L, max)
  2 * C * unname(stats::quantile(mx, level))
}
lambda_grid <- function(lam_target, taus) lam_target * sqrt(taus / taus[length(taus)])

## Cross-validation loss --------------------------------------------------------
## type = "ph" (default): held-out Peng-Huang objective.  For each grid point
##   l_k(b) = mean_i [ d_i |y_i - Z_i'b| + (d_i - 2 a_ik) Z_i'b ],
##   with a_ik the compensator built from the fitted path on the held-out
##   data.  l_k is convex and its gradient is 2 x the Peng-Huang estimating
##   function, so it is minimised at beta_0(tau_k): a proper loss.
## type = "martingale": mean_i |N_i(tau_k) - Lambda_i(tau_k)| (Zheng/Fei et al.)
## type = "deviance":   mean absolute deviance residual (the draft's choice).
## The last two are not minimised at the truth; in pilot runs they choose
## very large penalties at lower quantiles (empty screens at tau = 0.25).
cv_loss <- function(y, delta, Z, B, taus, type = c("ph", "martingale", "deviance")) {
  type <- match.arg(type)
  fit <- Z %*% B; n <- length(y); m <- length(taus)
  comp <- rep(H(taus[1L]), n); dH <- diff(H(taus)); L <- numeric(m)
  for (k in seq_len(m)) {
    if (k > 1L) comp <- comp + at_risk(y, fit[, k - 1L]) * dH[k - 1L]
    if (type == "ph") {
      L[k] <- mean(delta * abs(y - fit[, k]) + (delta - 2 * comp) * fit[, k])
      next
    }
    M <- delta * had_event_by(y, fit[, k]) - comp
    L[k] <- if (type == "martingale") mean(abs(M)) else {
      r <- -2 * (M + delta * log(pmax(delta - M, 1e-300)))
      mean(sqrt(pmax(r, 0)))
    }
  }
  mean(L)
}

## Folds stratified by the event indicator ---------------------------------------
make_folds <- function(delta, K) {
  fold <- integer(length(delta))
  for (lev in 0:1) {
    idx <- which(delta == lev)
    fold[idx] <- sample(rep_len(seq_len(K), length(idx)))
  }
  fold
}

## Tuning -----------------------------------------------------------------------
## rule = "cv":       K-fold CV over multipliers of the pivotal penalty
## rule = "pivotal":  fixed multiplier, no CV (fast; used for MDS)
## Returns the grid-wise penalty vector for the full (training) sample.
tune_lambda <- function(y, delta, Z, taus, rule = c("cv", "pivotal"),
                        mult = c(0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1, 1.25),
                        fixed_mult = 0.5, K = 5L, loss = "ph",
                        nsim = 500L, cores = 1L, foldid = NULL, lambda = NULL) {
  rule <- match.arg(rule)
  tau <- taus[length(taus)]
  lam0 <- pivotal_lambda(Z[, -1L, drop = FALSE], tau, nsim)
  if (!is.null(lambda))                         # user-supplied penalty at the target tau
    return(list(lambda = lambda_grid(lambda, taus), mult = lambda / lam0,
                pivotal = lam0, cv = NULL))
  ## the multiplier grid is returned also under a fixed multiplier, so that the
  ## "path" and "size" screens can move down the grid from fixed_mult
  if (rule == "pivotal")
    return(list(lambda = lambda_grid(fixed_mult * lam0, taus),
                mult = fixed_mult, grid = mult, pivotal = lam0, cv = NULL))
  fold <- if (is.null(foldid)) make_folds(delta, K) else as.integer(foldid)
  fold <- match(fold, sort(unique(fold)))       # relabel 1..K
  K <- max(fold)
  ## pivotal penalty on each training fold (same random draws per fold)
  lam_fold <- vapply(seq_len(K), function(k)
    pivotal_lambda(Z[fold != k, -1L, drop = FALSE], tau, nsim), numeric(1))
  jobs <- expand.grid(k = seq_len(K), l = seq_along(mult))
  run <- function(j) {
    k <- jobs$k[j]; l <- jobs$l[j]
    tr <- fold != k
    f <- hdcqr_path(y[tr], delta[tr], Z[tr, , drop = FALSE], taus,
                    lambda_grid(mult[l] * lam_fold[k], taus))
    if (f$status != "ok") return(NA_real_)
    cv_loss(y[!tr], delta[!tr], Z[!tr, , drop = FALSE], f$B, taus, loss)
  }
  res <- if (cores > 1L && .Platform$OS.type == "unix")
           unlist(parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = cores,
                                     mc.set.seed = FALSE))
         else vapply(seq_len(nrow(jobs)), run, numeric(1))
  score <- matrix(res, K, length(mult))
  mean_score <- colMeans(score)                 # NA if any fold failed
  if (all(is.na(mean_score))) stop("every CV candidate failed")
  best <- which(mean_score == min(mean_score, na.rm = TRUE))[1L]  # smallest multiplier on ties
  list(lambda = lambda_grid(mult[best] * lam0, taus), mult = mult[best], grid = mult,
       pivotal = lam0, cv = mean_score)
}
