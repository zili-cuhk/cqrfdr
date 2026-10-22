## ---------------------------------------------------------------------------
## Data-generating processes for the censored-quantile-regression simulations
## ---------------------------------------------------------------------------
##
## Two linear quantile models are provided.
##
##  model = "location"  (the design used in the current draft)
##      log T = alpha + Z' b* + a_F { E - F^{-1}(tau*) },   b* = beta(tau*)
##    This is a location-shift model: within one data set every quantile has
##    the same slope vector b*; only the intercept moves with u.  Changing the
##    target tau* changes the *data-generating process*, not the quantile
##    being estimated within one data set.
##
##  model = "hetero"  (recommended additional design)
##      log T = alpha + Z' b + sigma (1 + gamma * Z_h) E,   Z_h = Phi(W_h) in (0,1)
##    Q(u | Z) = alpha + sigma q(u) + Z' b + gamma sigma q(u) Z_h,
##    so the coefficient of Z_h is gamma*sigma*q(u): it is zero at the median
##    of a symmetric F and non-zero at other quantiles.  The active set
##    therefore genuinely depends on the target quantile.
##
## Covariates: Toeplitz/AR(1) correlation rho^|j-k|, Gaussian or
## multivariate t (scaled to unit variance).  Censoring: log C ~ N(c, sd_c)
## independent of (T, Z), with c calibrated to a target censoring rate.

## Error distributions --------------------------------------------------------
## edist: "normal", "t" (degrees of freedom edf; edf = 1 is Cauchy), "exp".
## Scaling (err_scale):
##   "sd"  : unit variance (the draft's convention; needs edf > 2)
##   "iqr" : same interquartile range as N(0,1); defined for every edf, so a
##           sweep over edf (including Cauchy) changes only the tails.
## Legacy name "t3" is accepted and means edist = "t", edf = 3.
.err_scale <- function(dist, df, scale) {
  if (dist != "t") return(1)
  if (scale == "sd") {
    if (df <= 2) stop("t errors with df <= 2 have no variance; use err_scale = 'iqr'")
    1 / sqrt(df / (df - 2))
  } else stats::qnorm(0.75) / stats::qt(0.75, df)
}
.norm_dist <- function(dist, df) if (dist == "t3") list("t", 3) else list(dist, df)
err_q <- function(u, dist, df = 3, scale = "sd") {
  nd <- .norm_dist(dist, df); dist <- nd[[1]]; df <- nd[[2]]
  switch(dist,
         normal = stats::qnorm(u),
         t      = stats::qt(u, df) * .err_scale(dist, df, scale),
         exp    = stats::qexp(u) - 1,
         stop("unknown error distribution: ", dist))
}
err_r <- function(n, dist, df = 3, scale = "sd") {
  nd <- .norm_dist(dist, df); dist <- nd[[1]]; df <- nd[[2]]
  switch(dist,
         normal = stats::rnorm(n),
         t      = stats::rt(n, df) * .err_scale(dist, df, scale),
         exp    = stats::rexp(n) - 1,
         stop("unknown error distribution: ", dist))
}

## Covariates -----------------------------------------------------------------
## AR(1) recursion gives Cov(Z_j, Z_k) = rho^|j-k| exactly; O(np).
## xdist = "t": multivariate t_xdf = Gaussian row / sqrt(chi2_xdf / xdf), i.e.
## one common random scale per subject.  For xdf > 2 it is rescaled to unit
## variance (Cov = Toeplitz); for xdf <= 2 (xdf = 1: multivariate Cauchy) no
## covariance exists and the Toeplitz matrix is only the scale matrix.
x_scale <- function(df) if (is.finite(df) && df > 2) sqrt((df - 2) / df) else 1
sim_covariates <- function(n, p, rho, xdist = "gaussian", df = 3) {
  nd <- .norm_dist(xdist, df); xdist <- nd[[1]]; df <- nd[[2]]
  X <- matrix(stats::rnorm(n * p), n, p)
  if (p > 1L && rho != 0) {
    s <- sqrt(1 - rho^2)
    for (j in 2:p) X[, j] <- rho * X[, j - 1L] + s * X[, j]
  }
  if (xdist == "t") X <- X * (x_scale(df) / sqrt(stats::rchisq(n, df) / df))
  X
}
## CDF of one covariate (used for the probability-integral transform)
x_cdf <- function(x, xdist, df) {
  nd <- .norm_dist(xdist, df)
  if (nd[[1]] == "gaussian") stats::pnorm(x) else stats::pt(x / x_scale(nd[[2]]), nd[[2]])
}

## Scenario defaults --------------------------------------------------------
default_spec <- function(...) {
  spec <- list(
    n = 800L, p = 800L, s = 10L, rho = 0.5, tau = 0.5,
    xdist = "gaussian", edist = "normal", cens = 0.30,
    xdf = 5, edf = 5, err_scale = "sd",
    model = "location",
    alpha = 0.5, beta_level = 0.8, beta_slope = 0.8, err_sd = 0.5,
    ## hetero model only
    gamma = 3, sigma = 0.5, hetero_index = NA_integer_,
    censor_sd = 1, calib_n = 200000L
  )
  args <- list(...)
  spec[names(args)] <- args
  if (spec$model == "hetero" && is.na(spec$hetero_index))
    spec$hetero_index <- spec$s + 1L
  spec
}

eq <- function(u, spec) err_q(u, spec$edist, spec$edf, spec$err_scale)
er <- function(n, spec) err_r(n, spec$edist, spec$edf, spec$err_scale)

## True coefficient vector (intercept first) of Q_{log T}(u | Z) in this DGP
true_beta <- function(spec, u = spec$tau) {
  p <- spec$p; s <- spec$s
  b <- numeric(p)
  if (spec$model == "location") {
    sgn   <- rep(c(1, -1), length.out = s)
    slope <- rep(c(1, 1, -1, -1), length.out = s)
    b[seq_len(s)] <- spec$beta_level * sgn +
      spec$beta_slope * (spec$tau - 0.5) * slope           # b* = beta(tau*)
    aF <- spec$err_sd                                      # unit-sd errors
    intercept <- spec$alpha + aF * (eq(u, spec) - eq(spec$tau, spec))
  } else if (spec$model == "hetero") {
    b[seq_len(s)] <- spec$beta_level * rep(c(1, -1), length.out = s)
    qu <- eq(u, spec)
    b[spec$hetero_index] <- spec$gamma * spec$sigma * qu
    intercept <- spec$alpha + spec$sigma * qu
  } else stop("unknown model")
  c(intercept, b)
}

## Active set at the target quantile (indices 1..p, no intercept)
true_active <- function(spec, tol = 1e-12) which(abs(true_beta(spec)[-1L]) > tol)

## Latent log event times --------------------------------------------------
sim_latent <- function(n, spec) {
  X <- sim_covariates(n, spec$p, spec$rho, spec$xdist, spec$xdf)
  if (spec$model == "location") {
    b <- true_beta(spec)[-1L]
    e <- er(n, spec)
    logT <- spec$alpha + drop(X %*% b) + spec$err_sd * (e - eq(spec$tau, spec))
  } else {
    h <- spec$hetero_index
    ## probability-integral transform keeps the correlation structure but
    ## gives a positive scale covariate in (0,1)
    X[, h] <- x_cdf(X[, h], spec$xdist, spec$xdf)
    b <- numeric(spec$p)
    b[seq_len(spec$s)] <- spec$beta_level * rep(c(1, -1), length.out = spec$s)
    e <- er(n, spec)
    logT <- spec$alpha + drop(X %*% b) + spec$sigma * (1 + spec$gamma * X[, h]) * e
  }
  list(X = X, logT = logT)
}

## Censoring calibration ------------------------------------------------------
## Returns the location c of log C ~ N(c, sd_c) giving the target censoring
## rate, plus a descriptive diagnostic at the target quantile:
##   at_risk = P(C >= Q_T(tau | Z) | Z) summarised over Z
## (the theory, Condition 8/9, needs it bounded away from 0 on the support of
## Z; with Gaussian Z and a large signal variance it is not).
calibrate_censoring <- function(spec) {
  n <- spec$calib_n
  ## only the columns entering the linear predictor are needed; simulate
  ## the full AR chain up to the last relevant index
  last <- if (spec$model == "hetero") max(spec$s, spec$hetero_index) else spec$s
  sp <- spec; sp$p <- last
  lat <- sim_latent(n, sp)
  g <- stats::rnorm(n)
  f <- function(c0) mean(c0 + spec$censor_sd * g < lat$logT) - spec$cens
  rng <- range(lat$logT) + c(-10, 10) * spec$censor_sd
  c0 <- stats::uniroot(f, rng, tol = 1e-10)$root
  ## identifiability diagnostic
  bt <- true_beta(sp)
  q_tau <- bt[1L] + drop(lat$X %*% bt[-1L])
  at_risk <- stats::pnorm((c0 - q_tau) / spec$censor_sd)   # P(log C >= Q_tau(Z) | Z)
  list(c0 = c0, achieved = mean(c0 + spec$censor_sd * g < lat$logT),
       at_risk_mean = mean(at_risk),
       at_risk_q05 = unname(stats::quantile(at_risk, 0.05)),
       at_risk_q01 = unname(stats::quantile(at_risk, 0.01)))
}

## One observed data set ------------------------------------------------------
sim_data <- function(spec, c0) {
  lat <- sim_latent(spec$n, spec)
  logC <- c0 + spec$censor_sd * stats::rnorm(spec$n)
  delta <- as.integer(lat$logT <= logC)
  list(y = pmin(lat$logT, logC), delta = delta, X = lat$X,
       active = true_active(spec), beta = true_beta(spec),
       cens_rate = 1 - mean(delta))
}

## Exported simulator ------------------------------------------------------------
sim_cqr <- function(n = 400L, p = 200L, s = 10L, beta = 0.8, beta_slope = beta,
                    tau = 0.5, rho = 0.5, cens = 0.3,
                    xdist = c("gaussian", "t"), edist = c("normal", "t", "exp"),
                    model = c("location", "hetero"), xdf = 5, edf = 5,
                    err_sd = 0.5, alpha = 0.5, gamma = 3, sigma = 0.5,
                    censor_sd = 1, c0 = NULL, calib_n = 200000L) {
  spec <- default_spec(n = as.integer(n), p = as.integer(p), s = as.integer(s),
                       beta_level = beta, beta_slope = beta_slope, tau = tau, rho = rho,
                       cens = cens, xdist = match.arg(xdist), edist = match.arg(edist),
                       model = match.arg(model), xdf = xdf, edf = edf, err_sd = err_sd,
                       alpha = alpha, gamma = gamma, sigma = sigma, censor_sd = censor_sd,
                       calib_n = as.integer(calib_n))
  if (spec$s < 1L || spec$s > spec$p) stop("need 1 <= s <= p")
  if (spec$model == "hetero" && spec$hetero_index > spec$p) stop("model = \"hetero\" needs p > s")
  cal <- NULL
  if (is.null(c0)) {
    if (cens <= 0) c0 <- Inf else {
      cal <- calibrate_censoring(spec); c0 <- cal$c0
    }
  }
  d <- sim_data(spec, c0)
  colnames(d$X) <- paste0("X", seq_len(spec$p))
  b <- d$beta; names(b) <- c("(Intercept)", colnames(d$X))
  list(x = d$X, y = d$y, delta = d$delta, time = exp(d$y), active = d$active,
       beta = b, cens_rate = d$cens_rate, c0 = c0, calibration = cal, spec = spec)
}

true_coef <- function(spec, tau = spec$tau) {
  b <- true_beta(spec, tau)
  names(b) <- c("(Intercept)", paste0("X", seq_len(spec$p)))
  b
}
