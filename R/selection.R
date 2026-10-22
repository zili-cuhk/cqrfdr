## ---------------------------------------------------------------------------
## Variable selection with FDR control: SS-MS, MDS and L-HDCQR
## ---------------------------------------------------------------------------

## Mirror statistic and data-driven threshold --------------------------------
mirror_stat <- function(b1, b2) sign(b1 * b2) * (abs(b1) + abs(b2))

## smallest t > 0 with (offset + #{M <= -t}) / max(#{M >= t}, 1) <= q
## offset = 1: U_+ (finite-sample bound); offset = 0: U (Dai et al., 2023)
mirror_threshold <- function(M, q, offset = 1) {
  if (length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1) stop("'q' must be one number in (0, 1)")
  ts <- sort(unique(abs(M[M != 0])))
  if (!length(ts)) return(Inf)
  neg <- vapply(ts, function(t) sum(M <= -t), numeric(1))
  pos <- vapply(ts, function(t) sum(M >= t), numeric(1))
  ok <- which((offset + neg) / pmax(pos, 1) <= q)
  if (length(ok)) ts[ok[1L]] else Inf
}

## Standardisation on a subset of rows -----------------------------------------
## std = "sd": mean / standard deviation; "robust": median / MAD.
## Returns the design with an intercept column; centre and scale as attributes.
standardize_by <- function(X, rows, std = c("sd", "robust")) {
  std <- match.arg(std)
  Xr <- X[rows, , drop = FALSE]
  if (std == "sd") {
    ctr <- colMeans(Xr)
    scl <- sqrt(colSums(sweep(Xr, 2L, ctr)^2) / (length(rows) - 1))
  } else {
    ctr <- apply(Xr, 2L, stats::median)
    scl <- apply(Xr, 2L, stats::mad)
  }
  bad <- which(!is.finite(scl) | scl <= 1e-10)
  if (length(bad))
    stop(sprintf("covariate(s) %s are constant on the fitting subset",
                 paste(utils::head(bad, 5L), collapse = ", ")))
  Z <- cbind(`(Intercept)` = 1, sweep(sweep(X, 2L, ctr), 2L, scl, "/"))
  attr(Z, "center") <- ctr; attr(Z, "scale") <- scl
  Z
}

fail <- function(status, ...)
  c(list(status = status, selected = integer(0), failed = TRUE), list(...))

## Screening cap ----------------------------------------------------------------
## "d2": floor(d2 / log d2), d2 = events in split 2 (the paper's simulations);
## "d1": the same with the events in split 1 (makes the screen D1-measurable);
## an integer: fixed cap; "none": no cap.
screen_cap_value <- function(screen_cap, d1, d2) {
  f <- function(d) if (d > 2) max(1L, as.integer(floor(d / log(d)))) else 1L
  if (identical(screen_cap, "d2") || identical(screen_cap, "theory")) return(f(d2))
  if (identical(screen_cap, "d1")) return(f(d1))
  if (identical(screen_cap, "none") || is.null(screen_cap)) return(NULL)
  if (is.numeric(screen_cap) && length(screen_cap) == 1L && screen_cap >= 1)
    return(as.integer(screen_cap))
  stop("'screen_cap' must be \"d2\", \"d1\", \"none\" or a positive integer")
}

## Algorithm 1: one split, mirror statistics -----------------------------------
## Internal engine shared by ssms() and mds().  X is the raw n x p matrix.
.ssms_engine <- function(y, delta, X, taus, q = 0.2, offset = 1,
                         lambda_rule = "cv", mult = NULL, fixed_mult = 0.5, K = 5L,
                         loss = "ph", screen_cap = "d2", coef_tol = 1e-6,
                         cores = 1L, split = NULL, std = "sd", foldid = NULL,
                         screen_rule = c("path", "cv", "size"), lambda = NULL,
                         nsim = 500L) {
  screen_rule <- match.arg(screen_rule)
  n <- length(y); p <- ncol(X); m <- length(taus)
  i1 <- if (is.null(split)) sort(sample.int(n, floor(n / 2))) else sort(as.integer(split))
  i2 <- setdiff(seq_len(n), i1)
  Z <- standardize_by(X, i1, std)                # split-1 centring/scaling
  targs <- list(y = y[i1], delta = delta[i1], Z = Z[i1, , drop = FALSE], taus = taus,
                rule = lambda_rule, fixed_mult = fixed_mult, K = K, loss = loss,
                cores = cores, nsim = nsim, lambda = lambda)
  if (!is.null(foldid)) targs$foldid <- foldid[i1]
  if (!is.null(mult)) targs$mult <- mult
  tun <- tryCatch(do.call(tune_lambda, targs), error = function(e) e)
  if (inherits(tun, "error"))
    return(fail("tuning_failed", message = conditionMessage(tun), split = i1))
  f1 <- hdcqr_path(y[i1], delta[i1], Z[i1, , drop = FALSE], taus, tun$lambda)
  if (f1$status != "ok")
    return(fail(paste0("stage1_", f1$status), tau_fail = taus[f1$k], split = i1))
  b1 <- f1$B[-1L, m]
  screen <- which(abs(b1) > coef_tol)
  cap <- screen_cap_value(screen_cap, sum(delta[i1]), sum(delta[i2]))
  mult_screen <- tun$mult
  ## "size": move down the multiplier grid until at least `cap` candidates are
  ## nonzero, then keep the `cap` largest |b1|.
  if (screen_rule == "size" && !is.null(cap) && length(screen) < cap && !is.null(tun$grid)) {
    for (mm in rev(tun$grid[tun$grid < tun$mult])) {
      f <- hdcqr_path(y[i1], delta[i1], Z[i1, , drop = FALSE], taus,
                      lambda_grid(mm * tun$pivotal, taus))
      if (f$status != "ok") break
      f1 <- f; b1 <- f1$B[-1L, m]; mult_screen <- mm
      screen <- which(abs(b1) > coef_tol)
      if (length(screen) >= cap) break
    }
  }
  ## "path" (default): candidates are ranked by the order in which they enter
  ## the lasso path as the multiplier decreases from the CV choice (ties:
  ## larger |coefficient| first); each keeps its coefficient from the fit at
  ## which it entered.  The CV support is always kept.  b1 stays a function of
  ## split 1 only, so the mirror-statistic argument is unchanged.
  if (screen_rule == "path" && !is.null(cap) && length(screen) < cap && !is.null(tun$grid)) {
    ord <- screen[order(-abs(b1[screen]))]
    b1e <- numeric(p); b1e[screen] <- b1[screen]
    for (mm in rev(tun$grid[tun$grid < tun$mult])) {
      f <- hdcqr_path(y[i1], delta[i1], Z[i1, , drop = FALSE], taus,
                      lambda_grid(mm * tun$pivotal, taus))
      if (f$status != "ok") break
      bb <- f$B[-1L, m]; mult_screen <- mm
      new <- setdiff(which(abs(bb) > coef_tol), ord)
      if (length(new)) {
        new <- new[order(-abs(bb[new]))]
        b1e[new] <- bb[new]; ord <- c(ord, new)
      }
      if (length(ord) >= cap) break
    }
    screen <- sort(ord[seq_len(min(cap, length(ord)))]); b1 <- b1e
  }
  full_screen <- length(screen)
  if (!is.null(cap) && length(screen) > cap)
    screen <- sort(screen[order(-abs(b1[screen]))[seq_len(cap)]])
  out <- list(status = "ok", selected = integer(0), failed = FALSE,
              screen = screen, screen_full_size = full_screen, cap = cap,
              mult = tun$mult, mult_screen = mult_screen, pivotal = tun$pivotal,
              cv = tun$cv, mult_grid = tun$grid, split = i1,
              b1 = b1, b2 = numeric(p), M = numeric(p), threshold = Inf,
              events = c(split1 = sum(delta[i1]), split2 = sum(delta[i2])))
  out$by_q <- lapply(q, function(qq) list(q = qq, selected = integer(0), threshold = Inf,
                                          status = "empty_screen"))
  if (!length(screen)) { out$status <- "empty_screen"; return(out) }
  f2 <- hdcqr_path(y[i2], delta[i2], Z[i2, c(1L, screen + 1L), drop = FALSE], taus, 0)
  if (f2$status != "ok")
    return(fail(paste0("stage2_", f2$status), tau_fail = taus[f2$k], screen = screen,
                split = i1))
  out$b2[screen] <- f2$B[-1L, m]
  out$M[screen] <- mirror_stat(b1[screen], out$b2[screen])
  out$by_q <- lapply(q, function(qq) {
    thr <- mirror_threshold(out$M, qq, offset)
    list(q = qq, selected = which(out$M >= thr), threshold = thr,
         status = if (is.finite(thr)) "ok" else "no_threshold")
  })
  out$selected <- out$by_q[[1L]]$selected
  out$threshold <- out$by_q[[1L]]$threshold
  out$status <- out$by_q[[1L]]$status
  out
}

## Input handling ------------------------------------------------------------------
.prep <- function(x, y, delta, tau, q, taus, nu, n_grid, log_time) {
  if (inherits(y, "Surv")) {
    if (!is.null(delta)) warning("'delta' ignored: the event indicator is taken from the Surv object")
    if (attr(y, "type") != "right") stop("only right-censored Surv objects are supported")
    delta <- y[, "status"]; y <- y[, "time"]
    if (is.null(log_time)) log_time <- TRUE
  }
  if (isTRUE(log_time)) {
    if (any(y <= 0)) stop("log_time = TRUE needs positive times")
    y <- log(y)
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("'x' must be a numeric matrix")
  y <- as.numeric(y); delta <- as.integer(delta)
  n <- nrow(x)
  if (length(y) != n || length(delta) != n) stop("'x', 'y' and 'delta' have different numbers of observations")
  if (anyNA(x) || anyNA(y) || anyNA(delta)) stop("missing values are not allowed")
  if (!all(delta %in% 0:1)) stop("'delta' must be 0 (censored) or 1 (event)")
  if (sum(delta) < 10L) stop("fewer than 10 events")
  if (is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
  if (!is.null(q) && (any(!is.finite(q)) || any(q <= 0 | q >= 1))) stop("'q' must lie in (0, 1)")
  if (is.null(taus)) {
    if (length(tau) != 1L || tau <= 0 || tau >= 1) stop("'tau' must be one number in (0, 1)")
    if (nu >= tau) stop("'nu' must be smaller than 'tau'")
    taus <- seq(nu, tau, length.out = n_grid)
  } else {
    taus <- sort(as.numeric(taus))
    if (any(taus <= 0 | taus >= 1)) stop("'taus' must lie in (0, 1)")
    tau <- taus[length(taus)]
  }
  list(x = x, y = y, delta = delta, taus = taus, tau = tau)
}

.lambda_args <- function(lambda) {
  if (is.numeric(lambda)) {
    if (length(lambda) != 1L || lambda <= 0) stop("a numeric 'lambda' must be one positive number")
    return(list(rule = "pivotal", value = lambda))
  }
  list(rule = match.arg(lambda, c("cv", "pivotal")), value = NULL)
}

## Exported: SS-MS -------------------------------------------------------------------
ssms <- function(x, y, delta = NULL, tau = 0.5, q = c(0.1, 0.2),
                 nu = 0.1, n_grid = 40L, taus = NULL, log_time = NULL,
                 lambda = c("cv", "pivotal"), mult_grid = NULL, fixed_mult = 0.5,
                 nfolds = 5L, foldid = NULL, loss = c("ph", "martingale", "deviance"),
                 screen = c("path", "cv", "size"), screen_cap = "d2",
                 offset = 1, split = NULL, std = c("sd", "robust"),
                 coef_tol = 1e-6, nsim = 500L, cores = 1L) {
  cl <- match.call()
  d <- .prep(x, y, delta, tau, q, taus, nu, n_grid, log_time)
  lam <- .lambda_args(lambda)
  if (!is.null(foldid) && length(foldid) != nrow(d$x)) stop("'foldid' must have one entry per observation")
  r <- .ssms_engine(d$y, d$delta, d$x, d$taus, q = q, offset = offset,
                    lambda_rule = lam$rule, mult = mult_grid, fixed_mult = fixed_mult,
                    K = nfolds, loss = match.arg(loss), screen_cap = screen_cap,
                    coef_tol = coef_tol, cores = cores, split = split,
                    std = match.arg(std), foldid = foldid, screen_rule = match.arg(screen),
                    lambda = lam$value, nsim = nsim)
  if (isTRUE(r$failed))
    warning("SS-MS failed (", r$status, "); no variables selected", call. = FALSE)
  vn <- colnames(d$x)
  res <- list(call = cl, method = "SS-MS", status = r$status, failed = isTRUE(r$failed),
              q = q, tau = d$tau, taus = d$taus, n = nrow(d$x), p = ncol(d$x),
              events = sum(d$delta), varnames = vn, offset = offset,
              selected = .named_sel(r$by_q, q, vn),
              threshold = if (!is.null(r$by_q)) vapply(r$by_q, `[[`, 0, "threshold") else rep(NA_real_, length(q)),
              M = stats::setNames(if (is.null(r$M)) rep(NA_real_, ncol(d$x)) else r$M, vn),
              b1 = if (!is.null(r$b1)) stats::setNames(r$b1, vn),
              b2 = if (!is.null(r$b2)) stats::setNames(r$b2, vn),
              screen = if (!is.null(r$screen)) stats::setNames(r$screen, vn[r$screen]),
              screen_rule = match.arg(screen), cap = r$cap, split = r$split,
              tuning = list(rule = lam$rule, mult = r$mult, mult_screen = r$mult_screen,
                            pivotal = r$pivotal, mult_grid = r$mult_grid, cv = r$cv),
              events_split = r$events, message = r$message)
  class(res) <- c("cqrfdr_ssms", "cqrfdr")
  res
}

.named_sel <- function(by_q, q, vn) {
  s <- if (is.null(by_q)) replicate(length(q), integer(0), simplify = FALSE)
       else lapply(by_q, function(b) b$selected)
  s <- lapply(s, function(v) stats::setNames(v, vn[v]))
  names(s) <- format(q)
  s
}

## Exported: DS, the single data-splitting procedure of Dai et al. (2023) -------
## Half-half split, penalty by K-fold CV (10 in their simulations) on the
## first half, candidate set = CV support without a cap, unpenalized
## Peng-Huang refit on the second half (their OLS step), cutoff U.
ds <- function(x, y, delta = NULL, tau = 0.5, q = c(0.1, 0.2), nu = 0.1, n_grid = 40L,
               taus = NULL, log_time = NULL, nfolds = 10L, foldid = NULL,
               loss = c("ph", "martingale", "deviance"), split = NULL,
               std = c("sd", "robust"), coef_tol = 1e-6, nsim = 500L, cores = 1L) {
  r <- ssms(x, y, delta, tau = tau, q = q, nu = nu, n_grid = n_grid, taus = taus,
            log_time = log_time, lambda = "cv", nfolds = nfolds, foldid = foldid,
            loss = match.arg(loss), screen = "cv", screen_cap = "none", offset = 0,
            split = split, std = match.arg(std), coef_tol = coef_tol, nsim = nsim,
            cores = cores)
  r$call <- match.call(); r$method <- "DS"
  class(r) <- c("cqrfdr_ds", class(r))
  r
}

## Multiple splits with inclusion-rate aggregation (Dai et al., 2023) ----------
## Shared by mds() (each split: CV-support screen, cutoff U, as in Dai et al.)
## and msms() (each split: SS-MS, i.e. path-ordered screen with the cap and
## cutoff U+).  The penalty multiplier is chosen once by CV on one random half
## (tune = "cv_once") and reused, with each split's own pivotal level.
.multi_split <- function(cl, method, x, y, delta, tau, q, B, nu, n_grid, taus, log_time,
                         tune, ties, fixed_mult, mult_grid, nfolds, foldid, loss,
                         screen_rule, screen_cap, offset, std, coef_tol, nsim, cores, verbose) {
  d <- .prep(x, y, delta, tau, q, taus, nu, n_grid, log_time)
  X <- d$x; y <- d$y; delta <- d$delta; taus <- d$taus
  p <- ncol(X); n <- nrow(X)
  if (!is.null(foldid) && length(foldid) != n) stop("'foldid' must have one entry per observation")
  if (length(B) != 1L || B < 1) stop("'B' must be a positive integer")
  mult_used <- fixed_mult
  if (tune == "cv_once") {                       # one CV on a random half, reused
    i1 <- sort(sample.int(n, floor(n / 2)))
    Z <- standardize_by(X, i1, std)
    targs <- list(y = y[i1], delta = delta[i1], Z = Z[i1, , drop = FALSE], taus = taus,
                  rule = "cv", K = nfolds, loss = loss, cores = cores, nsim = nsim)
    if (!is.null(foldid)) targs$foldid <- foldid[i1]
    if (!is.null(mult_grid)) targs$mult <- mult_grid
    tun <- tryCatch(do.call(tune_lambda, targs), error = function(e) e)
    if (inherits(tun, "error")) stop("tuning failed: ", conditionMessage(tun))
    mult_used <- tun$mult
  }
  incl <- matrix(0, p, length(q), dimnames = list(colnames(X), format(q)))
  nsel <- matrix(0L, B, length(q), dimnames = list(NULL, format(q)))
  nf <- 0L; status <- character(B); failed <- logical(B); mults <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    ## "cv_each": every split is a complete single-split fit with its own
    ## K-fold CV on its first half (Dai et al.: MDS repeats DS)
    eargs <- list(y, delta, X, taus, q = q, offset = offset,
                  lambda_rule = if (tune == "cv_each") "cv" else "pivotal",
                  fixed_mult = mult_used, K = nfolds, cores = cores, loss = loss,
                  screen_cap = screen_cap, coef_tol = coef_tol, std = std,
                  screen_rule = screen_rule, nsim = nsim)
    if (!is.null(mult_grid)) eargs$mult <- mult_grid
    r <- do.call(.ssms_engine, eargs)
    status[b] <- r$status
    if (tune == "cv_each" && !is.null(r$mult)) mults[b] <- r$mult
    if (verbose) message(sprintf("split %d/%d: %s", b, B, r$status))
    if (isTRUE(r$failed)) { nf <- nf + 1L; failed[b] <- TRUE; next }
    for (j in seq_along(q)) {
      sel <- r$by_q[[j]]$selected
      nsel[b, j] <- length(sel)
      if (length(sel)) incl[sel, j] <- incl[sel, j] + 1 / length(sel)
    }
  }
  if (nf == B) stop("all ", B, " splits failed (", paste(unique(status), collapse = ", "), ")")
  incl <- incl / (B - nf)
  sel <- lapply(seq_along(q), function(j) {
    I <- incl[, j]
    o <- if (ties == "dai") order(I) else order(I, stats::runif(p))
    cs <- cumsum(I[o]); l <- max(c(0L, which(cs <= q[j])))
    s <- if (ties == "dai") {
      cut <- if (l == 0L) -Inf else I[o][l]
      which(I > cut & I > 0)
    } else {
      keep <- if (l < p) o[(l + 1L):p] else integer(0)
      sort(keep[I[keep] > 0])
    }
    stats::setNames(s, colnames(X)[s])
  })
  names(sel) <- format(q)
  ok <- !failed
  res <- list(call = cl, method = method, status = "ok", failed = FALSE, q = q,
              tau = d$tau, taus = taus, n = n, p = p, events = sum(delta),
              varnames = colnames(X), selected = sel, inclusion = incl, B = B,
              n_failed_splits = nf, split_status = status,
              split_sizes = nsel[ok, , drop = FALSE],
              share_empty = colMeans(nsel[ok, , drop = FALSE] == 0),
              mult = if (tune == "cv_each") stats::median(mults, na.rm = TRUE) else mult_used,
              split_mult = if (tune == "cv_each") mults, tune = tune, ties = ties,
              screen_rule = screen_rule, offset = offset)
  class(res) <- c(if (method == "MS-MS") "cqrfdr_msms", "cqrfdr_mds", "cqrfdr")
  res
}

## Exported: multiple data splitting (Dai et al.) ---------------------------------
mds <- function(x, y, delta = NULL, tau = 0.5, q = c(0.1, 0.2), B = 50L,
                nu = 0.1, n_grid = 40L, taus = NULL, log_time = NULL,
                tune = c("cv_each", "cv_once", "pivotal"), ties = c("dai", "random"),
                fixed_mult = 0.5, mult_grid = NULL, nfolds = 10L, foldid = NULL,
                loss = c("ph", "martingale", "deviance"), screen_cap = "none",
                offset = 0, std = c("sd", "robust"), coef_tol = 1e-6, nsim = 500L,
                cores = 1L, verbose = FALSE) {
  .multi_split(match.call(), "MDS", x, y, delta, tau, q, B, nu, n_grid, taus, log_time,
               match.arg(tune), match.arg(ties), fixed_mult, mult_grid, nfolds, foldid,
               match.arg(loss), "cv", screen_cap, offset, match.arg(std), coef_tol, nsim,
               cores, verbose)
}

## Exported: multiple splits of SS-MS ------------------------------------------------
msms <- function(x, y, delta = NULL, tau = 0.5, q = c(0.1, 0.2), B = 50L,
                 nu = 0.1, n_grid = 40L, taus = NULL, log_time = NULL,
                 tune = c("cv_once", "cv_each", "pivotal"), ties = c("random", "dai"),
                 fixed_mult = 0.5, mult_grid = NULL, nfolds = 5L, foldid = NULL,
                 loss = c("ph", "martingale", "deviance"),
                 screen = c("path", "cv", "size"), screen_cap = "d2",
                 offset = 1, std = c("sd", "robust"), coef_tol = 1e-6, nsim = 500L,
                 cores = 1L, verbose = FALSE) {
  .multi_split(match.call(), "MS-MS", x, y, delta, tau, q, B, nu, n_grid, taus, log_time,
               match.arg(tune), match.arg(ties), fixed_mult, mult_grid, nfolds, foldid,
               match.arg(loss), match.arg(screen), screen_cap, offset, match.arg(std),
               coef_tol, nsim, cores, verbose)
}

## Exported: L-HDCQR fit -----------------------------------------------------------
hdcqr <- function(x, y, delta = NULL, tau = 0.5, nu = 0.1, n_grid = 40L, taus = NULL,
                  log_time = NULL, lambda = c("cv", "pivotal", "zheng"), mult_grid = NULL,
                  fixed_mult = 0.5, nfolds = 5L, foldid = NULL,
                  loss = c("ph", "martingale", "deviance"), std = c("sd", "robust"),
                  coef_tol = 1e-6, nsim = 500L, cores = 1L,
                  increments = c(0, exp(seq(log(0.02), log(20), length.out = 12)))) {
  cl <- match.call()
  d <- .prep(x, y, delta, tau, NULL, taus, nu, n_grid, log_time)
  if (is.character(lambda)) lambda <- match.arg(lambda, c("cv", "pivotal", "zheng"))
  zheng <- identical(lambda, "zheng")
  Z <- standardize_by(d$x, seq_len(nrow(d$x)), match.arg(std))
  if (zheng) {
    f <- .hdcqr_zheng(d$y, d$delta, Z, d$taus, K = nfolds, foldid = foldid,
                      increments = increments, nsim = nsim, cores = cores)
    if (f$status != "ok")
      stop(sprintf("the sequential fit failed at tau = %.3f (%s)", d$taus[f$k], f$status))
    tun <- list(lambda = f$lambda, mult = NA_real_, pivotal = f$lambda[1L], cv = NULL,
                grid = increments, increment = f$increment)
  } else {
    lam <- .lambda_args(lambda)
    targs <- list(y = d$y, delta = d$delta, Z = Z, taus = d$taus, rule = lam$rule,
                  fixed_mult = fixed_mult, K = nfolds, loss = match.arg(loss),
                  cores = cores, foldid = foldid, nsim = nsim, lambda = lam$value)
    if (!is.null(mult_grid)) targs$mult <- mult_grid
    tun <- do.call(tune_lambda, targs)
    f <- hdcqr_path(d$y, d$delta, Z, d$taus, tun$lambda)
    if (f$status != "ok")
      stop(sprintf("the path stopped at tau = %.3f (%s): %s", d$taus[f$k], f$status, f$message))
  }
  m <- length(d$taus)
  vn <- colnames(d$x)
  B <- f$B; rownames(B) <- c("(Intercept)", vn); colnames(B) <- format(d$taus, digits = 4)
  sel <- which(abs(B[-1L, m]) > coef_tol)
  res <- list(call = cl, coef_std = B, center = attr(Z, "center"), scale = attr(Z, "scale"),
              selected = stats::setNames(sel, vn[sel]), tau = d$tau, taus = d$taus,
              lambda = tun$lambda, mult = tun$mult, pivotal = tun$pivotal, cv = tun$cv,
              mult_grid = tun$grid, increment = tun$increment,
              tuning = if (is.numeric(lambda)) "fixed" else lambda,
              n = nrow(d$x), p = ncol(d$x), events = sum(d$delta), varnames = vn)
  class(res) <- "hdcqr"
  res
}
