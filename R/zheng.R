## L-HDCQR as specified by Zheng, Peng and He (2018, Section 5.1) --------------
## Sequential L1-penalized Peng-Huang estimation on the grid tau_0 < ... < tau_L
## with a penalty that grows along the grid:
##   * lambda_0: pivotal choice at tau_0 (Belloni and Chernozhukov, 2011);
##   * k = 1, ..., L: lambda_k = lambda_{k-1} + increment, the increment chosen
##     by K-fold CV from `increments` * sqrt(n log(p v n)) * eps_n, eps_n the
##     grid spacing; the held-out criterion at tau_k is the Peng-Huang
##     objective of the held-out data (Zheng et al. do not state the criterion);
##   * the compensator starts at tau_0, as in their equation (7);
##   * each step is solved as the augmented median (LAD) problem of their
##     Section 3.1.
## Used by hdcqr(lambda = "zheng").

.z18_state <- function(y, delta, Z, rows, a0) {
  ev <- rows[delta[rows] == 1]
  list(rows = rows, Ze = Z[ev, , drop = FALSE], ye = y[ev],
       zminus = -colSums(Z[ev, , drop = FALSE]), a = rep(a0, length(rows)), b = NULL)
}

.z18_solve <- function(st, Z, lam, R0 = 1e4, R_growth = 10, R_tries = 5L) {
  P <- ncol(Z)
  zplus <- 2 * colSums(Z[st$rows, , drop = FALSE] * st$a)
  pen <- c(0, rep(lam, P - 1L))
  init <- if (is.null(st$b)) integer(0) else which(st$b != 0)
  R <- R0
  for (t in seq_len(R_tries)) {
    fit <- tryCatch(
      if (P > 100L) lad_lasso_ws(rbind(st$Ze, st$zminus, zplus), c(st$ye, R, R), pen, init = init)
      else lad_lasso(rbind(st$Ze, st$zminus, zplus), c(st$ye, R, R), pen),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    b <- fit$coef
    if (all(R - c(sum(st$zminus * b), sum(zplus * b)) > 0)) return(b)
    R <- R * R_growth
  }
  NULL
}

.z18_advance <- function(st, y, Z, b, dH) {
  st$a <- st$a + at_risk(y[st$rows], drop(Z[st$rows, , drop = FALSE] %*% b)) * dH
  st$b <- b
  st
}

## Z: n x (p + 1) standardized design with intercept column first.
## Returns list(status, B, lambda, increment) or list(status, k) on failure.
.hdcqr_zheng <- function(y, delta, Z, taus, K = 5L, foldid = NULL,
                         increments = c(0, exp(seq(log(0.02), log(20), length.out = 12))),
                         nsim = 500L, cores = 1L) {
  n <- length(y); m <- length(taus); p <- ncol(Z) - 1L
  fold <- if (is.null(foldid)) make_folds(delta, K) else match(foldid, sort(unique(foldid)))
  K <- max(fold)
  dH <- diff(H(taus)); a0 <- taus[1L]
  unit <- sqrt(n * log(max(p, n))) * (if (m > 1L) taus[2L] - taus[1L] else 0)
  lam <- numeric(m)
  lam[1L] <- pivotal_lambda(Z[, -1L, drop = FALSE], taus[1L], nsim)
  full <- .z18_state(y, delta, Z, seq_len(n), a0)
  tr <- lapply(seq_len(K), function(k) .z18_state(y, delta, Z, which(fold != k), a0))
  te_rows <- lapply(seq_len(K), function(k) which(fold == k))
  te_comp <- lapply(te_rows, function(r) rep(a0, length(r)))
  ho_err <- function(b, kk) {                     # held-out Peng-Huang objective of fold kk
    r <- te_rows[[kk]]; fit <- drop(Z[r, , drop = FALSE] %*% b)
    sum(delta[r] * abs(y[r] - fit) + (delta[r] - 2 * te_comp[[kk]]) * fit)
  }
  B <- matrix(NA_real_, ncol(Z), m)
  choice <- rep(NA_real_, m)
  for (k in seq_len(m)) {
    if (k == 1L) {
      cand <- lam[1L]; sols <- lapply(tr, function(st) list(.z18_solve(st, Z, cand)))
      best <- 1L
    } else {
      cand <- lam[k - 1L] + increments * unit
      jobs <- expand.grid(f = seq_len(K), c = seq_along(cand))
      run <- function(j) .z18_solve(tr[[jobs$f[j]]], Z, cand[jobs$c[j]])
      res <- if (cores > 1L && .Platform$OS.type != "windows")
               parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = cores, mc.set.seed = FALSE)
             else lapply(seq_len(nrow(jobs)), run)
      err <- vapply(seq_along(cand), function(cc) {
        e <- vapply(seq_len(K), function(f) {
          b <- res[[which(jobs$f == f & jobs$c == cc)]]
          if (is.null(b)) Inf else ho_err(b, f) }, 0)
        sum(e) }, 0)
      if (all(!is.finite(err))) return(list(status = "no_root", k = k))
      best <- which(err == min(err))[1L]          # smallest increment on ties
      sols <- lapply(seq_len(K), function(f) list(res[[which(jobs$f == f & jobs$c == best)]]))
      choice[k] <- increments[best]
    }
    lam[k] <- cand[best]
    bf <- .z18_solve(full, Z, lam[k])
    if (is.null(bf)) return(list(status = "no_root", k = k))
    B[, k] <- bf
    if (k < m) {
      full <- .z18_advance(full, y, Z, bf, dH[k])
      for (f in seq_len(K)) {
        bt <- sols[[f]][[1L]]
        if (is.null(bt)) bt <- bf                 # failed fold: continue from the full fit
        te_comp[[f]] <- te_comp[[f]] + at_risk(y[te_rows[[f]]], drop(Z[te_rows[[f]], , drop = FALSE] %*% bt)) * dH[k]
        tr[[f]] <- .z18_advance(tr[[f]], y, Z, bt, dH[k])
      }
    }
  }
  list(status = "ok", B = B, lambda = lam, increment = choice)
}
