## Selection checks
library(cqrfdr)

## mirror statistics and threshold
M <- c(5, 4, 3.5, 3, 2.5, 2, -1, 0.5, -0.4, 0)
stopifnot(mirror_threshold(M, 0.2, 1) == 2,
          is.infinite(mirror_threshold(c(1, -1, 0.5), 0.2, 1)),
          identical(mirror_stat(c(1, -2, 3), c(2, 1, -1)), c(3, -3, -4)))
## U_+ is never smaller than U
set.seed(1)
for (r in 1:50) {
  M <- c(rnorm(10, 3), rnorm(40))
  stopifnot(mirror_threshold(M, 0.2, 1) >= mirror_threshold(M, 0.2, 0))
}
stopifnot(all.equal(fdp_power(c(1, 2, 11), 1:10),
                    c(selected = 3, false = 1, true = 2, FDP = 1 / 3, power = 0.2)))

## simulator
set.seed(4)
d <- sim_cqr(n = 2000, p = 20, cens = 0.3, calib_n = 50000)
stopifnot(abs(d$cens_rate - 0.3) < 0.04, identical(d$active, 1:10),
          ncol(d$x) == 20, all(d$delta %in% 0:1))
h <- sim_cqr(n = 50, p = 20, s = 5, model = "hetero", tau = 0.75, cens = 0)
stopifnot(all(h$delta == 1), 6 %in% h$active)

## SS-MS end to end: strong signal, recovers the active set; reproducible
set.seed(5)
d <- sim_cqr(n = 240, p = 60, beta = 1, cens = 0.25, calib_n = 20000)
set.seed(9); a <- ssms(d$x, d$y, d$delta, q = c(0.1, 0.2), n_grid = 10, nfolds = 3)
set.seed(9); b <- ssms(d$x, d$y, d$delta, q = c(0.1, 0.2), n_grid = 10, nfolds = 3)
stopifnot(identical(a$selected, b$selected), !a$failed,
          fdp_power(selected(a, 0.2), paste0("X", d$active))[["power"]] >= 0.8,
          length(a$screen) <= a$cap,
          all(selected(a, 0.1) %in% selected(a, 0.2)))   # thresholds are nested in q
s <- summary(a); stopifnot(nrow(s) == length(a$screen), !is.unsorted(rev(s$M)))
## the "path" screen contains the CV support of the same split
set.seed(9); c0 <- ssms(d$x, d$y, d$delta, q = 0.2, n_grid = 10, nfolds = 3, screen = "cv")
stopifnot(all(c0$screen %in% a$screen))
## Surv input and fixed split
if (requireNamespace("survival", quietly = TRUE)) {
  set.seed(9)
  e <- ssms(d$x, survival::Surv(d$time, d$delta), q = 0.2, n_grid = 10,
            lambda = "pivotal", split = 1:120)
  stopifnot(identical(e$split, 1:120))
}
## input errors
stopifnot(inherits(try(ssms(d$x, d$y, d$delta[-1]), silent = TRUE), "try-error"),
          inherits(try(ssms(d$x, d$y, d$delta, q = 1.2), silent = TRUE), "try-error"),
          inherits(try(ssms(d$x, d$y, d$delta, tau = 0.05), silent = TRUE), "try-error"))

## MDS and L-HDCQR
set.seed(6); m <- mds(d$x, d$y, d$delta, B = 3, n_grid = 10, nfolds = 3)
stopifnot(!m$failed, all(abs(colSums(m$inclusion) - 1) < 1e-8 | colSums(m$inclusion) == 0))
set.seed(7); h <- hdcqr(d$x, d$y, d$delta, n_grid = 10, nfolds = 3)
cf <- coef(h)
stopifnot(length(cf) == 61, all(d$active %in% h$selected),
          max(abs(predict(h, d$x[1:3, ]) - drop(cbind(1, d$x[1:3, ]) %*% cf))) < 1e-10)
## standardized and original scales agree
cs <- coef(h, scale = "standardized")
z <- (d$x[1, ] - h$center) / h$scale
stopifnot(abs(sum(c(1, z) * cs) - sum(c(1, d$x[1, ]) * cf)) < 1e-10)
cat("selection tests passed\n")

## MS-MS: SS-MS on several splits, aggregated by inclusion rates
set.seed(5)
d <- sim_cqr(n = 240, p = 60, beta = 1, cens = 0.25, calib_n = 20000)
set.seed(3); m1 <- msms(d$x, d$y, d$delta, q = c(0.1, 0.2), B = 4, n_grid = 10, nfolds = 3)
set.seed(3); m2 <- msms(d$x, d$y, d$delta, q = c(0.1, 0.2), B = 4, n_grid = 10, nfolds = 3)
stopifnot(identical(m1$selected, m2$selected), m1$method == "MS-MS",
          inherits(m1, "cqrfdr_mds"), m1$offset == 1, m1$screen_rule == "path",
          nrow(m1$split_sizes) == 4L - m1$n_failed_splits,
          all(m1$share_empty >= 0 & m1$share_empty <= 1),
          abs(sum(m1$inclusion[, "0.2"]) - mean(m1$split_sizes[, "0.2"] > 0)) < 1e-8,
          fdp_power(selected(m1, 0.2), paste0("X", d$active))[["power"]] >= 0.8)
invisible(capture.output(print(m1), summary(m1)))

## DS of Dai et al.: CV support without a cap, cutoff U; nested in q
set.seed(5)
d <- sim_cqr(n = 240, p = 60, beta = 1, cens = 0.25, calib_n = 20000)
set.seed(4); g <- ds(d$x, d$y, d$delta, q = c(0.1, 0.2), n_grid = 10, nfolds = 3)
stopifnot(g$method == "DS", g$offset == 0, is.null(g$cap), g$screen_rule == "cv",
          all(selected(g, 0.1) %in% selected(g, 0.2)))
invisible(capture.output(print(g)))
## MDS with the design of Dai et al. (every split tuned by its own CV)
set.seed(6); h <- mds(d$x, d$y, d$delta, q = 0.2, B = 3, n_grid = 10, nfolds = 3)
stopifnot(h$tune == "cv_each", h$ties == "dai", length(h$split_mult) == 3L)

## hdcqr with the penalty sequence of Zheng, Peng and He (2018)
set.seed(11)
dz <- sim_cqr(n = 200, p = 30, s = 5, beta = 1, cens = 0.2, calib_n = 20000)
fz <- hdcqr(dz$x, dz$y, dz$delta, tau = 0.5, n_grid = 8, nfolds = 3, lambda = "zheng")
stopifnot(inherits(fz, "hdcqr"), identical(fz$tuning, "zheng"),
          length(fz$lambda) == 8L, all(diff(fz$lambda) >= 0),
          all(paste0("X", dz$active) %in% names(fz$selected)))
