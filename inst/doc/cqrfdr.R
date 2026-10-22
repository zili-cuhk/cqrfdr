### R code from vignette source 'cqrfdr.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width = 80, digits = 4)


###################################################
### code chunk number 2: sim
###################################################
library(cqrfdr)
set.seed(2026)
d <- sim_cqr(n = 300, p = 200, s = 10, beta = 0.8, rho = 0.5, cens = 0.3,
             calib_n = 50000)
c(censoring = d$cens_rate, events = sum(d$delta))
d$beta[2:13]


###################################################
### code chunk number 3: ssms
###################################################
set.seed(1)
fit <- ssms(d$x, d$y, d$delta, tau = 0.5, q = c(0.1, 0.2), n_grid = 20)
fit


###################################################
### code chunk number 4: eval
###################################################
truth <- paste0("X", d$active)
rbind(`q = 0.1` = fdp_power(selected(fit, q = 0.1), truth),
      `q = 0.2` = fdp_power(selected(fit, q = 0.2), truth))


###################################################
### code chunk number 5: summary
###################################################
head(summary(fit), 12)


###################################################
### code chunk number 6: plot
###################################################
plot(fit)


###################################################
### code chunk number 7: msms
###################################################
set.seed(2)
fs <- msms(d$x, d$y, d$delta, tau = 0.5, q = c(0.1, 0.2), B = 10, n_grid = 20)
fs
rbind(`q = 0.1` = fdp_power(selected(fs, 0.1), truth),
      `q = 0.2` = fdp_power(selected(fs, 0.2), truth))
fs$share_empty      # share of splits whose own SS-MS selection was empty


###################################################
### code chunk number 8: comparators
###################################################
set.seed(3)
fd <- ds(d$x, d$y, d$delta, tau = 0.5, q = c(0.1, 0.2), n_grid = 20)
fdp_power(selected(fd, 0.2), truth)
fz <- hdcqr(d$x, d$y, d$delta, tau = 0.5, n_grid = 20, lambda = "zheng")
fdp_power(selected(fz), truth)
fh <- hdcqr(d$x, d$y, d$delta, tau = 0.5, n_grid = 20)
round(coef(fh, tau = c(0.25, 0.5))[1:6, ], 3)


###################################################
### code chunk number 9: mds (eval = FALSE)
###################################################
## ## each of the B splits runs its own 10-fold CV, so this takes several minutes
## fm <- mds(d$x, d$y, d$delta, tau = 0.5, q = c(0.1, 0.2), B = 50)


###################################################
### code chunk number 10: own (eval = FALSE)
###################################################
## library(survival)
## X <- as.matrix(my_data[, candidate_columns])
## fit <- msms(X, Surv(my_data$time, my_data$status), tau = 0.5, q = 0.1, cores = 4)
## selected(fit, q = 0.1)


###################################################
### code chunk number 11: plasmode (eval = FALSE)
###################################################
## plasmode_rep <- function(X, cens_rate, s = 10, b = 0.1, tau = 0.5, q = 0.1) {
##   Z <- scale(X); n <- nrow(Z)
##   act <- sort(sample(ncol(Z), s))
##   logT <- drop(Z[, act] %*% (b * sample(c(-1, 1), s, TRUE))) + rnorm(n)
##   g <- rnorm(n)
##   c0 <- uniroot(function(c) mean(c + g < logT) - cens_rate, range(logT) + c(-10, 10))$root
##   y <- pmin(logT, c0 + g); delta <- as.integer(logT <= c0 + g)
##   f <- msms(Z, y, delta, tau = tau, q = q)
##   fdp_power(selected(f, q), colnames(Z)[act])
## }


###################################################
### code chunk number 12: mc (eval = FALSE)
###################################################
## reps <- 100
## res <- t(vapply(seq_len(reps), function(r) {
##   set.seed(r)
##   dr <- sim_cqr(n = 300, p = 200, beta = 0.8, cens = 0.3, c0 = d$c0)
##   f <- ssms(dr$x, dr$y, dr$delta, q = 0.2)
##   fdp_power(selected(f, 0.2), paste0("X", dr$active))
## }, numeric(5)))
## colMeans(res[, c("FDP", "power")])   # FDR and average power


