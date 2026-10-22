## ---------------------------------------------------------------------------
## S3 methods and accessors
## ---------------------------------------------------------------------------

selected <- function(object, ...) UseMethod("selected")

selected.cqrfdr <- function(object, q = object$q[1L], names = TRUE, ...) {
  key <- format(object$q)[match(q, object$q)]
  if (length(q) != 1L || is.na(key))
    stop("'q' must be one of the levels used in the fit: ", paste(format(object$q), collapse = ", "))
  s <- object$selected[[key]]
  if (names) names(s) else unname(s)
}

selected.hdcqr <- function(object, names = TRUE, ...) {
  if (names) names(object$selected) else unname(object$selected)
}

.show_names <- function(v, width = 0.9 * getOption("width")) {
  if (!length(v)) return("(none)")
  out <- paste(v, collapse = " ")
  if (nchar(out) > width) {
    k <- max(1L, which(cumsum(nchar(v) + 1L) <= width - 12L))
    out <- paste(paste(v[seq_len(max(k))], collapse = " "), sprintf("... (+%d)", length(v) - max(k)))
  }
  strwrap(out, width = width, exdent = 2L)
}

print.cqrfdr <- function(x, ...) {
  lab <- if (x$method == "SS-MS") "single split, mirror statistics" else
    if (x$method == "DS") "single data splitting of Dai et al., cutoff U" else
    if (x$method == "MS-MS") sprintf("SS-MS on B = %d splits, aggregated by inclusion rates", x$B) else
    sprintf("multiple data splitting, B = %d", x$B)
  cat(sprintf("%s: %s\n", x$method, lab))
  cat(sprintf("  censored quantile regression at tau = %g; n = %d (%d events), p = %d\n",
              x$tau, x$n, x$events, x$p))
  cat(sprintf("  grid of %d quantile levels on [%g, %g]\n", length(x$taus), x$taus[1L], x$tau))
  if (x$method %in% c("SS-MS", "DS")) {
    if (x$failed) {
      cat("  status:", x$status, "- no selection\n"); return(invisible(x))
    }
    cat(sprintf("  screen: %d candidates (rule \"%s\", cap %s); penalty multiplier %.3g",
                length(x$screen), x$screen_rule, if (is.null(x$cap)) "none" else x$cap,
                x$tuning$mult))
    if (!is.null(x$tuning$mult_screen) && !isTRUE(all.equal(x$tuning$mult_screen, x$tuning$mult)))
      cat(sprintf(" (screen extended to %.3g)", x$tuning$mult_screen))
    cat("\n")
    tab <- data.frame(q = x$q, selected = lengths(x$selected),
                      threshold = signif(x$threshold, 4))
  } else {
    cat(sprintf("  penalty multiplier %.3g (%s); %d of %d splits failed\n", x$mult,
                x$tune, x$n_failed_splits, x$B))
    tab <- data.frame(q = x$q, selected = lengths(x$selected))
    if (!is.null(x$share_empty)) tab$empty_splits <- sprintf("%.0f%%", 100 * x$share_empty)
  }
  cat("\n"); print(tab, row.names = FALSE)
  for (j in seq_along(x$q)) {
    cat(sprintf("\nSelected at q = %g:\n", x$q[j]))
    cat(paste0("  ", .show_names(names(x$selected[[j]]))), sep = "\n")
  }
  invisible(x)
}

summary.cqrfdr <- function(object, ...) {
  vn <- object$varnames
  if (object$method %in% c("SS-MS", "DS")) {
    s <- unname(object$screen)
    if (!length(s)) return(data.frame())
    out <- data.frame(variable = vn[s], index = s, b1 = unname(object$b1[s]),
                      b2 = unname(object$b2[s]), M = unname(object$M[s]),
                      stringsAsFactors = FALSE)
    for (j in seq_along(object$q))
      out[[paste0("sel_q", format(object$q[j]))]] <- s %in% object$selected[[j]]
    out <- out[order(-out$M), ]
  } else {
    I <- object$inclusion
    keep <- which(rowSums(I) > 0)
    if (!length(keep)) return(data.frame())
    out <- data.frame(variable = vn[keep], index = keep, stringsAsFactors = FALSE)
    for (j in seq_along(object$q)) {
      out[[paste0("incl_q", format(object$q[j]))]] <- I[keep, j]
      out[[paste0("sel_q", format(object$q[j]))]] <- keep %in% object$selected[[j]]
    }
    out <- out[order(-out[[3L]]), ]
  }
  rownames(out) <- NULL
  out
}

plot.cqrfdr <- function(x, q = x$q, label = TRUE, ...) {
  op <- graphics::par(mar = c(4.2, 4.2, 1.5, 1)); on.exit(graphics::par(op))
  if (x$method %in% c("SS-MS", "DS")) {
    s <- unname(x$screen)
    if (!length(s)) stop("empty screen: nothing to plot")
    M <- x$M[s]; o <- order(-M); M <- M[o]
    yl <- range(c(0, M)); if (label) yl[2L] <- yl[2L] + 0.22 * diff(yl)
    graphics::plot(seq_along(M), M, type = "h", lwd = 2, xlab = "Screened covariates (ordered)",
                   ylab = "Mirror statistic", col = ifelse(M > 0, "grey20", "grey60"),
                   ylim = yl, ...)
    graphics::abline(h = 0, col = "grey70")
    lt <- c(2, 3, 4, 5)
    for (j in seq_along(q)) {
      k <- match(q[j], x$q)
      if (is.finite(x$threshold[k])) graphics::abline(h = x$threshold[k], lty = lt[(j - 1L) %% 4L + 1L])
    }
    fin <- is.finite(x$threshold[match(q, x$q)])
    if (any(fin)) graphics::legend("topright", bty = "n", cex = 0.85, lty = lt[seq_along(q)][fin],
                                   legend = sprintf("threshold, q = %g", q[fin]))
    if (label) {
      top <- seq_len(min(10L, length(M)))
      graphics::text(top, M[top], names(M)[top], srt = 90, adj = c(-0.15, 0.5), cex = 0.7)
    }
  } else {
    k <- match(q[1L], x$q)
    I <- x$inclusion[, k]; o <- order(-I); I <- I[o]; keep <- I > 0
    if (!any(keep)) stop("no covariate was ever selected")
    I <- I[keep]
    yl <- c(0, max(I) * if (label) 1.22 else 1)
    graphics::plot(seq_along(I), I, type = "h", lwd = 2, xlab = "Covariates (ordered)", ylim = yl,
                   ylab = sprintf("Inclusion rate (q = %g)", q[1L]),
                   col = ifelse(names(I) %in% names(x$selected[[k]]), "grey20", "grey60"), ...)
    if (label) {
      top <- seq_len(min(10L, length(I)))
      graphics::text(top, I[top], names(I)[top], srt = 90, adj = c(-0.15, 0.5), cex = 0.7)
    }
  }
  invisible(x)
}

print.hdcqr <- function(x, ...) {
  cat(sprintf("L-HDCQR fit at tau = %g: n = %d (%d events), p = %d, %d grid points on [%g, %g]\n",
              x$tau, x$n, x$events, x$p, length(x$taus), x$taus[1L], x$tau))
  if (identical(x$tuning, "zheng"))
    cat(sprintf("  penalty of Zheng, Peng and He (2018): %.3g at tau_0, %.3g at tau; %d nonzero coefficients at tau\n",
                x$lambda[1L], x$lambda[length(x$lambda)], length(x$selected)))
  else
    cat(sprintf("  penalty multiplier %.3g; %d nonzero coefficients at tau\n",
                x$mult, length(x$selected)))
  cat(paste0("  ", .show_names(names(x$selected))), sep = "\n")
  invisible(x)
}

coef.hdcqr <- function(object, tau = object$tau, scale = c("original", "standardized"), ...) {
  scale <- match.arg(scale)
  k <- vapply(tau, function(t) which.min(abs(object$taus - t)), 1L)
  B <- object$coef_std[, k, drop = FALSE]
  if (scale == "original") {
    sl <- B[-1L, , drop = FALSE] / object$scale
    B <- rbind(B[1L, ] - colSums(sl * object$center), sl)
    rownames(B)[1L] <- "(Intercept)"
  }
  colnames(B) <- format(object$taus[k], digits = 4)
  if (length(tau) == 1L) B[, 1L] else B
}

predict.hdcqr <- function(object, newx, tau = object$tau, ...) {
  newx <- as.matrix(newx)
  if (ncol(newx) != object$p) stop("'newx' must have ", object$p, " columns")
  B <- coef(object, tau = tau)
  B <- as.matrix(B)
  drop(cbind(1, newx) %*% B)
}
