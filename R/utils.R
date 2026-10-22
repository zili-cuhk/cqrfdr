## ---------------------------------------------------------------------------
## Evaluation helper
## ---------------------------------------------------------------------------
fdp_power <- function(selected, active) {
  if (is.character(selected) || is.character(active)) {
    selected <- as.character(selected); active <- as.character(active)
  } else {
    selected <- as.integer(selected); active <- as.integer(active)
  }
  selected <- unique(selected)
  if (!length(active)) stop("'active' is empty")
  V <- length(setdiff(selected, active)); R <- length(selected)
  c(selected = R, false = V, true = R - V, FDP = V / max(R, 1),
    power = (R - V) / length(unique(active)))
}
