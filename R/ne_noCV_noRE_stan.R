#' Bayesian negative exponential function with Stan,
#' without cross validation, without random effects
#'
#' @export
#' @param x Numeric vector of input values.
#' @param y Numberic vector of output values.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
ne_noCV_noRE_stan <- function(x, y, ...) {
  standata <- list(x = x, y = y, N = length(y))
  out <- rstan::sampling(stanmodels$ne_noCV_noRE, data = standata, ...)
  return(out)
}
