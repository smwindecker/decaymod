#' Bayesian linear regression with Stan
#'
#' @export
#' @param x Numeric vector of input values.
#' @param y Numberic vector of output values.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
negexp_cv_re_stan <- function(x, y, ...) {
  standata <- list(x = x, y = y, N = length(y))
  out <- rstan::sampling(stanmodels$negexp_cv_re, data = standata, ...)
  return(out)
}