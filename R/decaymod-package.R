#' The 'decaymod' package.
#'
#' @description This package contains the compiled stan code for Weibull and negative exponential decay models with options for random effects, fixed effects on each decay parameter, and cross validation.
#'
#' @docType package
#' @name decaymod-package
#' @aliases decaymod
#' @useDynLib decaymod, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.3. http://mc-stan.org
#'
NULL
