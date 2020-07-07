#'Empirical Bayes multi-state Cox model
#' 
#'This package implements an empirical Bayes,
#'multi-state Cox model. Different groups of regression
#'coefficients can be defined, with coefficients of the
#'same group sharing the same Gaussian prior. It takes
#'as input a data set in 'long format' and generates
#'estimates of relative hazards, cumulative hazard
#'functions and transition probabilities.
#' 
#' @name ebsurv-package
#' @docType package
#' @author Rui Costa, Moritz Gerstung
#' @details
#' \tabular{ll}{
#' Package: \tab ebsurv\cr
#' Type: \tab Package\cr
#' Version: \tab 0.0.73\cr
#' Date: \tab 2020-01-21\cr
#' License: \tab GPL 3\cr
#' }
#' @keywords package
#' @import parallel
#' @import survival
#' @import mstate
#' @import HDInterval
NA
#' @useDynLib ebsurv, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL