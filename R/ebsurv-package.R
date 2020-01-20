#'Empirical Bayes multi-state Cox model
#' 
#'This package implements an empirical Bayes,
#'multi-state Cox model.Different groups of regression
#'coefficients can be defined, with coefficients of the
#'same group sharing the same Gaussian prior. It takes
#'as input a data set in 'long format' and generates
#'estimates of relative hazards, cumulative hazard
#'functions and transition probabilities.
#' 
#' @name ebsurv-package
#' @docType package
#' @author Rui Costa, Moritz Gerstung (EMBL-EBI)
#' @keywords package
#' @import parallel
#' @import survival
#' @import RColorBrewer
#' @importFrom mice mice
#' @importFrom glmnet glmnet
#' @importFrom mvtnorm rmvnorm
#' @import MASS
NA
#' @useDynLib ebsurv, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL