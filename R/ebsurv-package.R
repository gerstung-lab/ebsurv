#' High-dimensional extensions of the Cox proportional hazards model
#' 
#' This package includes two extensions of the Cox proportional hazards model for high-dimensional regression and tools for non-parametric 
#' simulation of survival data. The main models are a random effects model and complementary pairs stability selection.#' 
#' @name ebsurv-package
#' @docType package
#' @author Moritz Gerstung, European Bioinformatics Institute EMBL-EBI
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