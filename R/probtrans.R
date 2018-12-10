#' Compute subject-specific or overall transition probabilities
#'
#'This function computes subject-specific or overall transition
#'probabilities in multi-state models.
#'
#' @param object An object describing the fit of
#' a multi-state Cox model.
#' @param ... other arguments
#' @return An object of class "msfit"
#' @export
probtrans_mstate <- function(object, ...){
  UseMethod("probtrans_mstate")
}
#' 
#' Compute subject-specific or overall transition
#' probabilities under a random effects Cox model
#' 
#' @export
probtrans_mstate.default<-function(object, predt,
                                  direction=c("forward","fixedhorizon"),
                                  method=c("aalen","greenwood"), 
                                  variance=TRUE, covariance=FALSE){
  return(probtrans(object,predt,direction,method,variance,covariance))
}
#'
#' Compute subject-specific or overall transition
#' probabilities under a fixed effects Cox model
#'
#' @export
probtrans_mstate.coxrfx<-function(object, predt, 
                                   direction=c("forward",
                                               "fixedhorizon")){
  return(probtrans(object,predt,direction,method="aalen",variance=FALSE
                   ,covariance=FALSE))
}
