
#' Empirical Bayes, multi-state Cox model
#' 
#' This function estimates a multi-state Cox model with one or more Gaussian priors
#' imposed on the regression coefficients (see Therneau et al., 2003).
#' Multiple groups of coefficients can be defined: coefficients within a group share 
#' the same (possibly unknown) mean and variance. The parameters and hyperparameters are
#' efficiently estimated by an EM-type algorithm.
#' @param Z A data frame consisting of the covariate columns of a data set in 'long format'.
#' @param surv A `survival' object created with \code{survival::Surv}.
#' @param groups A character or numeric vector whose \eqn{i}th element gives the group of the regression
#' coefficient associated with the \eqn{i}th covariate column of Z (coefficients belonging to the same group 
#' share the same Gaussian prior).
#' @param which.mu A vector with names or numbers of coefficient groups (see 
#' argument \code{groups}). If the name or number of a group of coefficients is
#' given in this argument, \code{CoxRFX} will estimate the mean of its Gaussian distribution;
#' otherwise the mean will be fixed at zero.
#' @param tol Convergence criterium of the EM algorithm. The algorithm stops unless
#'  there is at least one parameter (or hyperparameter) for which it holds that the
#'  current estimate differs in absolute terms by more than \code{tol} from the
#'  previous estimate. 
#' @param max.iter The maximum number of iterations in the EM algorithm.
#' @param sigma0 A vector with the initial value of the variance hyperparameter for each group of coefficients.
#' Or a single value, in case the initial value of the variance hyperparameter is meant to be the same for all groups.
#' @param sigma.hat Which estimator to use for the variance hyperparameters (see details).
#' @param verbose Gives more output.
#' @details The argument \code{Z} must be of class \code{c(data.frame,msdata)}. 
#' 
#' Different estimators exist for the variance hyperparameters: the default is "df", as used by Perperoglou (2014) and introduced by Schall (1991). 
#' Alternatives are MLE, REML, and BLUP, as defined by Therneau et al. (2003). 
#' Simulations suggest that the 'df' method is the most accurate.
#' 
#' The model can also be fitted using package \code{coxme}; the \code{coxme}
#' routine numerically optimises the integrated partial likelihood, which may
#' be more accurate, but is computationally expensive.
#' 
#' @references Terry M Therneau, Patricia M Grambsch & V. Shane Pankratz (2003) Penalized Survival Models and Frailty, Journal of Computational and Graphical Statistics, 12:1, 156-175, http://dx.doi.org/10.1198/1061860031365
#' 
#' A. Perperoglou (2014). Cox models with dynamic ridge penalties on time-varying effects of the covariates. Stat Med, 33:170-80. http://dx.doi.org/10.1002/sim.5921
#' 
#' R. Schall (1991). Estimation in generalized linear models with random effects. Biometrika, 78:719-727. http://dx.doi.org/10.1093/biomet/78.4.719

#' @return A coxph object (see \code{survival::coxph.object}) with a few extra fields: the inputs $groups, $Z, and $surv;
#' and the hyperparameters $sigma2 (variances) and $mu (means). 
#' 
#' @author Moritz Gerstung & Rui Costa
#' @seealso \code{\link{coxph.object}}; \code{\link{coxme}}; \code{\link{Surv}}.
#' @export
# @example inst/example/CoxRFX-example.R
CoxRFX <- function(Z, surv, groups = rep(1, ncol(Z)), which.mu = unique(groups), tol=1e-3, max.iter=50, sigma0 = 0.1, sigma.hat=c("df","MLE","REML","BLUP"), verbose=FALSE, ...){
  ##
  strata<-Z$strata
  Z$strata<-NULL
  namesZ<-names(Z)
  ##
  Z = as.matrix(Z)
	Z.df <- TRUE
	if(is.null(colnames(Z)))
		colnames(Z) <- make.names(1:ncol(Z))
	sigma.hat = match.arg(sigma.hat)
	o <- order(groups)
	Z <- Z[,o]
	groups <- factor(groups[o])
	uniqueGroups <- levels(groups)
	ZZ <- lapply(uniqueGroups, function(i) Z[,groups==i, drop=FALSE])
	names(ZZ) <- uniqueGroups
	sumZ <- sapply(which.mu, function(i) rowSums(ZZ[[i]]))
	nGroups = length(uniqueGroups)
	sigma2 <- sigma0ld <- rep(ifelse(sigma0>0, sigma0,1), nGroups)
	iter = 1
	mu <- mu0ld <- rep(0, nGroups)
	names(mu) <- uniqueGroups
	beta = rep(1,ncol(Z)+length(which.mu))
	beta0ld = rep(0,ncol(Z)+length(which.mu))
	sigma2.mu = sigma0
	nu=0
	penalize.mu = FALSE
	if(!is.null(which.mu)) 
		if(!penalize.mu)
			sumTerm <- "sumZ" 
		else
			sumTerm <- "ridge(sumZ, theta=1/sigma2.mu, scale=FALSE)"
	else sumTerm <- character(0)
	while((max(abs(beta-beta0ld)) > tol | max(abs(mu - mu0ld)) > tol | max(abs(sigma2 - sigma0ld)) > tol) & iter < max.iter){
		beta0ld = beta
		sigma0ld <- sigma2
		mu0ld <- mu
		formula <- formula(paste("surv ~", paste(c(sapply(1:nGroups, function(i) paste("ridge(ZZ[[",i,"]], theta=1/sigma2[",i,"], scale=FALSE)", sep="")), 
								sumTerm,"strata(strata)"), 
						collapse=" + ")))
		fit <- ebsurv:::coxph(formula, ...)
		if(any(is.na(coef(fit)))){
			warning(paste("NA during estimation (iter: ", iter, ", coef: ", paste(which(is.na(coef(fit)[order(o)])), sep=","), ")", sep=""))
			break
		}
		if(!is.null(which.mu))
			mu[which.mu] <- coef(fit)[-(1:ncol(Z))]
		if(verbose) cat("mu", mu, "\n", sep="\t")
		names(fit$df) <- c(uniqueGroups, rep("Offset", length(which.mu)>0))
		if(verbose) cat("df", fit$df,"\n", sep="\t")
		sigma2 = sapply(uniqueGroups, function(i){
					index <- which(groups==i) #& fit$coefficients > beta.thresh
					if(sigma.hat=="BLUP")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ))/(nu + length(index))
					else if(sigma.hat=="df")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ))/(nu + fit$df[i])
					else if(sigma.hat == "MLE")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ) + sum(diag(solve(solve(fit$var)[index,index]))))/(nu + length(index))
					else if(sigma.hat == "REML")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ) + sum(diag(fit$var)[index]))/(nu + length(index))
				})
		if(verbose) {
			cat("sigma2", sigma2, "\n", sep="\t")
			cat("loglik:", fit$loglik - c(0,fit$penalty[2] + 1/2 * sum(log(sigma2[groups]))),"\n", sep="\t")
		}
		if(penalize.mu){
			if(sigma.hat=="BLUP")
				sigma2.mu = (sigma0 * nu + sum((mu-0)^2)) / (nu + length(mu))
			else if(sigma.hat=="df")
				sigma2.mu = (sigma0 * nu + sum((mu-0)^2)) / (nu + fit$df["Offset"])
			else if(sigma.hat == "MLE")
				sigma2.mu = (nu * sigma0 + sum((mu-0)^2 ) + sum(diag(solve(solve(fit$var)[-(1:ncol(Z)),-(1:ncol(Z))]))))/(nu + length(mu))
			else if(sigma.hat == "REML")
				sigma2.mu = (nu * sigma0 + sum((mu-0)^2 ) + sum(diag(fit$var)[-(1:ncol(Z))]))/(nu + length(mu))
		}
		
		beta = fit$coefficients
		iter = iter+1
	}
	if(iter == max.iter)
		warning("Did not converge after ", max.iter, " iterations.")
	fit$iter[1] <- iter
	fit$sigma2 = sigma0ld
	names(fit$sigma2) <- uniqueGroups
	#fit$sigma2.mu = sigma2.mu
	fit$mu = mu
	fit$Z = Z[,order(o)]
	fit$surv = surv
	C <- rbind(diag(1, ncol(Z)),t(as.matrix(MakeInteger(groups)[which.mu]))) ## map from centred to uncentred coefficients 
	fit$groups = groups[order(o)]
	var = fit$var
	var2 = fit$var2
	colnames(var) <- rownames(var) <- colnames(var2) <- rownames(var2) <- rownames(C) <- c(colnames(Z), which.mu)
	colnames(C) <- colnames(Z)
	p <- ncol(Z)
	i <- c(order(o), (1:ncol(var))[-p:-1])
	j <- order(o)
	fit$C <- C[i,j]
	fit$Hinv <- var[i,i] ## Hinv 
	fit$V <- var2[i,i] ## Hinv I Hinv
	fit$z <- (fit$coefficients / sqrt(diag(var)))[i] ## z-scores of centred coefficients
	fit$z2 <- (fit$coefficients / sqrt(diag(var2)))[i] ## z-scores of centred coefficients (var2)
	fit$var = (t(C) %*% var %*% C)[j,j] ## covariance of uncentred coef
	fit$var2 = (t(C) %*% var2 %*% C)[j,j] ## covariance of uncentred coef (var2)
	fit$mu.var = var[-(1:p),-(1:p)] ## covariance of mean
	fit$mu.var2 = var2[-(1:p),-(1:p)] ## covariance of mean (var2)
	fit$means = fit$means[1:p][j]
	fit$coefficients <- (fit$coefficients %*% C)[j]
	names(fit$means) <- names(fit$coefficients) <- namesZ
	fit$terms <- fit$terms[1:length(uniqueGroups)]
	fit$penalized.loglik <- fit$loglik[2] - fit$penalty[2] - 1/2 * sum(log(fit$sigma2[groups]))
	## Fake call for predict.coxph and survfit.coxph
	call <- match.call()
	if(Z.df){
		call["data"] <- call["Z"]
		formula <- as.formula(paste(as.character(call["surv"]),"~",paste(colnames(Z)[j], collapse="+"),"+strata(strata)"))
	}else{
		formula <- as.formula(paste(as.character(call["surv"]),"~",as.character(call["Z"])))
	}
	attr(formula,".Environment") <- parent.frame()
	fit$formula <- formula
	call["formula"] <- call("foo",formula=formula)["formula"]
	fit$terms <- terms(formula)
	attr(fit$terms,"specials")<-list(strata=NULL,cluster=NULL)
	attr(fit$terms,"specials")$strata<-length(attr(fit$terms,"term.labels"))+1
	fit$call <- call
	fit$strata<-strata
	class(fit) <- c("coxrfx", class(fit))
	return(fit)
}


#' A summary method for CoxRFX models
#' 
#' This function prints the point estimates of parameters and 
#' hyperparameters contained in a \code{coxrfx} object.
#' 
#' 
#' @param object A \code{coxrfx} object 
#' (obtained by running the function \code{CoxRFX}).
#' @return NULL
#' 
#' @author Rui Costa
#' @export
#' @method summary coxrfx
summary.coxrfx <- function(object){
	cat("Hyperparameters:\n")
	show(format(data.frame(mean=object$mu, variance=object$sigma2), digits=2))
	cat("\n Parameters: \n")
	show(format(data.frame(group=object$groups, estimate=object$coefficients,row.names = names(object$coefficients)), digits=2))
	}

#' Print method for CoxRFX objects
#' 
#' This function implicitly calls summary.coxrfx().
#' @param x A \code{coxrfx} object
#' @return NULL
#' 
#' @author Moritz Gerstung & Rui Costa
#' @export
#' @method print coxrfx

print.coxrfx <- function(x){
	summary.coxrfx(x)
}

#' Print method for \code{msfit} objects
#'  generated by \code{msfit_generic}
#' 
#' This method is a simple call to \code{print.default}.
#' Its main purpose is to override \code{print.coxrfx}
#' when printing an
#' object of double class \code{msfit} and \code{coxrfx}.
#' @param x An object of class \code{msfit} or double class \code{msfit} 
#' and \code{coxrfx}.
#' @return NULL
#' 
#' @author Rui Costa
#' @export
#' @method print msfit

print.msfit <- function(x){
  #mstate:::summary.msfit(x)
  print.default(x)
}


#' Convert factor to integer.
#' @param F A factor vector.
#' @return A data.frame with columns corresponding to levels in the factor. 
#' @details An internal function of \code{CoxRFX}, not meant
#' to called directly by the user.
#' @author Moritz Gerstung
#' @seealso \code{\link{CoxRFX}}
MakeInteger <- function(F){
  res <- as.data.frame(lapply(levels(F), `==`, F))
  colnames(res) <- levels(F)
  res + 0
}

coxph<-function (formula, data, weights, subset, na.action, init, control, 
          ties = c("efron", "breslow", "exact"), singular.ok = TRUE, 
          robust, model = FALSE, x = FALSE, y = TRUE, tt, method = ties, 
          id, cluster, istate, statedata, ...) 
{
  ties <- match.arg(ties)
  Call <- match.call()
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(coxph.control))
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
    if (any(indx == 0L)) 
      stop(gettextf("Argument %s not matched", names(extraArgs)[indx == 
                                                                  0L]), domain = NA)
  }
  if (missing(control)) 
    control <- coxph.control(...)
  if (missing(formula)) 
    stop("a formula argument is required")
  ss <- c("cluster", "offset")
  if (is.list(formula)) 
    Terms <- if (missing(data)) 
      terms(formula[[1]], specials = ss)
  else terms(formula[[1]], specials = ss, data = data)
  else Terms <- if (missing(data)) 
    terms(formula, specials = ss)
  else terms(formula, specials = ss, data = data)
  tcl <- attr(Terms, "specials")$cluster
  if (length(tcl) > 1) 
    stop("a formula cannot have multiple cluster terms")
  if (length(tcl) > 0) {
    factors <- attr(Terms, "factors")
    if (any(factors[tcl, ] > 1)) 
      stop("cluster() cannot be in an interaction")
    if (attr(Terms, "response") == 0) 
      stop("formula must have a Surv response")
    temp <- attr(Terms, "term.labels")
    oo <- attr(Terms, "specials")$offset
    if (!is.null(oo)) {
      ooterm <- rownames(factors)[oo]
      if (oo < tcl) 
        temp <- c(ooterm, temp)
      else temp <- c(temp, ooterm)
    }
    if (is.null(Call$cluster)) 
      Call$cluster <- attr(Terms, "variables")[[1 + tcl]][[2]]
    else warning("cluster appears both in a formula and as an argument, formula term ignored")
    if (is.list(formula)) 
      formula[[1]][[3]] <- reformulate(temp[1 - tcl])[[2]]
    else formula[[3]] <- reformulate(temp[1 - tcl])[[2]]
    Call$formula <- formula
  }
  indx <- match(c("formula", "data", "weights", "subset", 
                  "na.action", "cluster", "id", "istate"), names(Call), 
                nomatch = 0)
  if (indx[1] == 0) 
    stop("A formula argument is required")
  tform <- Call[c(1, indx)]
  tform[[1L]] <- quote(stats::model.frame)
  if (is.list(formula)) {
    multiform <- TRUE
    dformula <- formula[[1]]
    if (missing(statedata)) 
      covlist <- parsecovar1(formula[-1])
    else {
      if (!inherits(statedata, "data.frame")) 
        stop("statedata must be a data frame")
      if (is.null(statedata$state)) 
        stop("statedata data frame must contain a 'state' variable")
      covlist <- parsecovar1(formula[-1], names(statedata))
    }
    tlab <- unlist(lapply(covlist$rhs, function(x) attr(terms.formula(x$formula), 
                                                        "term.labels")))
    tlab <- c(attr(terms.formula(dformula), "term.labels"), 
              tlab)
    newform <- reformulate(tlab, dformula[[2]])
    environment(newform) <- environment(dformula)
    formula <- newform
    tform$na.action <- na.pass
  }
  else {
    multiform <- FALSE
    covlist <- NULL
    dformula <- formula
  }
  special <- c("strata", "tt")
  tform$formula <- if (missing(data)) 
    terms(formula, special)
  else terms(formula, special, data = data)
  if (!is.null(attr(tform$formula, "specials")$tt)) {
    coxenv <- new.env(parent = environment(formula))
    assign("tt", function(x) x, envir = coxenv)
    environment(tform$formula) <- coxenv
  }
  mf <- eval(tform, parent.frame())
  if (nrow(mf) == 0) 
    stop("No (non-missing) observations")
  Terms <- terms(mf)
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
  type <- attr(Y, "type")
  multi <- FALSE
  if (type == "mright" || type == "mcounting") 
    multi <- TRUE
  else if (type != "right" && type != "counting") 
    stop(paste("Cox model doesn't support \"", type, "\" survival data", 
               sep = ""))
  data.n <- nrow(Y)
  if (!multi && multiform) 
    stop("formula is a list but the response is not multi-state")
  if (control$timefix) 
    Y <- aeqSurv(Y)
  if (length(attr(Terms, "variables")) > 2) {
    ytemp <- terms.inner(formula[1:2])
    xtemp <- terms.inner(formula[-2])
    if (any(!is.na(match(xtemp, ytemp)))) 
      warning("a variable appears on both the left and right sides of the formula")
  }
  strats <- attr(Terms, "specials")$strata
  hasinteractions <- FALSE
  dropterms <- NULL
  if (length(strats)) {
    stemp <- untangle.specials(Terms, "strata", 1)
    if (length(stemp$vars) == 1) 
      strata.keep <- mf[[stemp$vars]]
    else strata.keep <- strata(mf[, stemp$vars], shortlabel = TRUE)
    istrat <- as.integer(strata.keep)
    for (i in stemp$vars) {
      if (any(attr(Terms, "order")[attr(Terms, "factors")[i, 
                                                          ] > 0] > 1)) 
        hasinteractions <- TRUE
    }
    if (!hasinteractions) 
      dropterms <- stemp$terms
  }
  else istrat <- NULL
  if (hasinteractions && multi) 
    stop("multi-state coxph does not support strata*covariate interactions")
  timetrans <- attr(Terms, "specials")$tt
  if (missing(tt)) 
    tt <- NULL
  if (length(timetrans)) {
    if (multi) 
      stop("the tt() transform is not implemented for multi-state models")
    timetrans <- untangle.specials(Terms, "tt")
    ntrans <- length(timetrans$terms)
    if (is.null(tt)) {
      tt <- function(x, time, riskset, weights) {
        obrien <- function(x) {
          r <- rank(x)
          (r - 0.5)/(0.5 + length(r) - r)
        }
        unlist(tapply(x, riskset, obrien))
      }
    }
    if (is.function(tt)) 
      tt <- list(tt)
    if (is.list(tt)) {
      if (any(!sapply(tt, is.function))) 
        stop("The tt argument must contain function or list of functions")
      if (length(tt) != ntrans) {
        if (length(tt) == 1) {
          temp <- vector("list", ntrans)
          for (i in 1:ntrans) temp[[i]] <- tt[[1]]
          tt <- temp
        }
        else stop("Wrong length for tt argument")
      }
    }
    else stop("The tt argument must contain a function or list of functions")
    if (ncol(Y) == 2) {
      if (length(strats) == 0) {
        sorted <- order(-Y[, 1], Y[, 2])
        newstrat <- rep.int(0L, nrow(Y))
        newstrat[1] <- 1L
      }
      else {
        sorted <- order(istrat, -Y[, 1], Y[, 2])
        newstrat <- as.integer(c(1, 1 * (diff(istrat[sorted]) != 
                                           0)))
      }
      if (storage.mode(Y) != "double") 
        storage.mode(Y) <- "double"
      counts <- .Call(Ccoxcount1, Y[sorted, ], as.integer(newstrat))
      tindex <- sorted[counts$index]
    }
    else {
      if (length(strats) == 0) {
        sort.end <- order(-Y[, 2], Y[, 3])
        sort.start <- order(-Y[, 1])
        newstrat <- c(1L, rep(0, nrow(Y) - 1))
      }
      else {
        sort.end <- order(istrat, -Y[, 2], Y[, 3])
        sort.start <- order(istrat, -Y[, 1])
        newstrat <- c(1L, as.integer(diff(istrat[sort.end]) != 
                                       0))
      }
      if (storage.mode(Y) != "double") 
        storage.mode(Y) <- "double"
      counts <- .Call(Ccoxcount2, Y, as.integer(sort.start - 
                                                  1L), as.integer(sort.end - 1L), as.integer(newstrat))
      tindex <- counts$index
    }
    Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
    type <- "right"
    mf <- mf[tindex, ]
    istrat <- rep(1:length(counts$nrisk), counts$nrisk)
    weights <- model.weights(mf)
    if (!is.null(weights) && any(!is.finite(weights))) 
      stop("weights must be finite")
    tcall <- attr(Terms, "variables")[timetrans$terms + 
                                        2]
    pvars <- attr(Terms, "predvars")
    pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
    for (i in 1:ntrans) {
      newtt <- (tt[[i]])(mf[[timetrans$var[i]]], Y[, 1], 
                         istrat, weights)
      mf[[timetrans$var[i]]] <- newtt
      nclass <- class(newtt)
      if (any(nclass %in% pmethod)) {
        dummy <- as.call(list(as.name(class(newtt)[1]), 
                              tcall[[i]][[2]]))
        ptemp <- makepredictcall(newtt, dummy)
        pvars[[timetrans$terms[i] + 2]] <- ptemp
      }
    }
    attr(Terms, "predvars") <- pvars
  }
  xlevels <- .getXlevels(Terms, mf)
  cluster <- model.extract(mf, "cluster")
  id <- model.extract(mf, "id")
  if (!is.null(id) && is.null(cluster) && (missing(robust) || 
                                           robust)) 
    cluster <- id
  if (missing(robust)) 
    robust <- !is.null(cluster)
  else if (robust && is.null(cluster)) {
    if (ncol(Y) == 2) 
      cluster <- seq.int(1, nrow(mf))
    else stop("one of cluster or id is needed")
  }
  contrast.arg <- NULL
  attr(Terms, "intercept") <- 1
  id <- model.extract(mf, "id")
  if (multi) {
    if (length(dropterms)) {
      Terms2 <- Terms[-dropterms]
      dformula <- formula(Terms2)
    }
    else Terms2 <- Terms
    if (length(id) == 0) 
      stop("an id statement is required for multi-state models")
    istate <- model.extract(mf, "istate")
    mcheck <- survcheck2(Y, id, istate)
    if (mcheck$flag["overlap"] > 0) 
      stop("data set has overlapping intervals for one or more subjects")
    transitions <- mcheck$transitions
    istate <- mcheck$istate
    states <- mcheck$states
    if (missing(statedata)) 
      covlist2 <- parsecovar2(covlist, NULL, dformula = dformula, 
                              Terms2, transitions, states)
    else covlist2 <- parsecovar2(covlist, statedata, dformula = dformula, 
                                 Terms2, transitions, states)
    tmap <- covlist2$tmap
    if (!is.null(covlist)) {
      good.tran <- bad.tran <- rep(FALSE, nrow(Y))
      termname <- rownames(attr(Terms, "factors"))
      trow <- (!is.na(match(rownames(tmap), termname)))
      termiss <- matrix(0L, nrow(mf), ncol(mf))
      for (i in 1:ncol(mf)) {
        xx <- is.na(mf[[i]])
        if (is.matrix(xx)) 
          termiss[, i] <- apply(xx, 1, any)
        else termiss[, i] <- xx
      }
      for (i in levels(istate)) {
        rindex <- which(istate == i)
        j <- which(covlist2$mapid[, 1] == match(i, states))
        for (jcol in j) {
          k <- which(trow & tmap[, jcol] > 0)
          bad.tran[rindex] <- (bad.tran[rindex] | apply(termiss[rindex, 
                                                                k, drop = FALSE], 1, any))
          good.tran[rindex] <- (good.tran[rindex] | 
                                  apply(!termiss[rindex, k, drop = FALSE], 
                                        1, all))
        }
      }
      n.partially.used <- sum(good.tran & bad.tran & !is.na(Y))
      omit <- (!good.tran & bad.tran) | is.na(Y)
      if (all(omit)) 
        stop("all observations deleted due to missing values")
      temp <- setNames(seq(omit)[omit], attr(mf, "row.names")[omit])
      attr(temp, "class") <- "omit"
      mf <- mf[!omit, , drop = FALSE]
      attr(mf, "na.action") <- temp
      Y <- Y[!omit]
      id <- id[!omit]
      if (length(istate)) 
        istate <- istate[!omit]
    }
  }
  if (length(dropterms)) {
    Terms2 <- Terms[-dropterms]
    X <- model.matrix(Terms2, mf, constrasts = contrast.arg)
    temp <- attr(X, "assign")
    shift <- sort(dropterms)
    temp <- temp + 1 * (shift[1] <= temp)
    if (length(shift) == 2) 
      temp + 1 * (shift[2] <= temp)
    attr(X, "assign") <- temp
  }
  else X <- model.matrix(Terms, mf, contrasts = contrast.arg)
  Xatt <- attributes(X)
  if (hasinteractions) 
    adrop <- c(0, untangle.specials(Terms, "strata")$terms)
  else adrop <- 0
  xdrop <- Xatt$assign %in% adrop
  X <- X[, !xdrop, drop = FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts
  offset <- model.offset(mf)
  if (is.null(offset) | all(offset == 0)) 
    offset <- rep(0, nrow(mf))
  else if (any(!is.finite(offset) | !is.finite(exp(offset)))) 
    stop("offsets must lead to a finite risk score")
  weights <- model.weights(mf)
  if (!is.null(weights) && any(!is.finite(weights))) 
    stop("weights must be finite")
  assign <- attrassign(X, Terms)
  contr.save <- attr(X, "contrasts")
  if (sum(Y[, ncol(Y)]) == 0) {
    ncoef <- ncol(X)
    ctemp <- rep(NA, ncoef)
    names(ctemp) <- colnames(X)
    concordance = c(concordant = 0, discordant = 0, tied.x = 0, 
                    tied.y = 0, tied.xy = 0, concordance = NA, std = NA, 
                    timefix = FALSE)
    rval <- list(coefficients = ctemp, var = matrix(0, ncoef, 
                                                    ncoef), loglik = c(0, 0), score = 0, iter = 0, linear.predictors = offset, 
                 residuals = rep(0, data.n), means = colMeans(X), 
                 method = method, n = data.n, nevent = 0, terms = Terms, 
                 assign = assign, concordance = concordance, wald.test = 0, 
                 y = Y, call = Call)
    class(rval) <- "coxph"
    return(rval)
  }
  if (multi) {
    cmap <- parsecovar3(tmap, colnames(X), attr(X, "assign"))
    xstack <- stacker(cmap, as.integer(istate), X, Y, strata = istrat, 
                      states = states)
    rkeep <- unique(xstack$rindex)
    transitions <- survcheck2(Y[rkeep, ], id[rkeep], istate[rkeep])$transitions
    X <- xstack$X
    Y <- xstack$Y
    istrat <- xstack$strata
    if (length(offset)) 
      offset <- offset[xstack$rindex]
    if (length(weights)) 
      weights <- weights[xstack$rindex]
    if (length(cluster)) 
      cluster <- cluster[xstack$rindex]
    t2 <- tmap[-1, , drop = FALSE]
    r2 <- row(t2)[!duplicated(as.vector(t2))]
    c2 <- col(t2)[!duplicated(as.vector(t2))]
    a2 <- lapply(seq(along = r2), function(i) {
      cmap[1 + assign[[r2[i]]], c2[i]]
    })
    tab <- table(r2)
    count <- tab[r2]
    names(a2) <- ifelse(count == 1, row.names(t2)[r2], paste(row.names(t2)[r2], 
                                                             colnames(cmap)[c2], sep = "_"))
    assign <- a2
  }
  if (!all(is.finite(X))) 
    stop("data contains an infinite predictor")
  if (missing(init)) 
    init <- NULL
  else {
    if (length(init) != ncol(X)) 
      stop("wrong length for init argument")
    temp <- X %*% init - sum(colMeans(X) * init) + offset
    if (any(exp(temp) > .Machine$double.xmax) || all(exp(temp) == 
                                                     0)) 
      stop("initial values lead to overflow or underflow of the exp function")
  }
  pterms <- sapply(mf, inherits, "coxph.penalty")
  if (any(pterms)) {
    pattr <- lapply(mf[pterms], attributes)
    pname <- names(pterms)[pterms]
    ord <- attr(Terms, "order")[match(pname, attr(Terms, 
                                                  "term.labels"))]
    if (any(ord > 1)) 
      stop("Penalty terms cannot be in an interaction")
    pcols <- assign[match(pname, names(assign))]
    fit <- coxpenal.fit(X, Y, istrat, offset, init = init, 
                        control, weights = weights, method = method, row.names(mf), 
                        pcols, pattr, assign)
  }
  else {
    if (method == "breslow" || method == "efron") {
      if (grepl("right", type)) 
        fitter <- get("coxph.fit")
      else fitter <- get("agreg.fit")
    }
    else if (method == "exact") {
      if (type == "right") 
        fitter <- get("coxexact.fit")
      else fitter <- get("agexact.fit")
    }
    else stop(paste("Unknown method", method))
    fit <- fitter(X, Y, istrat, offset, init, control, weights = weights, 
                  method = method, row.names(mf))
  }
  if (is.character(fit)) {
    fit <- list(fail = fit)
    class(fit) <- "coxph"
  }
  else {
    if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
      vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
      msg <- paste("X matrix deemed to be singular; variable", 
                   paste(vars, collapse = " "))
      if (!singular.ok) 
        stop(msg)
    }
    fit$n <- data.n
    fit$nevent <- sum(Y[, ncol(Y)])
    fit$terms <- Terms
    fit$assign <- assign
    class(fit) <- fit$class
    fit$class <- NULL
    if (robust && !is.null(fit$coefficients) && !all(is.na(fit$coefficients))) {
      fit$naive.var <- fit$var
      fit2 <- c(fit, list(x = X, y = Y, weights = weights))
      if (length(istrat)) 
        fit2$strata <- istrat
      if (length(cluster)) {
        temp <- residuals.coxph(fit2, type = "dfbeta", 
                                collapse = cluster, weighted = TRUE)
        if (is.null(init)) 
          fit2$linear.predictors <- 0 * fit$linear.predictors
        else fit2$linear.predictors <- c(X %*% init)
        temp0 <- residuals.coxph(fit2, type = "score", 
                                 collapse = cluster, weighted = TRUE)
      }
      else {
        temp <- residuals.coxph(fit2, type = "dfbeta", 
                                weighted = TRUE)
        fit2$linear.predictors <- 0 * fit$linear.predictors
        temp0 <- residuals.coxph(fit2, type = "score", 
                                 weighted = TRUE)
      }
      fit$var <- t(temp) %*% temp
      u <- apply(as.matrix(temp0), 2, sum)
      fit$rscore <- coxph.wtest(t(temp0) %*% temp0, u, 
                                control$toler.chol)$test
    }
    if (length(fit$coefficients) && is.null(fit$wald.test)) {
      nabeta <- !is.na(fit$coefficients)
      if (is.null(init)) 
        temp <- fit$coefficients[nabeta]
      else temp <- (fit$coefficients - init[1:length(fit$coefficients)])[nabeta]
      fit$wald.test <- coxph.wtest(fit$var[nabeta, nabeta], 
                                   temp, control$toler.chol)$test
    }
    if (length(cluster)) 
      temp <- concordancefit(Y, fit$linear.predictors, 
                             istrat, weights, cluster = cluster, reverse = TRUE, 
                             timefix = FALSE)
    else temp <- concordancefit(Y, fit$linear.predictors, 
                                istrat, weights, reverse = TRUE, timefix = FALSE)
    if (is.matrix(temp$count)) 
      fit$concordance <- c(colSums(temp$count), concordance = temp$concordance, 
                           std = sqrt(temp$var))
    else fit$concordance <- c(temp$count, concordance = temp$concordance, 
                              std = sqrt(temp$var))
    na.action <- attr(mf, "na.action")
    if (length(na.action)) 
      fit$na.action <- na.action
    if (model) {
      if (length(timetrans)) {
        stop("'model=TRUE' not supported for models with tt terms")
      }
      fit$model <- mf
    }
    if (x) {
      fit$x <- X
      if (length(timetrans)) 
        fit$strata <- istrat
      else if (length(strats)) 
        fit$strata <- strata.keep
    }
    if (y) 
      fit$y <- Y
    fit$timefix <- control$timefix
  }
  if (!is.null(weights) && any(weights != 1)) 
    fit$weights <- weights
  names(fit$means) <- names(fit$coefficients)
  if (multi) {
    fit$transitions <- transitions
    fit$states <- states
    fit$cmap <- cmap
    fit$resid <- rowsum(fit$resid, xstack$rindex)
    names(fit$coefficients) <- seq(along = fit$coefficients)
    if (x) 
      fit$strata <- istrat
    class(fit) <- c("coxphms", class(fit))
  }
  fit$formula <- formula(Terms)
  if (length(xlevels) > 0) 
    fit$xlevels <- xlevels
  fit$contrasts <- contr.save
  if (any(offset != 0)) 
    fit$offset <- offset
  fit$call <- Call
  fit
}



coxpenal.fit<-function (x, y, strata, offset, init, control, weights, method, 
          rownames, pcols, pattr, assign) 
{
  eps <- control$eps
  n <- nrow(y)
  if (is.matrix(x)) 
    nvar <- ncol(x)
  else if (length(x) == 0) 
    stop("Must have an X variable")
  else nvar <- 1
  if (missing(offset) || is.null(offset)) 
    offset <- rep(0, n)
  if (missing(weights) || is.null(weights)) 
    weights <- rep(1, n)
  else {
    if (any(weights <= 0)) 
      stop("Invalid weights, must be >0")
  }
  if (ncol(y) == 3) {
    if (length(strata) == 0) {
      sorted <- cbind(order(-y[, 2], y[, 3]), order(-y[, 
                                                       1])) - 1L
      newstrat <- as.integer(n)
    }
    else {
      sorted <- cbind(order(strata, -y[, 2], y[, 3]), 
                      order(strata, -y[, 1])) - 1L
      newstrat <- as.integer(cumsum(table(strata)))
    }
    status <- y[, 3]
    andersen <- TRUE
  }
  else {
    if (length(strata) == 0) {
      sorted <- order(-y[, 1], y[, 2]) - 1L
      newstrat <- as.integer(n)
    }
    else {
      sorted <- order(strata, -y[, 1], y[, 2]) - 1L
      newstrat <- as.integer(cumsum(table(strata)))
    }
    status <- y[, 2]
    andersen <- FALSE
  }
  n.eff <- sum(y[, ncol(y)])
  npenal <- length(pattr)
  if (npenal == 0 || length(pcols) != npenal) 
    stop("Invalid pcols or pattr arg")
  sparse <- sapply(pattr, function(x) !is.null(x$sparse) && 
                     x$sparse)
  if (sum(sparse) > 1) 
    stop("Only one sparse penalty term allowed")
  pterms <- rep(0, length(assign))
  names(pterms) <- names(assign)
  pindex <- rep(0, npenal)
  for (i in 1:npenal) {
    temp <- unlist(lapply(assign, function(x, y) (length(x) == 
                                                    length(y) && all(x == y)), pcols[[i]]))
    if (sparse[i]) 
      pterms[temp] <- 2
    else pterms[temp] <- 1
    pindex[i] <- (seq(along.with = temp))[temp]
  }
  if ((sum(pterms == 2) != sum(sparse)) || (sum(pterms > 0) != 
                                            npenal)) 
    stop("pcols and assign arguments disagree")
  if (any(pindex != sort(pindex))) {
    temp <- order(pindex)
    pindex <- pindex[temp]
    pcols <- pcols[temp]
    pattr <- pattr[temp]
  }
  ptype <- any(sparse) + 2 * (any(!sparse))
  f.expr1 <- function(coef) NULL
  f.expr2 <- function(coef) NULL
  if (any(sparse)) {
    sparse.attr <- (pattr[sparse])[[1]]
    fcol <- unlist(pcols[sparse])
    if (length(fcol) > 1) 
      stop("Sparse term must be single column")
    xx <- x[, -fcol, drop = FALSE]
    for (i in 1:length(assign)) {
      j <- assign[[i]]
      if (j[1] > fcol) 
        assign[[i]] <- j - 1
    }
    for (i in 1:npenal) {
      j <- pcols[[i]]
      if (j[1] > fcol) 
        pcols[[i]] <- j - 1
    }
    frailx <- x[, fcol]
    frailx <- match(frailx, sort(unique(frailx)))
    nfrail <- max(frailx)
    nvar <- nvar - 1
    pfun1 <- sparse.attr$pfun
    f.expr1 <- function(coef) {
      coxlist1$coef <- coef
      if (is.null(extra1)) 
        temp <- pfun1(coef, theta1, n.eff)
      else temp <- pfun1(coef, theta1, n.eff, extra1)
      if (!is.null(temp$recenter)) 
        coxlist1$coef <- coxlist1$coef - as.double(temp$recenter)
      if (!temp$flag) {
        coxlist1$first <- -as.double(temp$first)
        coxlist1$second <- as.double(temp$second)
      }
      coxlist1$penalty <- -as.double(temp$penalty)
      coxlist1$flag <- as.logical(temp$flag)
      if (any(sapply(coxlist1, length) != c(rep(nfrail, 
                                                3), 1, 1))) 
        stop("Incorrect length in coxlist1")
      coxlist1
    }
    if (!is.null(getOption("survdebug"))) 
      debug(f.expr1)
    coxlist1 <- list(coef = double(nfrail), first = double(nfrail), 
                     second = double(nfrail), penalty = 0, flag = FALSE)
  }
  else {
    xx <- x
    frailx <- 0
    nfrail <- 0
  }
  if (sum(!sparse) > 0) {
    full.imat <- !all(unlist(lapply(pattr, function(x) x$diag)))
    ipenal <- (1:length(pattr))[!sparse]
    f.expr2 <- function(coef) {
      coxlist2$coef <- coef
      pentot <- 0
      for (i in ipenal) {
        pen.col <- pcols[[i]]
        coef <- coxlist2$coef[pen.col]
        if (is.null(extralist[[i]])) 
          temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], 
                                      n.eff)
        else temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], 
                                         n.eff, extralist[[i]])
        if (!is.null(temp$recenter)) 
          coxlist2$coef[pen.col] <- coxlist2$coef[pen.col] - 
          temp$recenter
        if (temp$flag) 
          coxlist2$flag[pen.col] <- TRUE
        else {
          coxlist2$flag[pen.col] <- FALSE
          coxlist2$first[pen.col] <- -temp$first
          if (full.imat) {
            tmat <- matrix(coxlist2$second, nvar, nvar)
            tmat[pen.col, pen.col] <- temp$second
            coxlist2$second <- c(tmat)
          }
          else coxlist2$second[pen.col] <- temp$second
        }
        pentot <- pentot - temp$penalty
      }
      coxlist2$penalty <- as.double(pentot)
      if (any(sapply(coxlist2, length) != length2)) 
        stop("Length error in coxlist2")
      coxlist2
    }
    if (!is.null(getOption("survdebug"))) 
      debug(f.expr2)
    if (full.imat) {
      coxlist2 <- list(coef = double(nvar), first = double(nvar), 
                       second = double(nvar * nvar), penalty = 0, flag = rep(FALSE, 
                                                                             nvar))
      length2 <- c(nvar, nvar, nvar * nvar, 1, nvar)
    }
    else {
      coxlist2 <- list(coef = double(nvar), first = double(nvar), 
                       second = double(nvar), penalty = 0, flag = rep(FALSE, 
                                                                      nvar))
      length2 <- c(nvar, nvar, nvar, 1, nvar)
    }
  }
  else full.imat <- FALSE
  if (nfrail > 0) 
    finit <- rep(0, nfrail)
  else finit <- 0
  if (!missing(init) && !is.null(init)) {
    if (length(init) != nvar) {
      if (length(init) == (nvar + nfrail)) {
        finit <- init[-(1:nvar)]
        init <- init[1:nvar]
      }
      else stop("Wrong length for inital values")
    }
  }
  else init <- double(nvar)
  cfun <- lapply(pattr, function(x) x$cfun)
  parmlist <- lapply(pattr, function(x, eps) c(x$cparm, eps2 = eps), 
                     sqrt(eps))
  extralist <- lapply(pattr, function(x) x$pparm)
  iterlist <- vector("list", length(cfun))
  thetalist <- vector("list", length(cfun))
  printfun <- lapply(pattr, function(x) x$printfun)
  for (i in 1:length(cfun)) {
    temp <- (cfun[[i]])(parmlist[[i]], iter = 0)
    if (sparse[i]) {
      theta1 <- temp$theta
      extra1 <- extralist[[i]]
    }
    thetalist[[i]] <- temp$theta
    iterlist[[i]] <- temp
  }
  temp1 <- c("x", "coef", "plik", "loglik", "status", "neff", 
             "df", "trH")
  temp2 <- c("frailx", "coxfit$fcoef", "loglik1", "coxfit$loglik", 
             "status", "n.eff")
  temp3 <- c("xx[,pen.col]", "coxfit$coef[pen.col]", "loglik1", 
             "coxfit$loglik", "status", "n.eff")
  calls <- vector("expression", length(cfun))
  cargs <- lapply(pattr, function(x) x$cargs)
  for (i in 1:length(cfun)) {
    tempchar <- paste("(cfun[[", i, "]])(parmlist[[", i, 
                      "]], iter,", "iterlist[[", i, "]]")
    temp2b <- c(temp2, paste("pdf[", i, "]"), paste("trH[", 
                                                    i, "]"))
    temp3b <- c(temp3, paste("pdf[", i, "]"), paste("trH[", 
                                                    i, "]"))
    if (length(cargs[[i]]) == 0) 
      calls[i] <- parse(text = paste(tempchar, ")"))
    else {
      temp <- match(cargs[[i]], temp1)
      if (any(is.na(temp))) 
        stop(paste((cargs[[i]])[is.na(temp)], "not matched"))
      if (sparse[i]) 
        temp4 <- paste(temp2b[temp], collapse = ",")
      else temp4 <- paste(temp3b[temp], collapse = ",")
      calls[i] <- parse(text = paste(paste(tempchar, temp4, 
                                           sep = ","), ")"))
    }
  }
  need.df <- any(!is.na(match(c("df", "trH"), unlist(cargs))))
  varnames <- dimnames(xx)[[2]]
  for (i in 1:length(cfun)) {
    if (!is.null(pattr[[i]]$varname)) 
      varnames[pcols[[i]]] <- pattr[[i]]$varname
  }
  rho <- environment()
  storage.mode(y) <- storage.mode(weights) <- "double"
  storage.mode(xx) <- storage.mode(offset) <- "double"
  if (andersen) 
    coxfit <- .C(Cagfit5a, as.integer(n), as.integer(nvar), 
                 y, xx, offset, weights, newstrat, sorted, means = double(nvar), 
                 coef = as.double(init), u = double(nvar), loglik = double(1), 
                 as.integer(method == "efron"), as.integer(ptype), 
                 as.integer(full.imat), as.integer(nfrail), as.integer(frailx), 
                 f.expr1, f.expr2, rho)
  else coxfit <- .C(Ccoxfit5a, as.integer(n), as.integer(nvar), 
                    y, xx, offset, weights, newstrat, sorted, means = double(nvar), 
                    coef = as.double(init), u = double(nvar), loglik = double(1), 
                    as.integer(method == "efron"), as.integer(ptype), as.integer(full.imat), 
                    as.integer(nfrail), as.integer(frailx), f.expr1, f.expr2, 
                    rho)
  loglik0 <- coxfit$loglik
  means <- coxfit$means
  iter2 <- 0
  iterfail <- NULL
  thetasave <- unlist(thetalist)
  for (outer in 1:control$outer.max) {
    if (andersen) 
      coxfit <- .C(Cagfit5b, iter = as.integer(control$iter.max), 
                   as.integer(n), as.integer(nvar), as.integer(newstrat), 
                   coef = as.double(init), u = double(nvar + nfrail), 
                   hmat = double(nvar * (nvar + nfrail)), hinv = double(nvar * 
                                                                          (nvar + nfrail)), loglik = double(1), flag = integer(1), 
                   as.double(control$eps), as.double(control$toler.chol), 
                   as.integer(method == "efron"), as.integer(nfrail), 
                   fcoef = as.double(finit), fdiag = double(nfrail + 
                                                              nvar), f.expr1, f.expr2, rho)
    else coxfit <- .C(Ccoxfit5b, iter = as.integer(control$iter.max), 
                      as.integer(n), as.integer(nvar), as.integer(newstrat), 
                      coef = as.double(init), u = double(nvar + nfrail), 
                      hmat = double(nvar * (nvar + nfrail)), hinv = double(nvar * 
                                                                             (nvar + nfrail)), loglik = double(1), flag = integer(1), 
                      as.double(control$eps), as.double(control$toler.chol), 
                      as.integer(method == "efron"), as.integer(nfrail), 
                      fcoef = as.double(finit), fdiag = double(nfrail + 
                                                                 nvar), f.expr1, f.expr2, rho)
    iter <- outer
    iter2 <- iter2 + coxfit$iter
    if (coxfit$iter >= control$iter.max) 
      iterfail <- c(iterfail, iter)
    temp <- rep(FALSE, nvar + nfrail)
    if (nfrail > 0) 
      temp[1:nfrail] <- coxlist1$flag
    if (ptype > 1) 
      temp[nfrail + 1:nvar] <- coxlist2$flag
    fdiag <- ifelse(temp, 0, coxfit$fdiag)
    if (need.df) {
      if (nfrail > 0) 
        temp1 <- coxlist1$second
      else temp1 <- 0
      if (ptype > 1) 
        temp2 <- coxlist2$second
      else temp2 <- 0
      dftemp <- coxpenal.df(matrix(coxfit$hmat, ncol = nvar), 
                            matrix(coxfit$hinv, ncol = nvar), fdiag, assign, 
                            ptype, nvar, temp1, temp2, pindex[sparse])
      df <- dftemp$df
      var <- dftemp$var
      var2 <- dftemp$var2
      pdf <- df[pterms > 0]
      trH <- dftemp$trH[pterms > 0]
    }
    if (nfrail > 0) 
      penalty <- -coxlist1$penalty
    else penalty <- 0
    if (ptype > 1) 
      penalty <- penalty - coxlist2$penalty
    loglik1 <- coxfit$loglik + penalty
    if (iter == 1) 
      penalty0 <- penalty
    done <- TRUE
    for (i in 1:length(cfun)) {
      pen.col <- pcols[[i]]
      temp <- eval(calls[i])
      if (sparse[i]) 
        theta1 <- temp$theta
      thetalist[[i]] <- temp$theta
      iterlist[[i]] <- temp
      done <- done & temp$done
    }
    if (done) 
      break
    if (iter == 1) {
      init <- coefsave <- coxfit$coef
      finit <- fsave <- coxfit$fcoef
      thetasave <- cbind(thetasave, unlist(thetalist))
    }
    else {
      temp <- as.vector(unlist(thetalist))
      coefsave <- cbind(coefsave, coxfit$coef)
      fsave <- cbind(fsave, coxfit$fcoef)
      howclose <- apply((thetasave - temp)^2, 2, sum)
      which <- min((1:iter)[howclose == min(howclose)])
      if (nvar > 0) 
        init <- coefsave[, which]
      if (nfrail > 0) 
        finit <- fsave[, which]
      thetasave <- cbind(thetasave, temp)
    }
  }
  if (nfrail > 0) {
    lp <- offset + coxfit$fcoef[frailx]
    if (nvar > 0) 
      lp <- lp + x[, -fcol, drop = FALSE] %*% coxfit$coef - 
        sum(means * coxfit$coef)
  }
  else lp <- offset + as.vector(x %*% coxfit$coef) - sum(means * 
                                                           coxfit$coef)
  if (andersen) {
    .C(Cagfit5c, as.integer(nvar))
    if (length(strata) < nrow(y)) 
      strata <- rep(1L, nrow(y))
    resid <- NA # CHANGE!
  }
  else {
    expect <- .C(Ccoxfit5c, as.integer(n), as.integer(nvar), 
                 as.integer(newstrat), as.integer(method == "efron"), 
                 expect = double(n))$expect
    resid <- status - expect
  }
  names(resid) <- NULL # CHANGE!
  if (!need.df) {
    if (nfrail > 0) 
      temp1 <- coxlist1$second
    else temp1 <- 0
    if (ptype > 1) 
      temp2 <- coxlist2$second
    else temp2 <- 0
    dftemp <- coxpenal.df(matrix(coxfit$hmat, ncol = nvar), 
                          matrix(coxfit$hinv, ncol = nvar), fdiag, assign, 
                          ptype, nvar, temp1, temp2, pindex[sparse])
    df <- dftemp$df
    trH <- dftemp$trH
    var <- dftemp$var
    var2 <- dftemp$var2
  }
  if (control$iter.max > 1 && length(iterfail) > 0) 
    warning(paste("Inner loop failed to coverge for iterations", 
                  paste(iterfail, collapse = " ")))
  which.sing <- (fdiag[nfrail + 1:nvar] == 0)
  coef <- coxfit$coef
  names(coef) <- varnames
  coef[which.sing] <- NA
  names(iterlist) <- names(pterms[pterms > 0])
  if (nfrail > 0) {
    if (nvar > 0) {
      list(coefficients = coef, var = var, var2 = var2, 
           loglik = c(loglik0, loglik1), iter = c(iter, 
                                                  iter2), linear.predictors = as.vector(lp), 
           residuals = resid, means = means, method = method, 
           class = c("coxph.penal", "coxph"), frail = coxfit$fcoef, 
           fvar = dftemp$fvar, df = df, df2 = dftemp$df2, 
           penalty = c(penalty0, penalty), pterms = pterms, 
           assign2 = assign, history = iterlist, coxlist1 = coxlist1, 
           printfun = printfun)
    }
    else {
      list(loglik = c(loglik0, loglik1), iter = c(iter, 
                                                  iter2), linear.predictors = as.vector(lp), residuals = resid, 
           means = means, method = method, class = c("coxph.penal", 
                                                     "coxph"), frail = coxfit$fcoef, fvar = dftemp$fvar, 
           df = df, df2 = dftemp$df2, penalty = c(penalty0, 
                                                  penalty), pterms = pterms, assign2 = assign, 
           history = iterlist, printfun = printfun)
    }
  }
  else {
    list(coefficients = coef, var = var, var2 = var2, loglik = c(loglik0, 
                                                                 loglik1), iter = c(iter, iter2), linear.predictors = lp, 
         residuals = resid, means = means, method = method, 
         class = c("coxph.penal", "coxph"), df = df, df2 = dftemp$df2, 
         penalty = c(penalty0, penalty), pterms = pterms, 
         assign2 = assign, history = iterlist, coxlist2 = coxlist2, 
         printfun = printfun)
  }
}

terms.inner<-function (x) 
{
  if (inherits(x, "formula")) {
    if (length(x) == 3) 
      c(terms.inner(x[[2]]), terms.inner(x[[3]]))
    else terms.inner(x[[2]])
  }
  else if (inherits(x, "call") && (x[[1]] != as.name("$") && 
                                   x[[1]] != as.name("["))) {
    if (x[[1]] == "+" || x[[1]] == "*" || x[[1]] == "-") {
      if (length(x) == 3) 
        c(terms.inner(x[[2]]), terms.inner(x[[3]]))
      else terms.inner(x[[2]])
    }
    else if (x[[1]] == as.name("Surv")) 
      unlist(lapply(x[-1], terms.inner))
    else terms.inner(x[[2]])
  }
  else (deparse(x))
}

coxpenal.df<-function (hmat, hinv, fdiag, assign.list, ptype, nvar, pen1, 
                       pen2, sparse) 
{
  if (ptype == 1 & nvar == 0) {
    hdiag <- 1/fdiag
    list(fvar2 = (hdiag - pen1) * fdiag^2, df = sum((hdiag - 
                                                       pen1) * fdiag), fvar = fdiag, trH = sum(fdiag))
  }
  else if (ptype == 2) {
    hmat.full <- t(hmat) %*% (ifelse(fdiag == 0, 0, 1/fdiag) * 
                                hmat)
    hinv.full <- hinv %*% (fdiag * t(hinv))
    if (length(pen2) == length(hmat.full)) 
      imat <- hmat.full - pen2
    else imat <- hmat.full - diag(pen2)
    var <- hinv.full %*% imat %*% hinv.full
    if (length(assign.list) == 1) 
      list(var2 = var, df = sum(imat * hinv.full), trH = sum(diag(hinv.full)), 
           var = hinv.full)
    else {
      df <- trH <- NULL
      d2 <- diag(hinv.full)
      for (i in assign.list) {
        temp <- coxph.wtest(hinv.full[i, i], var[i, 
                                                 i])$solve
        if (is.matrix(temp)) 
          df <- c(df, sum(diag(temp)))
        else df <- c(df, sum(temp))
        trH <- c(trH, sum(d2[i]))
      }
      list(var2 = var, df = df, trH = trH, var = hinv.full)
    }
  }
  else {
    nf <- length(fdiag) - nvar
    nr1 <- 1:nf
    nr2 <- (nf + 1):(nf + nvar)
    d1 <- fdiag[nr1]
    d2 <- fdiag[nr2]
    temp <- t(hinv[nr1, ])
    temp2 <- t(hinv[nr2, , drop = FALSE])
    A.diag <- d1 + c(rep(1, nvar) %*% (temp^2 * d2))
    B <- hinv[nr1, ] %*% (d2 * temp2)
    C <- hinv[nr2, ] %*% (d2 * temp2)
    var2 <- C - t(B) %*% (pen1 * B)
    if (ptype == 3) {
      hmat.22 <- t(hmat) %*% (ifelse(fdiag == 0, 0, 1/fdiag) * 
                                hmat)
      temp <- C - coxph.wtest(hmat.22, diag(nvar))$solve
      if (nvar == 1) {
        var2 <- var2 - C * pen2 * C
        temp2 <- c(temp * pen2)
      }
      else if (length(pen2) == nvar) {
        var2 <- var2 - C %*% (pen2 * C)
        temp2 <- sum(diag(temp) * pen2)
      }
      else {
        var2 <- var2 - C %*% matrix(pen2, nvar) %*% 
          C
        temp2 <- sum(diag(temp * pen2))
      }
    }
    else temp2 <- 0
    df <- trH <- NULL
    cdiag <- diag(C)
    for (i in 1:length(assign.list)) {
      if (sparse == i) {
        df <- c(df, nf - (sum(A.diag * pen1) + temp2))
        trH <- c(trH, sum(A.diag))
      }
      else {
        j <- assign.list[[i]]
        temp <- coxph.wtest(C[j, j], var2[j, j])$solve
        if (is.matrix(temp)) 
          df <- c(df, sum(diag(temp)))
        else df <- c(df, sum(temp))
        trH <- c(trH, sum(cdiag[j]))
      }
    }
    list(var = C, df = df, trH = trH, fvar = A.diag, var2 = var2)
  }
}


