
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
  stratum<-Z$stratum
  Z$stratum<-NULL
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
								sumTerm,"strata(stratum)"), 
						collapse=" + ")))
		fit <- coxph(formula, ...)
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
		formula <- as.formula(paste(as.character(call["surv"]),"~",paste(colnames(Z)[j], collapse="+"),"+strata(stratum)"))
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
	fit$stratum<-stratum
	class(fit) <- c("coxrfx", class(fit))
	return(fit)
}


#' A summary method for CoxRFX models
#' 
#' This model prints the means and variances for each groups of covariates, as well as the variance components.
#' For the means a Wald test with 1 df is computed testing the null-hypothesis of being zero.
#' 
#' The null-hypothesis of zero variance is tested using a combined Wald test that all coefficients in the group are
#' identical to the mean. Gray (1992) suggests to use \eqn{\beta H^{-1} \beta} as a test statistic in a chi-square
#' test with \eqn{\mathrm{tr}[H^{-1} I]}{tr[H^{-1} I]} df, where H is the Hessian of the penalised model
#' and I is the Hessian of the unpenalised coxph model. Note that all variables taken over the subset of interest only. 
#' As noted by Therneau (2003) this test may be somewhat optimistic. Here, we are using \eqn{z^2 = \sum_i\beta_i^2/H_{ii}},
#' which appears to be a more conservative choice, but the consequences remain to be thoroughly evaluated.  
#' 
#' @note Note that there is a lively debate about testing random effects in generalised linear models. See for example
#' http://glmm.wikidot.com/faq
#' 
#' @references  R. J. Gray (1992). Flexible Methods for Analyzing Survival Data Using Splines, with Applications to Breast Cancer Prognosis. Journal of the American Statistical Association, 87:942-951. http://dx.doi.org/10.1080/01621459.1992.10476248
#' T. M. Therneau, P. M. Grambsch, and V. S. Pankratz (2003). Penalized Survival Models and Frailty. Journal of Computational and Graphical Statistics, 12:156-175. http://dx.doi.org/10.1198/1061860031365
#' @param object A CoxRFX model
#' @param ... Currently unused
#' @return NULL
#' 
#' @author Moritz Gerstung
#' @export
#' @method summary coxrfx
summary.coxrfx <- function(object, ...){
	which.mu <- names(object$mu)[object$mu!=0]
	p <- z <- s <- object$mu
	z[which.mu] <- object$mu[which.mu]/sqrt(diag(as.matrix(object$mu.var2)))
	s[which.mu] <- sqrt(diag(as.matrix(object$mu.var2)))
	p <- pchisq(z^2,1,lower.tail=FALSE)
	p[!names(p) %in% which.mu] <- NA
	cat("Means:\n")
	show(format(data.frame(mean=object$mu, sd=s, z=z, p.val=p, sig=sig2star(p)), digits=2))
	cat("\nVariances - p-values only indicative:\n")
	v <- object$sigma2
	c <- coef(object) - object$mu[object$groups] ## centred coefficients
	chisq <- sapply(split(c^2/diag(object$Hinv)[1:length(c)], object$groups), sum)
	df <- object$df[-(nlevels(object$groups)+1)]
	p <- pchisq(chisq, df, lower.tail=FALSE)
#	f <- as.numeric(table(x$groups)/x$df[-(nlevels(x$groups)+1)])
#	u <- sapply(split(coef(x), x$groups), function(x) sum((x-mean(x))^2)/qchisq(0.025, length(x)))
#	l <- sapply(split(coef(x), x$groups), function(x) sum((x-mean(x))^2)/qchisq(0.975, length(x)))
	show(format(data.frame(sigma2=v, chisq=chisq, df = df, p.val=p, sig=sig2star(p)), digits=2))
	cat("\nPartial log hazard:\n")
	newZ <- object$Z[setdiff(1:nrow(object$Z), object$na.action),]
	p <- PartialRisk(object, newZ = newZ)
	v <- VarianceComponents(object, newZ = newZ)
	e <- colMeans(PartialRiskVar(object, newZ = newZ))
	show(format(data.frame(`Cov[g,g]`=c(diag(cov(p)), TOTAL=NaN), `Sum(Cov[,g])`=c(rowSums(cov(p)),TOTAL=sum(cov(p))), `MSE`=c(e, TOTAL=v[length(v)]), check.names = FALSE),  digits=2))
}

#' Print method for CoxRFX
#' 
#' This function implicitly calls summary.CoxRFX().
#' @param x CoxRFX
#' @param ... Currently unused
#' @return NULL
#' 
#' @author Moritz Gerstung
#' @method print CoxRFX
#' @export
print.CoxRFX <- function(x, ...){
	summary.CoxRFX(x)
}

#' Convert factor to integer.
#' @param F A factor vector.
#' @return A data.frame with columns corresponding to levels in the factor. 
#' @details An internal function of \code{CoxRFX}, not meant
#' to called directly by the user.
#' @author Moritz Gerstung
#' @seealso \code{\link{CoxRFX}}
#' @export
MakeInteger <- function(F){
  res <- as.data.frame(lapply(levels(F), `==`, F))
  colnames(res) <- levels(F)
  res + 0
}

