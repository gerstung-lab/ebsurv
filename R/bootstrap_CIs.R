#' Bootstrap confidence intervals for regression coefficients
#' 
#' This function computes non-parametric bootstrap confidence intervals for the regression coefficients estimated by CoxRFX.
#' @param mstate_data_expanded Data in `long format`, possibly with `expanded` covariates
#' (as obtained by running mstate::expand.covs).
#' @param which_group A character vector with the same meaning as the `groups` argument of `CoxRFX` but named (with the covariate names).
#' @param min_nr_samples The confidence interval of any coefficient is based on a number of bootstrap samples which at least as high as this argument. See details.
#' @param output Determines the sort of output. See value.
#' @return For each regression coefficient, the confidence intervals and the number of bootstrap samples on which they are based, if the `output` argument is equal to `CIs`;
#' if `output` is equal to `CIs_and_coxrfx_fits`, also the CoxRFX objects for each bootstrap sample.  
#' 
#' @author rc28
#' @export

boot_coxrfx<-function(mstate_data_expanded,which_group,min_nr_samples=100,output="CIs"){
  coxrfx_fits_boot<-vector("list")
  rownames(mstate.data.expanded)<-1:nrow(mstate.data.expanded)
  boot_matrix<-matrix(nrow=0,ncol = ncol(mstate.data.expanded)-8,dimnames = list(NULL,names(mstate.data.expanded)[-(1:8)]))
  j<-1
  repeat{
    boot_samples_trans_1<-sample(rownames(mstate.data.expanded[mstate.data.expanded$trans==1,]),replace = T)
    boot_samples_trans_2<-sample(rownames(mstate.data.expanded[mstate.data.expanded$trans==2,]),replace = T)
    boot_samples_trans_3<-sample(rownames(mstate.data.expanded[mstate.data.expanded$trans==3,]),replace = T)
    boot_samples<-c(boot_samples_trans_1,boot_samples_trans_2,boot_samples_trans_3) 
    
    mstate.data.expanded.boot<-mstate.data.expanded[boot_samples,]
    covariate_df<-mstate.data.expanded.boot[-(1:8)]
    covariate_df<-covariate_df[,which(apply(covariate_df, 2, var)>0)]
    groups2<-which_group[names(covariate_df)]
    covariate_df$transition<-mstate.data.expanded.boot$trans
    
    coxrfx_fits_boot[[j]]<-CoxRFX(covariate_df,Surv(mstate.data.expanded.boot$time,mstate.data.expanded.boot$status),groups =groups2,max.iter = 200 )
    
    if(coxrfx_fits_boot[[j]]$iter[1]!=as.list(coxrfx_fits_boot[[j]]$call)$max.iter & sum(is.na(coxrfx_fits_boot[[j]]$coefficients))==0){
      boot_matrix<-rbind(boot_matrix,rep(NA,ncol(boot_matrix)))
      boot_matrix[j,names(coxrfx_fits_boot[[j]]$coefficients)]<-coxrfx_fits_boot[[j]]$coefficients
      print(min(apply(boot_matrix, 2, function(x) sum(!is.na(x)))))
      j<-j+1
    } 
    if(min(apply(boot_matrix, 2, function(x) sum(!is.na(x))))==min_nr_samples) break
    
  }
  
  CIs<-apply(boot_matrix,2,HDInterval::hdi,credMass=0.95)
  CIs<-rbind(CIs,apply(boot_matrix, 2, function(x) sum(!is.na(x))))
  dimnames(CIs)[[1]][3]<-"n_samples"
  if(output=="CIs_and_coxrfx_fits"){
    return(list(CIs=CIs,coxrfx_fits_boot=coxrfx_fits_boot))
  }else if(output=="CIs"){
    return(CIs)
  }
}