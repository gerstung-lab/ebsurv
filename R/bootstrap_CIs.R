
#' Bootstrap confidence intervals for transition probabilities
#' 
#' Generates 95\% highest density bootstrap interval estimates for transition probabilities computed using \code{probtrans_ebsurv} (semi-Markov version).
#' 
#' @param coxrfx_fits_boot The list of CoxRFX objects obtained by running \code{boot_coxrfx}.
#' @param patient_data (Single) patient data in `long format`, possibly with `expanded` covariates
#' (as obtained by running \code{mstate::expand.covs}).
#' @param tmat Transition matrix for the multi-state model, as obtained by running
#' \code{mstate::transMat}.
#' @return Interval estimates for transition probabilities. 
#' @author Rui Costa
#' @seealso \code{\link{probtrans_ebsurv}}; \code{\link{boot_coxrfx}}; 
#' \code{\link[mstate]{transMat}}; \code{\link[mstate]{expand.covs}}
#' @export

boot_probtrans<-function(coxrfx_fits_boot,patient_data,tmat){
  msfit_objects_boot<-vector("list",length(coxrfx_fits_boot))
  probtrans_objects_boot<-vector("list",length(coxrfx_fits_boot))
  for(i in 1:length(coxrfx_fits_boot)){
    print(i)
    covariate_df<-as.data.frame(coxrfx_fits_boot[[i]]$Z)
    covariate_df$transition<-coxrfx_fits_boot[[i]]$transition
    mstate_data_expanded.boot<-list()
    mstate_data_expanded.boot$time<-coxrfx_fits_boot[[i]]$surv[,1]
    mstate_data_expanded.boot$status<-coxrfx_fits_boot[[i]]$surv[,2]
    patient_data2<-patient_data[names(patient_data)%in%names(covariate_df)]
    patient_data2$strata<-patient_data$strata
    
    environment(coxrfx_fits_boot[[1]]$formula)$covariate_df<-covariate_df
    
    msfit_objects_boot[[i]]<-msfit_generic(coxrfx_fits_boot[[i]],patient_data2,trans=tmat)
    probtrans_objects_boot[[i]]<-probtrans_ebsurv(msfit_objects_boot[[i]],"semiMarkov")[[1]]
    probtrans_objects_boot[[i]]<-probtrans_objects_boot[[i]][sapply(seq(from=0,to=max(probtrans_objects_boot[[i]]$time),length.out = 400),function(x) which.min(abs(probtrans_objects_boot[[i]]$time-x))),]
    
  }

  probtrans_CIs<-lapply(colnames(tmat),CIs_for_target_state)
  names(probtrans_CIs)<-colnames(tmat)
  return(probtrans_CIs)
}

#' Ancillary function to \code{boot_probtrans}.
#' 
#' Extracts the bootstrap estimates of transition probabilities for
#' target state `tstate` from a list
#' with bootstrap estimates of transition probabilities into multiple states.
#' This function is not meant to be called by the user.
#' 
#' @param list_object A list in which each individual element is a single
#' bootstrap estimate of the probability of transition
#' into different states.
#' @param tstate The state whose bootstrap estimates of transition probabilities we wish to extract
#' from \code{list_object}. 
#' @return Bootstrap estimates of transition probabilities into target state `tstate`. 
#' @details This function is an ancillary function of \code{CIs_for_target_state}, which
#' in turn is an ancillary function of \code{boot_probtrans}.
#' @author Rui Costa
#' @seealso \code{\link{CIs_for_target_state}}; \code{\link{boot_probtrans}} 
#' @export

extract_function<-function(list_object,tstate){
  as.vector(list_object[tstate])
}

#' Ancillary function of \code{boot_probtrans}.
#' 
#' Computes 95\% highest density bootstrap confidence 
#' intervals for the transition probabilities into \code{target_state}, 
#' given a list object with boostrap estimates of transition probabilities into multiple states. This 
#' function is not meant to be called by the user.
#' 
#' @param target_state The target state for whose transition probabilties the confidence intervals
#' are computed.
#' @return 95\% highest density bootstrap confidence intervals for the transition
#' probabilities into \code{target_state}. 
#' @details Uses function \code{extract_function}.
#' @author Rui Costa
#' @seealso \code{\link{boot_probtrans}}; \code{\link{extract_function}}.
#' @export

CIs_for_target_state<-function(target_state){
  target_state_boot_samples<-as.data.frame(sapply(probtrans_objects_boot, extract_function,tstate=target_state))
  apply(target_state_boot_samples,1,HDInterval::hdi,credMass=0.95)
}

#' Bootstrap confidence intervals for regression coefficients
#' 
#' This function computes 95\% highest density bootstrap confidence intervals (non-parametric) for the regression coefficients estimated by CoxRFX.
#' 
#' @param mstate_data_expanded Data in `long format`, possibly with `expanded` covariates (as obtained by running mstate::expand.covs).
#' @param which_group A character vector with the same meaning as the `groups` argument of the function \code{CoxRFX} but named (with the covariate names).
#' @param min_nr_samples The confidence interval of any coefficient is based on a number of bootstrap samples at least as high as this argument. See details.
#' @param output Determines the sort of output. See value.
#' @param ... Further arguments to the CoxRFX function.
#' @return For each regression coefficient, the confidence intervals and the number of bootstrap samples on which they are based, if the `output` argument is equal to `CIs`; if `output` is equal to `CIs_and_coxrfx_fits`, also the \code{CoxRFX} objects for each bootstrap sample.  
#' @details In a given bootstrap sample there might not be enough information to generate 
#' estimates for all coefficients. If a covariate has little or no variation in a given bootstrap sample, 
#' no estimate of its coefficient will be computed. The present function will
#' keep taking bootstrap samples until every coefficient has been estimated
#' at least \code{min_nr_samples} times.
#' @author Rui Costa
#' @export

boot_coxrfx<-function(mstate_data_expanded,which_group,min_nr_samples=100,output="CIs",...){
  coxrfx_fits_boot<-vector("list")
  rownames(mstate_data_expanded)<-1:nrow(mstate_data_expanded)
  boot_matrix<-matrix(nrow=0,ncol = ncol(mstate_data_expanded)-8,dimnames = list(NULL,names(mstate_data_expanded)[-(1:8)]))
  j<-1
  repeat{
    boot_samples_trans_1<-sample(rownames(mstate_data_expanded[mstate_data_expanded$trans==1,]),replace = T)
    boot_samples_trans_2<-sample(rownames(mstate_data_expanded[mstate_data_expanded$trans==2,]),replace = T)
    boot_samples_trans_3<-sample(rownames(mstate_data_expanded[mstate_data_expanded$trans==3,]),replace = T)
    boot_samples<-c(boot_samples_trans_1,boot_samples_trans_2,boot_samples_trans_3) 
    
    mstate_data_expanded.boot<-mstate_data_expanded[boot_samples,]
    
    ##exclude covariates without variance and binary covariates with less than 0.05 cases
    # vars_to_exclude<-vector("list",3)
    # for(i in 1:3){
    #   string_split<-strsplit(names(mstate_data_expanded.boot),"[.]")
    #   var_indices<-sapply(string_split,function(x) x[length(x)]==as.character(i))
    #   dummy_dataset<-mstate_data_expanded.boot[mstate_data_expanded.boot$trans==i,var_indices]
    #   which_have_variance<-apply(dummy_dataset, 2, function(x) var(x)>0)
    #   vars_to_exclude[[i]]<-names(dummy_dataset)[!which_have_variance]
    #   dummy_dataset<-dummy_dataset[which_have_variance]
    #   non_categorical_vars<-paste0(c("age_log","hb","anc_log","plt_log","bm_blasts_logit","ring_sideroblasts_logit","ipss","date"),paste0(".",as.character(i)))
    #   percentage_of_ones<-apply(dummy_dataset[!names(dummy_dataset)%in%non_categorical_vars], 2, function(x) sum(x)/length(x))
    #   which_less_than_five_percent<-which(percentage_of_ones<0.05)
    #   vars_to_exclude[[i]]<-c(vars_to_exclude[[i]],names(percentage_of_ones)[which_less_than_five_percent])
    # }
    # mstate_data_expanded.boot<-mstate_data_expanded.boot[!names(mstate_data_expanded.boot)%in%unlist(vars_to_exclude)]
    # 
    
    covariate_df<-mstate_data_expanded.boot[-(1:8)]
    groups2<-which_group[names(covariate_df)]
    covariate_df$stratum<-mstate_data_expanded.boot$trans
    
    coxrfx_fits_boot[[j]]<-CoxRFX(covariate_df,Surv(mstate_data_expanded.boot$time,mstate_data_expanded.boot$status),groups =groups2,... )
    
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

