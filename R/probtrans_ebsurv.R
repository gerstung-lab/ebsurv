
#'Spline approximations of the cumulative hazard functions
#'
#'Creates a spline approximation for the vector of cumulative hazards of each transition.
#'
#'This function is used by the function \code{probtrans_by_convolution}. It is not meant to be called by the user.
#'
#'@param cumhaz An object of class \code{msfit}, created by 
#'\code{\link{msfit_generic}} or \code{\link{msfit}}.
#'@return A list of estimated cumulative hazard functions (one for each transition).
#'@seealso \code{\link{msfit_generic}}; \code{\link{msfit}}; \code{\link{probtrans_by_convolution}}.
#'@author Rui Costa
#'@export

cumhaz_splines<-function(cumhaz){
  spline_list<-vector("list",length(unique(cumhaz$Haz$trans)))
  for(i in unique(cumhaz$Haz$trans)){
    cumhaz_subset<-cumhaz$Haz[cumhaz$Haz$trans==i,]
    spline_list[[i]]<-splinefun(cumhaz_subset[,"time"], cumhaz_subset[,"Haz"], method="monoH.FC")
  }
  spline_list
}

#'Find all possible paths until absorption from a given starting state
#'
#'\code{unique_paths} finds all possible sequences of states until absorption
#'when the process has a tree-like structure.
#'
#'This function is used by the function \code{\link{probtrans_by_convolution}}. 
#'It is not meant to be called by the user.
#'
#'@param from_state Initial state.
#'@param tmat A transition matrix describing the states and transitions in 
#' the multi-state model, as can be obtained by running
#' \code{\link{transMat}}. 
#' See argument \code{trans} in \code{\link{msprep}} (\code{mstate}
#' package) for more detailed information.
#'@return A matrix where each column is a sequence of states taken by the process until absorption. 
#'There are as many columns as the number of possible paths until absorption.
#'
#'@author Rui Costa
#'@seealso \code{\link{probtrans_by_convolution}};
#' \code{\link{transMat}}.
#'@export

unique_paths<-function(from_state,tmat){
  M<-matrix(from_state)
  repeat{
    m<-NULL
    for(i in 1:ncol(M)){
      last_state_of_col_i<-M[dim(M)[1],i]
      if(is.na(last_state_of_col_i)){
        m<-cbind(m,c(M[,i],NA))
        next
      }
      target_states_of_last_state_of_col_i<-names(which(!is.na(tmat[last_state_of_col_i,])))
      if(length(target_states_of_last_state_of_col_i)>0){
        for(j in target_states_of_last_state_of_col_i){
          m<-cbind(m,c(M[,i],j))
        }
      }else{
        m<-cbind(m,c(M[,i],NA))
      }
      
    }
    M<-m
    if(sum(is.na(M[dim(M)[1],]))==dim(M)[2]){
      return(M)
    } 
  }
}

#'Find the unique sequence of states between two states
#'
#'\code{successful_transitions} finds the unique path between 
#'two states (if there is one) when the process has a tree-like structure.
#'
#'This function is used by 'probtrans_by_convolution_Markov' and 'probtrans_by_convolution_semiMarkov'.
#'It is not meant to be called by the user.
#'
#'@param unique_paths_object
#'@param to_state
#'@return A vector with the unique sequence of states between two states. 
#'
#'@author rc28
#'@export

successful_transitions<-function(unique_paths_object,to_state){
  row_of_paths_object_with_to_state<-which(apply(unique_paths_object,1,function(x) sum(na.omit(x==to_state))>0))
  reduced_unique_paths_object<-matrix(unique(unique_paths_object[1:row_of_paths_object_with_to_state,],MARGIN = 2),ncol = ncol(unique_paths_object))
  sucessful_path_column<-which(apply(reduced_unique_paths_object,2,function(x) sum(x==to_state)>0))
  successful_path<-reduced_unique_paths_object[,sucessful_path_column]
  if(length(successful_path)==1){
    return(NULL)
  }else{
    transitions<-NULL
    for(i in 1:(length(successful_path)-1)){
      transitions<-c(transitions,tmat[successful_path[i],successful_path[i+1]])
    }
    return(transitions)
  }
}

#'Compute the cumulative hazard of leaving a given state 
#'
#'\code{joint_cum_hazard_function} returns the cumulative
#'hazard of leaving state \code{i} to any state that can be
#'reached directly from \code{i}, at each of the time points in \code{t}.
#'
#'This function is not meant to be called by the user.
#'
#'@param t A vector of time points.
#'@param competing_transitions The transitions that can occur when the process
#'is in state \code{i}.
#'@param spline_list A list whose elements are spline functions 
#'approximating the cumulative hazard of making each possible transition from 
#'state \code{i}.
#'@return A vector with the cumulative hazard of leaving a given state evaluated at given time points.
#'
#'@author rc28
#'@export

joint_cum_hazard_function<-function(t,competing_transitions,spline_list){
  if(length(competing_transitions)>0){
    return(sum(sapply(competing_transitions,function(x) spline_list[[x]](t))))
  }else{
    return(rep(0,length(t)))
  }
}

#'Compute subject-specific transition probabilities
#'using convolution.
#'
#'@param cumhaz An 'msfit' object
#'@param model Either 'Markov' or 'semiMarkov'
#'
#'@return An object of class 'probtrans'
#'
#'@export
probtrans_ebsurv<-function(initial_state,cumhaz,model){
  probtrans_object<-lapply(initial_state, probtrans_by_convolution,tmat=cumhaz$trans,cumhaz=cumhaz,model=model)
  probtrans_object$trans<-cumhaz$trans
  probtrans_object$direction<-"forward"
  probtrans_object$predt<-0
  class(probtrans_object)<-"probtrans"
  return(probtrans_object) 
}

#'Compute all transition probabilities from a given state.
#'
#'@param tmat Transition matrix.
#'@param cumhaz 'msfit' object.
#'@param from_state Initial state.
#'@param model
#'
#'@export

probtrans_by_convolution<-function(tmat,cumhaz,from_state,model){
  spline_list<-cumhaz_splines(cumhaz)
  unique_paths_object<-unique_paths(from_state,tmat)
  all_states<-na.omit(unique(as.vector(unique_paths_object)))
  #all_target_states<-all_states[-which(all_states==from_state)]
  maximum_time<-max(cumhaz$Haz$time)
  time<-seq(0,maximum_time,length.out=10000)
  #time<-sort(c(0,sort(coxph.detail(fit)$time)-1,sort(coxph.detail(fit)$time)+1,maximum_time))
  if(model=="semiMarkov"){
    transprobs_for_all_states<-sapply(all_states, probtrans_by_convolution_semiMarkov,cumhaz=cumhaz,tmat=tmat,from_state=from_state,spline_list=spline_list,unique_paths_object=unique_paths_object,time=time)
  }else{
    transprobs_for_all_states<-sapply(all_states, probtrans_by_convolution_Markov,cumhaz=cumhaz,tmat=tmat,from_state=from_state,spline_list=spline_list,unique_paths_object=unique_paths_object,time=time)
  }
  non_reacheable_states<-rownames(tmat)[!rownames(tmat)%in%all_states]
  transprobs_for_all_states<-cbind(transprobs_for_all_states,matrix(0,nrow = length(time),ncol = length(non_reacheable_states),dimnames = list(NULL,non_reacheable_states)))
  
  #order_of_states_according_to_tmat<-rownames(tmat)[rownames(tmat)%in%colnames(transprobs_for_all_states)]
  transprobs_for_all_states<-transprobs_for_all_states[,rownames(tmat)]
  
  as.data.frame(cbind(time,transprobs_for_all_states))
  
}

#'Compute transition probabilities for a given starting state and target state.
#'
#'Compute transition probabilities for a given starting state and target state
#'under a inhomogeneous Markov model
#'
#'@param tmat Transition matrix.
#'@param cumhaz 'msfit' object.
#'@param from_state Initial state.
#'@param to_state Target state.
#'@param spline_list
#'@param unique_paths_object
#'@param time
#'
#'@export
probtrans_by_convolution_Markov<-function(tmat,cumhaz,from_state,to_state,spline_list,unique_paths_object,time){
  row_of_tmat_of_current_state<-which(rownames(tmat)==from_state)
  competing_transitions<-na.omit(tmat[row_of_tmat_of_current_state,])
  probtrans_vector_1<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
  if(from_state==to_state){
    return(probtrans_vector_1)
  }
  successful_transitions_object<-successful_transitions(unique_paths_object,to_state)
  #lagged_differences_vector<-diff(spline_list[[successful_transitions_object[1]]](time))
  #integrand_1<-survival_function*c(0,lagged_differences_vector)
  
  
  #successful_transitions_object<-successful_transitions(unique_paths_object,to_state)
  #row_of_tmat_of_current_state<-which(tmat==successful_transitions_object[1],arr.ind = T)[1]
  #competing_transitions<-na.omit(tmat[row_of_tmat_of_current_state,])
  #probtrans_vector_1<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
  
  
  for(i in 1:length(successful_transitions_object)){
    lagged_differences_vector<-c(diff(spline_list[[successful_transitions_object[i]]](time)))
    
    column_of_tmat_with_successful_trans<-which(tmat==successful_transitions_object[i],arr.ind = T)[2]
    name_of_next_state<-colnames(tmat)[column_of_tmat_with_successful_trans]
    if(length(na.omit(tmat[name_of_next_state,]))==0) {
      probtrans_vector_2<-rep(1,length(time))
    }else{
      competing_transitions<-na.omit(tmat[name_of_next_state,])
      probtrans_vector_2<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
    }
    
    probtrans_vector_1<-convolute_Markov(time,lagged_differences_vector,probtrans_vector_1,probtrans_vector_2)
    
  }
  
  return(probtrans_vector_1)
}

#'Compute transition probabilities for a given starting state and target state.
#'
#'Compute transition probabilities for a given starting state and target state
#'under a homogeneous semi-Markov model
#'
#'@param tmat Transition matrix.
#'@param cumhaz 'msfit' object.
#'@param from_state Initial state.
#'@param to_state Target state.
#'@param spline_list
#'@param unique_paths_object
#'@param time
#'
#'@export
probtrans_by_convolution_semiMarkov<-function(tmat,cumhaz,from_state,to_state,spline_list,unique_paths_object,time){
  row_of_tmat_of_current_state<-which(rownames(tmat)==from_state)
  competing_transitions<-na.omit(tmat[row_of_tmat_of_current_state,])
  survival_function<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
  if(from_state==to_state){
    return(survival_function)
  }
  successful_transitions_object<-successful_transitions(unique_paths_object,to_state)
  lagged_differences_vector<-diff(spline_list[[successful_transitions_object[1]]](time))
  integrand_1<-survival_function*c(lagged_differences_vector,0)
  for(i in 1:length(successful_transitions_object)){
    column_of_tmat_with_successful_trans<-which(tmat==successful_transitions_object[i],arr.ind = T)[2]
    name_of_next_state<-colnames(tmat)[column_of_tmat_with_successful_trans]
    if(length(na.omit(tmat[name_of_next_state,]))==0) {
      integrand_2<-rep(1,length(time))
    }else if(i==length(successful_transitions_object)){
      competing_transitions<-na.omit(tmat[name_of_next_state,])
      integrand_2<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
    }else{
      competing_transitions<-na.omit(tmat[name_of_next_state,])
      survival_function<-exp(-(sapply(time,joint_cum_hazard_function,competing_transitions=competing_transitions,spline_list=spline_list)))
      lagged_differences_vector<-diff(spline_list[[successful_transitions_object[i+1]]](time))
      integrand_2<-survival_function*c(lagged_differences_vector,0)
    }
    
    integrand_1<-convolute_semiMarkov(time,integrand_1,integrand_2)
    
  }
  
  return(integrand_1)
}