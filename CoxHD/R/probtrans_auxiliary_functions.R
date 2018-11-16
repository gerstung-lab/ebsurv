#' @useDynLib ebsurv
#' @importFrom Rcpp sourceCpp

#'Create a spline approximation for each vector of cumulative hazards
#'
#'
#'@param cumhaz
#'@return A list with a spline function for each vector of cumulative hazards 
#'
#'@author rc28
#'@export

cumhaz_splines<-function(cumhaz){
  spline_list<-vector("list",length(unique(cumhaz$Haz$trans)))
  for(i in unique(cumhaz$Haz$trans)){
    cumhaz_subset<-cumhaz$Haz[cumhaz$Haz$trans==i,]
    spline_list[[i]]<-splinefun(cumhaz_subset[,"time"], cumhaz_subset[,"Haz"], method="monoH.FC")
  }
  spline_list
}

#'Find the all possible paths until absorption from a given starting state
#'
#'
#'@param from_state
#'@param tmat
#'@return A matrix where each column is a sequence of states taken by the process until absorption. 
#'There are as many columns as the number of possible paths until absorption.
#'
#'@author rc28
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

#'Create a function that returns the cumulative hazard of leaving a given state 
#'
#'
#'@param t
#'@param competing_transitions
#'@param spline_list
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
