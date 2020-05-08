
#############################################################

# Calculate decay times
calc_absc_utiles_dt= function(Times_utiles,absc,smax){
  M=sqrt(length(smax))
  f= function(p){
    l = (p-1)%%M+1
    m =(p-1)%/%M+1
    if(is.list(absc)==TRUE){
      absc = as.matrix(absc[[m]][,1])
    }
    if((length(Times_utiles[[l]][,1])>0)&(length(absc)>0)){
      List_lm = calc_absc_utiles_dt_p(Times_utiles[[l]], absc, smax[p])
    }else{
      List_lm=list()
      }
    return(List_lm)
    }
  absc_utiles=lapply(1:M^2,f)
  return(absc_utiles)
}

# For each combination of variates
calc_absc_utiles_dt_p=function(Times_utiles_l, absc_m, smax_p){
  n_lm= length(absc_m)
  decay_times = list()
  decay_times_counts = list()
  
  decay = function(k){
    decay_times_k = absc_m[k] - Times_utiles_l[(Times_utiles_l[,1] < absc_m[k]) & ((absc_m[k]-Times_utiles_l[,1]) <= smax_p) & (Times_utiles_l[,2] > 0), 1]
   if (length(decay_times_k)==0){decay_times_k=c(0)}
    decay_times[[k]] = decay_times_k
    return(decay_times_k)
  }
  decay_times = lapply(1:n_lm, decay)
  
  decay_counts = function(k){
    decay_times_counts_k = Times_utiles_l[(Times_utiles_l[,1] < absc_m[k]) & ((absc_m[k]-Times_utiles_l[,1]) <= smax_p) & (Times_utiles_l[,2] >0), 2]
    if (length(decay_times_counts_k)==0){decay_times_k=c(0)}
    decay_times_counts[[k]] = decay_times_counts_k
    return(decay_times_counts_k)
  }
  decay_times_counts = lapply(1:n_lm, decay_counts)
  
  return(list(decay_times=decay_times, decay_times_counts=decay_times_counts))
}

#############################################################

# Calculate conditional intensity function
lambda_cond_multidim_dt = function(absc_utiles,h,nu){
  M = length(nu)
  f=function(m){
    y_m = lambda_cond_multidim_dt_m(absc_utiles,h,nu,m)
    return(y_m)}
  L=lapply(1:M,f)
  return(L)
}

# For each combination of variates
lambda_cond_multidim_dt_m = function(absc_utiles,h,nu,m){
  M =sqrt(length(h$s))
  g=function(l){
    p = (m-1)*M+l
    decay_times = absc_utiles[[p]][[1]]
    decay_times_counts = absc_utiles[[p]][[2]]
    if(length(decay_times) > 0){
      decay_vals = lapply(decay_times, function(x) h$alpha[[p]][x])
      decay_vals_total = Map('*',decay_vals,decay_times_counts)
    }else{
      decay_vals_total = list()
    }
    L_m = sapply(1:length(decay_vals_total), function(x) sum(decay_vals_total[[x]], na.rm=TRUE))
    
    return(L_m)
  }
  
  L_m=lapply(1:M,g)
  L_m=matrix(unlist(L_m), nrow=length(unlist(L_m)/M), ncol=M)
  if (is.vector(L_m)==TRUE){L_m=matrix(L_m,nrow=1,ncol=M)}
  R_m =rowSums(L_m)
  y_m=nu[m] + R_m
  return(y_m)
  
}

#############################################################

# Log likelihood 
log_likelihood_multidim_dt=function(data,h,nu,op='sum'){
  M = length(nu)
  eval_lambda = lambda_cond_multidim_dt(data$obs_utiles,h,nu)
  L= vapply(1:M,function(m){sum(dpois(data$Times_obs[[m]][,2],lambda=eval_lambda[[m]],log=TRUE))},1)
  if (op=='sum'){L = sum(L)};
  return(L)
}


# Log likelihood- intervention model
log_likelihood_multidim_inter_dt=function(data,h,nu,beta,op='sum'){
  
  M = length(nu)
  h_2 = h
  for (p in 1:M) {
    h_2$alpha[[p]] =  h$alpha[[p]] * beta
  }
  
  eval_lambda_1 = lambda_cond_multidim_dt(data$obs_utiles_1,h,nu)
  eval_lambda_2 = lambda_cond_multidim_dt(data$obs_utiles_2,h_2,nu)
  L= vapply(1:M,function(m){sum(dpois(data$Times_obs_1[[m]][,2],lambda=eval_lambda_1[[m]],log=TRUE)) + 
                            sum(dpois(data$Times_obs_2[[m]][,2],lambda=eval_lambda_2[[m]],log=TRUE))},1)
  if (op=='sum'){L = sum(L)};
  return(L)
}

