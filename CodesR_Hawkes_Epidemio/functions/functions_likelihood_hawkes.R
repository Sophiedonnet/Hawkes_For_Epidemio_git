#library(Rcpp)
#sourceCpp('functions/functions_C.cpp')

##########################################
jumps_integral_m = function(Times_utiles,s,Tinf_m,Tmax_m,smax,m){
  M=sqrt(length(s))
  jumps_m = lapply(1:M, function(l){p = (m-1)*M+l;c(Tinf_m,Tmax_m,as.vector(outer(s[[p]],Times_utiles[[l]],'+')))})
  jumps_m=sort(unique(unlist(jumps_m,use.names = FALSE)))
  jumps_m = jumps_m[jumps_m<=Tmax_m]
  jumps_m = jumps_m[jumps_m>=Tinf_m]
  return(jumps_m)
}

jumps_integral = function(Times_utiles,s,Tinf,Tmax,smax){
  M=sqrt(length(smax))
  jumps=lapply(1:M, function(m){jumps_integral_m(Times_utiles,s,Tinf[m],Tmax[m],smax,m)})
  return(jumps)
}

center_jumps=function(jumps){
  M = length(jumps)
  C = lapply(1:M,function(m){d_m=length(jumps[[m]]);if(d_m>1){res=0.5*(jumps[[m]][1:(d_m-1)]+jumps[[m]][2:(d_m)])}; if(d_m<=1){res=c()}; return(res)})
  return(C)
}

#######################################################################
##########################################
jumps_integral_diff_m = function(Times_utiles,s,s.prime,Tinf_m,Tmax_m,smax,m){
  M=sqrt(length(s))
  jumps_m = lapply(1:M, function(l){p = (m-1)*M+l;c(Tinf_m,Tmax_m,as.vector(outer(s[[p]],Times_utiles[[l]],'+')),as.vector(outer(s.prime[[p]],Times_utiles[[l]],'+')))})
  jumps_m=sort(unique(unlist(jumps_m,use.names = FALSE)))
  jumps_m = jumps_m[jumps_m<=Tmax_m]
  jumps_m = jumps_m[jumps_m>=Tinf_m]
  return(jumps_m)
}

jumps_integral_diff = function(Times_utiles,s,s.prime,Tinf,Tmax,smax){
  M=sqrt(length(smax))
  jumps=lapply(1:M, function(m){jumps_integral_diff_m(Times_utiles,s,s.prime,Tinf[m],Tmax[m],smax,m)})
  return(jumps)
}

center_jumps=function(jumps){
  M = length(jumps)
  C = lapply(1:M,function(m){d_m=length(jumps[[m]]);if(d_m>1){res=0.5*(jumps[[m]][1:(d_m-1)]+jumps[[m]][2:(d_m)])}; if(d_m<=1){res=c()}; return(res)})
  return(C)
}

#################################################################"



calc_mat_absc_utiles= function(Times_utiles,absc,smax){
  M=sqrt(length(smax))
  f= function(p){
    l = (p-1)%%M+1
    m =(p-1)%/%M+1
   if((length(Times_utiles[[l]])>0)&(length(absc[[m]])>0)){
     # List_lm = calc_mat_absc_utiles_p_C(Times_utiles[[l]], absc[[m]], smax[p], M,l, m)
      List_lm = calc_mat_absc_utiles_p(Times_utiles[[l]], absc[[m]], smax[p])
      
      size_lm=List_lm$size_mat
      a_m =length(absc[[m]])
      mat_lm = List_lm$mat
      if (is.vector(mat_lm)){mat_lm = matrix(mat_lm,nrow=size_lm,ncol=a_m)}
      if (dim(mat_lm)[2]!= a_m){mat_lm =t(mat_lm)}
    }else{mat_lm=c()}
  return(mat_lm)}
  mat_absc_utiles=lapply(1:M^2,f)
  return(mat_absc_utiles)
}
  # ### verif taille
  

calc_mat_absc_utiles_p=function(Times_utiles_l, absc_m, smax_p){
  n_lm= length(absc_m);
  g= function(k){
    ulm_k = absc_m[k]-Times_utiles_l[Times_utiles_l<absc_m[k]];
    ulm_k = ulm_k[ulm_k<smax_p]
    if (length(ulm_k)==0){ulm_k=c(0)}
    return(c(ulm_k))
  }
  List_lm = lapply(1:n_lm,g);
  len = vapply(1:n_lm,function(k){length(List_lm[[k]])},1)
  Lmax = max(len)+1
  v=rep(0,Lmax)
  mat_lm = vapply(1:n_lm,function(k){c(List_lm[[k]][1:len[k]],rep(0,Lmax-len[k]))},v)
  res_lm=list(mat=mat_lm,size_mat=dim(mat_lm)[1])
  return(res_lm)
}



calc_mat_absc_utiles_m= function(Times_utiles,absc,smax,m0){
  M=sqrt(length(smax))
  f= function(p){
    l = (p-1)%%M+1
    m =(p-1)%/%M+1
    if (m==m0){
       List_lm = calc_mat_absc_utiles_p_C(Times_utiles[[l]], absc[[m]], smax[p], M,l, m)
      #List_lm = calc_mat_absc_utiles_p(Times_utiles[[l]], absc[[m]], smax[p])
       size_lm=List_lm$size_mat
        a_m =length(absc[[m]])
        mat_lm = List_lm$mat
        if (is.vector(mat_lm)){mat_lm = matrix(mat_lm,nrow=size_lm,ncol=a_m)}
        if (dim(mat_lm)[2]!= a_m){mat_lm =t(mat_lm)}
        }
    if(m!=m0){mat_lm=c()}
    return(mat_lm)}
  mat_absc_utiles=lapply(1:M^2,f)
  
  # ### verif taille
  
  return(mat_absc_utiles)}







########################



lambda_cond_multidim_mat = function(mat_absc_utiles,h,nu){
  M = length(nu)
  f=function(m){
    y_m = lambda_cond_multidim_mat_m(mat_absc_utiles,h,nu,m)
    return(y_m)}
  L=lapply(1:M,f)
  return(L)
}





##############################

lambda_cond_multidim_mat_m = function(mat_absc_utiles,h,nu,m){
  M =sqrt(length(h$s))
  p_2=(m-1)*M+1;
  if (is.vector(mat_absc_utiles[[p_2]])==TRUE){mat_absc_utiles[[p_2]]=matrix(mat_absc_utiles[[p_2]],nrow=1)}
  g=function(l){
    p = (m-1)*M+l
    d=dim(mat_absc_utiles[[p]])
    # changement par Fabian
    idx = findInterval(as.vector(mat_absc_utiles[[p]]), h$s[[p]] + .Machine$double.eps)
    l_p = matrix(c(0,h$alpha[[p]],0)[idx + 1], dim(mat_absc_utiles[[p]]))
    #
    #h_p = stepfun(h$s[[p]],c(0,h$alpha[[p]],0),right=TRUE)
    #l_p = sapply(1:d[2],function(k){h_p(mat_absc_utiles[[p]][,k])})
    if(is.vector(l_p)){l_p = matrix(l_p,nrow=1)}
    L_p = colSums(l_p)
    return(L_p)
  }
  
  dim_mat = vapply(1:M,function(l){p = (m-1)*M+l; return(sum(dim(mat_absc_utiles[[p]])))},1)
  L_m=sapply(which(dim_mat>0),g)
  if (is.vector(L_m)==TRUE){L_m=matrix(L_m,nrow=1,ncol=M)}
  R_m =rowSums(L_m)
  y_m=0.5*(abs(nu[m]+R_m)+nu[m]+R_m)
  y_m[y_m<0]=0;
  return(y_m)

}

lambda_cond_multidim=function(t,Times,h,nu){
  ## calcul de la fonction d'intensites conditionnelle
  ## alpha positifs ou negatifs.
  M = length(Times)
  L = function(m){lambda_cond_multidim_m(t[[m]],Times,h,nu,m)}
  res = lapply(1:M,L)  
  return(res)  
}

############### function lambda naive
lambda_function = function(grid_t,Times,h,nu){
  M = length(grid_t); 
  res=lapply(1:M,function(m){
    res_m=rep(nu[m],length(grid_t[[m]]))
    for (l in 1:M){
      p = (m-1)*M+l; 
      hlm = stepfun(h$s[[p]],c(0,h$alpha[[p]],0),right=TRUE)
      Diff = outer(grid_t[[m]],Times[[l]],FUN=function(x,y){d = (x-y); return(hlm(d))})  
      res_m =res_m+ rowSums(Diff)
    }
    return(res_m)})
  return(res)  
}




###############################################
lambda_cond_multidim_m=function(t_m,Times,h,nu,m){
  
  M=sqrt(length(h$s))
  Mat_y=matrix(0,M+1,length(t_m)) 
  Mat_y[1,]=nu[m]
  
  Remp_M=function(l){
    p= (m-1)*M+l; 
    Times_l =Times[[l]]
    
    #Times_utiles_l_m = Times_l
    Mat_y_l= unlist(lapply(1:length(t_m), function(k){
      Times_utiles_l_m=Times_l[((t_m[k]-Times_l)<h$smax[p]) &((t_m[k]-Times_l)>0) ];
      f_truc <- stepfun(h$s[[p]],c(0, h$alpha[[p]], 0),right=T)
      sum(f_truc(t_m[k]-Times_utiles_l_m))}),use.names = FALSE)
  }
  Mat_y[2:(M+1),] = matrix(unlist(lapply(1:M,Remp_M),use.names = FALSE),M,byrow=TRUE)
  
  
  y = colSums(Mat_y)
  res = 0.5*(y+abs(y))
  return(res)
}




###############################
int_lambda_multidim=function(mat_centerjumps_utiles,jumps,h,nu){
  l = lambda_cond_multidim_mat(mat_centerjumps_utiles,h,nu)
  int_lambda_m = unlist(lapply(1:M,function(m){sum(l[[m]]*diff(jumps[[m]]))}),use.names = FALSE)
  return(int_lambda_m)
}

  
##########################################
log_likelihood_hawkes_multidim=function(mat_obs_utiles,mat_centerjumps_utiles,jumps,h,nu,op='sum'){
  
  I = int_lambda_multidim(mat_centerjumps_utiles,jumps,h,nu)
  eval_lambda = lambda_cond_multidim_mat(mat_obs_utiles,h,nu)
  
  g=function(m){if(length(eval_lambda[[m]])==0){return(1)} else{return(eval_lambda[[m]])}}
  eval_lambda=lapply(1:M,g)
  L = vapply(1:M, function(m){sum(log(eval_lambda[[m]]))},1)-I
  if (op=='sum'){L = sum(L)};
  return(L)
}









int_lambda_multidim_m=function(mat_centerjumps_utiles,jumps,h,nu,m){
  ## calcul de l'integrale de l'intensite sur Tinf Tmax
  d_m=length(jumps[[m]])
  l_m = lambda_cond_multidim_mat_m(mat_centerjumps_utiles,h,nu,m)
  int_lambda_m = sum(l_m*diff(jumps[[m]]))
  return(int_lambda_m)
}









log_likelihood_hawkes_multidim_m=function(mat_obs_utiles,mat_centerjumps_utiles,jumps,h,nu,m){
  I_m = int_lambda_multidim_m(mat_centerjumps_utiles,jumps,h,nu,m)
  eval_lambda_m = lambda_cond_multidim_mat_m(mat_obs_utiles,h,nu,m)
  L_m = sum(log(eval_lambda_m))-I_m
  return(L_m)
}

#
# calc_mat_absc_utiles_non_C= function(Times_utiles,absc,smax){
#   M=sqrt(length(smax))
#   f= function(p){
#     l = (p-1)%%M+1
#     m =(p-1)%/%M+1
#     mat_lm = calc_mat_absc_utiles_p(Times_utiles[[l]], absc[[m]], smax, M,l, m)
#     if (dim(mat_lm)[2]!= length(absc[[m]])){mat_lm =t(mat_lm) }
#     if (dim(mat_lm)[2]!= length(absc[[m]])){print('Erreur taille')}
#     return(mat_lm)}
#   mat_absc_utiles=lapply(1:M^2,f)
# 
#    # ### verif taille
# 
#   return(mat_absc_utiles)}


max_vect=function(v){
  if (length(v)==0){res=NA}else{res=max(v)}
  return(res)
}

min_vect=function(v){
  if (length(v)==0){res=NA}else{res=min(v)}
  return(res)
}

max_list = function(L){
  len_L=length(L); 
  m_L = vapply(1:len_L,function(k){max_vect(L[[k]])},1)
  return(m_L)
}


#

define_blocks=function(Times_obs,len,first=1){
  M = length(Times_obs)
  neuron_times_obs = unlist(lapply(1:M,function(m){rep(m,length(Times_obs[[m]]))}))
  Times_obs_ord =unlist(Times_obs)
  ord_obs = order(Times_obs_ord)
  Times_obs_ord = Times_obs_ord[ord_obs]
  neuron_times_obs_ord = neuron_times_obs[ord_obs]
  #
  
  neuron_vu=c()
  istart=1
  while(length(neuron_vu)<M){
    istart=istart+1
    neuron_vu=c(neuron_vu,neuron_times_obs_ord[istart])
    neuron_vu=unique(neuron_vu)
  }
  
  istart=max(istart,first)
  blocks=rep(1,istart)
  length_blocks = len
  U=(length(Times_obs_ord)-istart)%/%length_blocks
  blocks=c(blocks,c(sort(rep(2:(U+1),length_blocks)),rep(U+2,length(Times_obs_ord)-istart-U*length_blocks)))
  return(blocks)
}



