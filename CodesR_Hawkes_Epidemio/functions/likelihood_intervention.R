
log_likelihood_hawkes_multidim_positive_intensity_2=function(data,mat_obs_utiles_1,mat_obs_utiles_2,
                                                             idx_1,idx_2,J_1,J_2,h,nu,beta,Lag,op='sum'){
  
  # Log-likelihood with 2 additional parameters:
  # beta < 1: scalar multiplying the kernel for t > tau + L (tau is a change point)
  # Lag: Lag after time tau
  # need to add in data 2 mat_obs_utiles: before and after tau + L
  
  Tinter = min(data$tau+Lag,data$Tmax)
  
  h_beta = duplicate(h)
  for (p in  1:M) {h_beta$alpha[[p]] =  h$alpha[[p]] * beta}
  
  M = length(nu)
  I_1 = int_lambda_multidim_positive_intensity(data$Tinf,Tinter,h,J_1,nu)
  I_2 = int_lambda_multidim_positive_intensity(Tinter,data$Tmax,h_beta,J_2,nu)
  
  eval_lambda_1 = lambda_cond_multidim_mat_positive_intensity(mat_obs_utiles_1,idx_1,h,nu)
  eval_lambda_2 = lambda_cond_multidim_mat_positive_intensity(mat_obs_utiles_2,idx_2,h_beta,nu)
  
  #print(c(I_1, I_2, sum(log(eval_lambda_1[[1]])),sum(log(eval_lambda_2[[1]]))))
  
  L= vapply(1:M,function(m){sum(log(eval_lambda_1[[m]])) + sum(log(eval_lambda_2[[m]])) - I_1[m] - I_2[m]},1)
  if (op=='sum'){
    L = sum(L)
  };
  return(L)
}



# ##############################
# 
# func_idx_p=function(mat_absc_utiles,s,p){
#   # changement par Fabian
#   idx_p = findInterval(as.vector(mat_absc_utiles[[p]]), s[[p]] + .Machine$double.eps)
#   return(idx_p)
# }
# 
# func_idx=function(mat_absc_utiles,s){
#   P=length(s)
#   return(lapply(1:P,function(p){func_idx_p(mat_absc_utiles,s,p)})) 
# }
#   
# 
# lambda_cond_multidim_mat_positive_intensity_after_m = function(mat_absc_utiles,idx,h,nu,beta,m){
#   M =sqrt(length(h$s))
#   p_2=(m-1)*M+1;
#   if (is.vector(mat_absc_utiles[[p_2]])==TRUE){mat_absc_utiles[[p_2]]=matrix(mat_absc_utiles[[p_2]],nrow=1)}
#   g=function(l){
#     p = (m-1)*M+l
#     d=dim(mat_absc_utiles[[p]])
#     l_p = matrix(c(0,h$alpha[[p]],0)[idx[[p]] + 1], dim(mat_absc_utiles[[p]]))
#     if(is.vector(l_p)){l_p = matrix(l_p,nrow=1)}
#     L_p = colSums(l_p)
#     return(L_p)
#   }
#   
#   dim_mat = vapply(1:M,function(l){p = (m-1)*M+l; return(sum(dim(mat_absc_utiles[[p]])))},1)
#   L_m=sapply(which(dim_mat>0),g)
#   if (is.vector(L_m)==TRUE){L_m=matrix(L_m,nrow=1,ncol=M)}
#   R_m =rowSums(L_m)
#   y_m=nu[m]+ beta * R_m
#   return(y_m)
#   
# }
# 
# 
# lambda_cond_multidim_mat_positive_intensity = function(mat_absc_utiles,idx,h,nu){
#   M <- length(nu); 
#   return(lapply(1:M,function(m){lambda_cond_multidim_mat_positive_intensity_m(mat_absc_utiles,idx,h,nu,m)}))
# }
# 
# 
# lambda_cond_multidim_allh_m = function(absc_m,Times,h,nu,m){
#   
#   
#   F_l <- function(l){
#     p = (m-1)*M+l
#     D <- outer(absc_m,Times[[l]],'-')
#     D[D<0] <- Inf
#     return(rowSums(h[[p]](D)))
#   }
#   u = vapply(1:M,F_l,absc_m); 
#   if(is.vector(u)){u=matrix(u,nrow=1)}
#  L_m <- nu[m]+rowSums(u)
#  return(L_m)
# }
# 
# lambda_cond_multidim_allh = function(absc,Times,h,nu){
#   
#   L <- lapply(1:M,function(m){lambda_cond_multidim_allh_m(absc[[m]],Times,h,nu,m)})
#   return(L)
#   
# }
#   
# calc_Jlm = function(Times_l,slm,Tinf_m,Tmax_m){
#   
#   ### s est liste de taille l*m
#   Klm <- length(slm)-1; 
#   Jlm <- rep(0,Klm)
# 
#   J_lm <- vapply(2:(Klm+1),function(k){
#     b.g = Times_l+slm[k-1] # bord gauche de l'intervalle ? int?grer
#     b.d = Times_l+slm[k] # bord gauche de l'intervalle ? int?grer
#     Mklm <- as.numeric((b.g>Tinf_m)&(b.d<Tmax_m))
#     Nklm <- as.numeric((b.g<Tinf_m)&(b.d>Tmax_m))
#     Dklm <- as.numeric((b.d>Tmax_m)&((b.g>Tinf_m))&(b.g<Tmax_m))
#     Gklm <- as.numeric((b.g<Tinf_m)&((b.d>Tinf_m))&(b.d<Tmax_m))
#     Jlm_k<- sum(Mklm*(b.d-b.g)) + sum((b.d-Tinf_m)*Gklm) + sum((Tmax_m-b.g)*Dklm)+sum((Tmax_m-Tinf_m)*Nklm)
#     return(Jlm_k)},1)
#   return(J_lm)  
# }
# 
# calc_J  = function(Times_utiles,s,Tinf,Tmax){
#   M <- length(Times_utiles)
#   J <- lapply(1:M^2, function(p){
#     l <- (p-1)%%M+1
#     m <- (p-1)%/%M+1
#     Jlm=calc_Jlm(Times_utiles[[l]],s[[p]],Tinf[m],Tmax[m]); return(Jlm)})
#   return(J)
# }
# 
# 
# 
# 
# 
# 
# int_lambda_multidim_positive_intensity_m= function(Tinf,Tmax,h,J,nu,m){
#   M <- length(nu); 
#   I_m = nu[m]*(Tmax[m]-Tinf[m])
#   ### Jlm contains all the integrales (Jklm)_(k=1...Klm)
#   F_l = function(l){
#     p = (m-1)*M+l; 
#     Slm = h$alpha[[p]]*J[[p]]; 
#     return(sum(Slm))
#   }
#   I_m = I_m+sum(vapply(1:M,F_l,1))
#   return(I_m)
# }
# 
# 
# 
# 
# int_lambda_multidim_positive_intensity= function(Tinf,Tmax,h,J,nu){
#   I=vapply(1:M,function(m){int_lambda_multidim_positive_intensity_m(Tinf,Tmax,h,J,nu,m)},0)
#   return(I)
# }
# 
# 
# 
# log_likelihood_hawkes_multidim_positive_intensity_m=function(data,idx,J,h,nu,m){
#   I_m = int_lambda_multidim_positive_intensity_m(data$Tinf,data$Tmax,h,J,nu,m)
#   eval_lambda_m = lambda_cond_multidim_mat_positive_intensity_m(data$mat_obs_utiles,idx,h,nu,m)
#   L_m = sum(log(eval_lambda_m))-I_m
#   return(L_m)
# }
# 
# log_likelihood_hawkes_multidim_positive_intensity=function(data,idx,J,h,nu,op='sum'){
#   M = length(nu)
#   I = int_lambda_multidim_positive_intensity(data$Tinf,data$Tmax,h,J,nu)
#   eval_lambda = lambda_cond_multidim_mat_positive_intensity(data$mat_obs_utiles,idx,h,nu)
#   L= vapply(1:M,function(m){sum(log(eval_lambda[[m]]))-I[m]},1)
#   if (op=='sum'){L = sum(L)};
#   return(L)
# }


