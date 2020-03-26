
##############################################################################
######################  POUR Reversible Jumps ################################
# ################################################################################
# prob_moves=function(K_min,K_max,p_K,fixed_s=TRUE,op_echan_K=TRUE){
#   if(fixed_s==TRUE){
#     vec_k <- K_min:K_max; 
#     L <- length(vec_k)
#     r <- sum(vec_k*p_K)
#     d_k=b_k=vec_k; for(k in 1:L){b_k[k] <- min(1,r/vec_k[k])}
#     b_k[L] <- 0
#     for(k in 1:L){d_k[k] <- min(1,(vec_k[k]-1)/r)}
#     d_k[1] <- 0;
#     C <- 0.7/max(b_k+d_k)
#     death_k <- C*d_k
#     birth_k <- C*b_k
#     #rho_k <- 1-death_k-birth_k
#     prob <- cbind(birth_k,death_k)} ### hauteurs naissance mort
#   
#   if(fixed_s==FALSE){
#     if(op_echan_K==TRUE){
#       vec_k <- K_min:K_max; 
#       L <- length(vec_k)
#       r <- sum(vec_k*p_K)
#       d_k=b_k=vec_k; for(k in 1:L){b_k[k] <- min(1,r/vec_k[k])}
#       b_k[L] <- 0
#       for(k in 1:L){d_k[k] <- min(1,(vec_k[k]-1)/r)}
#       d_k[1] <- 0;
#       C <- 0.7/max(b_k+d_k)
#       death_k <- C*d_k
#       birth_k <- C*b_k
#       #rho_k <- 0.5*(1-death_k-birth_k)
#       #eta_k <- rho_k;
#       prob <- cbind(birth_k,death_k) #### largeur hauteur naissance mort
#     }
#     if(op_echan_K==FALSE){
#       prob <-matrix(0,ncol=2,nrow=length(K_min:K_max)) 
#     }
#   }
#   
#   return(prob)
# }

birth_regular_s = function(h_small,p,pi_moves,hyperParam_prior){

  
  h_big = h_small
  smax = h_small$smax  
  ############simulation de la nouvelle fonction avec un pallier supplementaire
  ### simulation de la place du nouveau pallier
  K_p = h_small$K[p] 
  alpha_p = h_small$alpha[[p]]
  h_big$K[p] = K_p + 1
  h_big$s[[p]] = seq(0,smax[p],length = K_p +2)

  k_new =sample(1:(K_p),1)
  u = c(0,alpha_p)
  xi = sample(c(0,1),1)
  mu_k_new =xi*u[k_new+1]+(1-xi)*u[k_new]
  sigma = min(abs(u[k_new+1]-u[k_new]),0.1);
  alpha_k_new = mu_k_new + sigma*rnorm(1);
  khi1 = alpha_k_new-u[k_new+1]
  khi2 = alpha_k_new-u[k_new]
  if(k_new==1){h_big$alpha[[p]]=c(alpha_k_new,alpha_p)}
  if(k_new>1){h_big$alpha[[p]]=c(alpha_p[1:(k_new-1)],alpha_k_new,alpha_p[(k_new):K_p])}
  Jacob = 1;  
  log_prob_proposal =log(dnorm(khi1,0,sigma)*0.5+dnorm(khi2,0,sigma)*0.5)+log(pi_moves[2])-log(pi_moves[1])+ log(1/(K_p)) 
  log_prob_passage = log(Jacob)-log_prob_proposal; 
  
  
  return(list(log_prob_passage=log_prob_passage,h=h_big,where=k_new))
  
  
  
  
}

death_regular_s = function(h_big,p,pi_moves,hyperParam_prior){
  K_p = h_big$K[p]
  alpha_p = h_big$alpha[[p]]
  h_small=h_big
  
  smax = h_big$smax  
  h_small$K[p] = K_p - 1
  h_small$s[[p]] = seq(0,smax[p],length = K_p)
  
  k_out =sample(1:max((K_p-1),1),1)
  k=k_out-1; 
  
  h_small$alpha[[p]] = alpha_p[-k_out]
  u = c(0,alpha_p)
  sigma = max(abs(u[k+3]-u[k+1]),0.1);
  khi1 = u[k+2]- u[k+3]
  khi2 = u[k+2]-u[k+1]
  d = log(dnorm(khi1,0,sigma)*0.5+0.5*dnorm(khi2,0,sigma));  
  
  Jacob=1
  log_prob_proposal =d+log(pi_moves[2])-log(pi_moves[1])+ log(1/(K_p-1)) 
  log_prob_passage = log(Jacob)-log_prob_proposal; 
  log_prob_proposal=-log_prob_proposal
  return(list(log_prob_passage=log_prob_passage,h=h_small,where=k_out))
 
}
#####################################################################################################################
####################### with moving s
###################################################################################################################


birth_moving_s = function(h_small,p,pi_moves,hyperParam_prior){

  
  h_big = h_small
  smax = h_small$smax  
  ############simulation de la nouvelle fonction avec un pallier supplementaire
  ### simulation de la place du nouveau pallier
  K_p = h_small$K[p] 
  alpha_p = h_small$alpha[[p]]
  h_big$K[p] = K_p + 1
  
  
  k_new =sample(1:(K_p),1)
  u = c(0,alpha_p)
  xi = sample(c(0,1),1)
  mu_k_new =xi*u[k_new+1]+(1-xi)*u[k_new]
  sigma = min(abs(u[k_new+1]-u[k_new]),0.1);
  alpha_k_new = mu_k_new + sigma*rnorm(1);
  khi1 = alpha_k_new-u[k_new+1]
  khi2 = alpha_k_new-u[k_new]
  if(k_new==1){h_big$alpha[[p]]=c(alpha_k_new,alpha_p)}
  if(k_new>1){h_big$alpha[[p]]=c(alpha_p[1:(k_new-1)],alpha_k_new,alpha_p[(k_new):K_p])}
  h_big$s[[p]] = sort(c(h_small$s[[p]],runif(1,h_small$s[[p]][k_new],h_small$s[[p]][k_new+1])))
  
  
  Jacob = 1;  
  log_prob_proposal =log(dnorm(khi1,0,sigma)*0.5+dnorm(khi2,0,sigma)*0.5)+log(pi_moves[2])-log(pi_moves[1])+ log(1/(K_p))-log(h_small$s[[p]][k_new+1]-h_small$s[[p]][k_new]) 
  log_prob_passage = log(Jacob)-log_prob_proposal; 
  
  
  return(list(log_prob_passage=log_prob_passage,h=h_big,where=k_new))
}


death_moving_s = function(h_big,p,pi_moves,hyperParam_prior){
  
  
  K_p = h_big$K[p]
  alpha_p = h_big$alpha[[p]]
  h_small=h_big
  
  smax = h_big$smax  
  h_small$K[p] = K_p - 1
  
  
  k_out =sample(1:max((K_p-1),1),1)
  k=k_out-1; 
  
  h_small$alpha[[p]] = alpha_p[-k_out]
  h_small$s[[p]] = h_big$s[[p]][-(k_out+1)]
  u = c(0,alpha_p)
  sigma = max(abs(u[k+3]-u[k+1]),0.1);
  khi1 = u[k+2]- u[k+3]
  khi2 = u[k+2]-u[k+1]
  d = log(dnorm(khi1,0,sigma)*0.5+0.5*dnorm(khi2,0,sigma));  
  
  Jacob=1
  log_prob_proposal =d+log(pi_moves[2])-log(pi_moves[1])+ log(1/(K_p-1)) - log(h_big$s[[p]][k_out+2]-h_big$s[[p]][k_out])
  log_prob_passage = log(Jacob)-log_prob_proposal; 
  log_prob_proposal=-log_prob_proposal
  return(list(log_prob_passage=log_prob_passage,h=h_small,where=k_out))
  
  
}






