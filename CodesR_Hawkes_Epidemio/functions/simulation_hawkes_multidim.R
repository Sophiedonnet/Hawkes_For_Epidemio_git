library(gtools)
#####################################
#####################################

rprior_delta <- function(M,pi_0){
  delta <- sample(c(0,1),M^2,rep=TRUE,prob=c(pi_0,1-pi_0))
  return(delta)
}
#-------------------- 

rprior_K=function(M,lambda_K){
  K_sim = rpois(M^2,lambda_K)+1
  return(K_sim)
}

#--------------------
rprior_lalpha=function(K,mu_lalpha,s_lalpha,p_Z){
  
  
  kstar= sample(1:K,1)
  Z=sample(c(0,1),K,replace=TRUE,prob=p_Z)
  Z[kstar]=1
  lalpha=(Z==0)*(-10000)+(Z==1)*rnorm(K,mean=mu_lalpha,sd=s_lalpha)
  return(lalpha)
}

#--------------------
rprior_s=function(K_p,smax_p,a_s){
  x = smax_p*rdirichlet(1,rep(a_s,K_p))
  s = c(0,cumsum(x)); 
  return(s)}  
#--------------------
dprior_delta <- function(delta,pi_0,log=TRUE){
  d = dbinom(delta,1,1-pi_0,log=TRUE)
  if(log==FALSE){d <-exp(d)}
  return(d)
}
#------------------------------
dprior_K <- function(K,lambda_K,log=TRUE){
  res=dpois(K-1,lambda_K,log=TRUE)
  if(log==FALSE){res=exp(res)}
  return(res)
}
#----------------------
#--------------------


dprior_lalpha=function(lalpha,mu_lalpha,s_lalpha,p_Z,log=TRUE){
  K=length(lalpha)
  mix1=function(x){
    if(length(x)>0){y = p_Z[2]*dnorm(x,mu_lalpha,s_lalpha)*as.numeric(x!=-10000)+p_Z[1]*as.numeric(x==-10000)}
    else{y=1}
    return(y)
  }
  mix2=function(x){y = dnorm(x,mu_lalpha,s_lalpha);return(y)}
  
  p=1/K*sum(vapply(1:K,function(kstar){prod(mix1(lalpha[-kstar]))*mix2(lalpha[kstar])},1))
  #p = prod(mix1(lalpha))
  if(log==TRUE){return(log(p))}else{return(p)}  
  }
#--------------------
dprior_s=function(s,a,log=TRUE){ ### pour s matricielle
  
  if (is.vector(s)==TRUE){
    K = length(s)-1
    n=1
    s = matrix(s,nrow=1,ncol=K+1)
  }
  if (is.matrix(s)==TRUE){
    K = dim(s)[2]-1
    n = dim(s)[1]
  }
  s_max = s[1,K+1]
  u = t(apply(s,1,diff))/s_max
  x = log(ddirichlet(u,matrix(a,n,K)))-(K-1)*log(s_max)
  if (log!=TRUE){
    x=exp(x) 
  }
  return(x)
}                                         

dprior_h_regular_s=function(h, hyperParam_prior){
  M=hyperParam_prior$M
  log_prior=list()
  log_prior$delta = dprior_delta(h$delta,hyperParam_prior$pi_0,log=TRUE)
  log_prior$K = dprior_K(h$K,hyperParam_prior$lambda_K,log=TRUE)
  log_prior$K[h$delta==0]=NA
  log_prior$lalpha = unlist(lapply(1:M^2,function(p){dprior_lalpha(h$lalpha[[p]],hyperParam_prior$mu_lalpha,hyperParam_prior$s_lalpha,hyperParam_prior$p_Z,log=TRUE)}))
  log_prior$lalpha[h$delta==0]=NA
  log_prior$h =  sum(log_prior$delta)+sum(log_prior$lalpha[h$delta==1])+sum(log_prior$K[h$delta==1])
  return(log_prior)
  
}

dprior_h_moving_s=function(h, hyperParam_prior){
  M=hyperParam_prior$M
  log_prior=list()
  log_prior$delta = dprior_delta(h$delta,hyperParam_prior$pi_0,log=TRUE)
  log_prior$K = dprior_K(h$K,hyperParam_prior$lambda_K,log=TRUE)
  log_prior$K[h$delta==0]=NA
  log_prior$lalpha = unlist(lapply(1:(M^2),function(p){dprior_lalpha(h$lalpha[[p]],hyperParam_prior$mu_lalpha,hyperParam_prior$s_lalpha,hyperParam_prior$p_Z,log=TRUE)}))
  log_prior$lalpha[h$delta==0]=NA
  log_prior$s <- unlist(lapply(1:(M^2),function(p){dprior_s(h$s[[p]],hyperParam_prior$a_s,log=TRUE)}))
  log_prior$s[h$delta==0] = NA
  log_prior$h = sum(log_prior$s[h$delta==1])+sum(log_prior$delta)+sum(log_prior$lalpha[h$delta==1])+sum(log_prior$K[h$delta==1])
  return(log_prior)
  
}
#########################################

r_h_prior_regular_s = function(hyperParam_prior,op_echan,h_vrai){
  M=hyperParam_prior$M
  smax=h_vrai$smax
  delta <- rprior_delta(M,hyperParam_prior$pi_0)
  delta[setdiff(1:M^2,op_echan$delta)]=h_vrai$delta[setdiff(1:M^2,op_echan$delta)]
  
  K =rprior_K(M,hyperParam_prior$lambda_K)
  K[delta==0]=1; 
  K[setdiff(1:M^2,op_echan$K)] = h_vrai$K[setdiff(1:M^2,op_echan$K)]

  s = lapply(1:(M^2), function(p){seq(0,smax[p],length=K[p]+1)})
  
  
  lalpha=lapply(1:(M^2),function(p){
    if(delta[p]>0){lalpha_p = rprior_lalpha(K[p],hyperParam_prior$mu_lalpha,hyperParam_prior$s_lalpha,hyperParam_prior$p_Z)}else{lalpha_p=c(-10000)}; return(lalpha_p)})
  for (p in setdiff(1:M^2,op_echan$alpha)){lalpha[[p]] = h_vrai$lalpha[[p]]}
  
  alpha=vapply(1:(M^2),function(p){exp(lalpha[[p]])})
    h = list(s=s,smax=smax,K=K,smax=smax,lalpha=lalpha, alpha=alpha,delta=delta)
  return(h)
}



# r_h_pseudoprior_regular_s = function(hyperParam_pseudoprior,op_echan,h_vrai){
#   M=hyperParam_pseudoprior$M
#   K_min = hyperParam_pseudoprior$K_min;
#   K_max = hyperParam_pseudoprior$K_max;
#   smax=h_vrai$smax
#   
#   K =rprior_K(M^2,K_min,K_max,hyperParam_pseudoprior$p_K)
#   for (p in setdiff(1:M^2,op_echan$K)){K[p] = h_vrai$K[p]}
#   
#   s = lapply(1:(M^2), function(p){seq(0,smax[p],length=K[p]+1)})
#   h = list(s=s,smax=smax,K=K)
#   
#   FF =   F_param_alpha_pseudo_all(h,hyperParam_pseudoprior)
#   
#   alpha=lapply(1:M^2,function(p){
#     alpha_p = FF$mean_pseudoprior_alpha[[p]]+ FF$sd_pseudoprior_alpha[[p]]*rnorm(K[p],0,1)
#     if(K[p]>1){alpha_p[K[p]]=0}
#     return(alpha_p)})
#   for (p in setdiff(1:M^2,op_echan$alpha)){alpha[[p]] = h_vrai$alpha[[p]]}
#   h$alpha=alpha
#   
#   return(h)
# }
# 
# 
# 
# 
# 
# 
# 
# d_h_pseudoprior_regular_s=function(h,hyperParam_pseudoprior){
#   M=hyperParam_pseudoprior$M
#   K_min = hyperParam_pseudoprior$K_min; 
#   K_max = hyperParam_pseudoprior$K_max
#   log_pseudoprior=list()
#   log_pseudoprior$K = dprior_K(h$K,K_min,K_max,hyperParam_pseudoprior$p_K,log=TRUE)
#   
#   FF =   F_param_alpha_pseudo_all(h,hyperParam_pseudoprior)
#  
#   log_pseudoprior$alpha = unlist(lapply(1:M^2,function(p){sum(dnorm(h$alpha[[p]],FF$mean_pseudoprior_alpha[[p]],FF$sd_pseudoprior_alpha[[p]],log=TRUE))}))
#   log_pseudoprior$h =  sum(log_pseudoprior$alpha)+sum(log_pseudoprior$K)
#   return(log_pseudoprior)
# }





#####################################
######################################

func_h=function(absc,alpha,s){
  f_truc <- stepfun(s,c(0, alpha, 0),right=T)
  return(f_truc(absc))
}


###############################################
###############################################

simulation_hawkes_multidim=function(Tmax,h,nu,op_affichage=0){

  ############### simulation 
  M = length(nu) # nombre de neurones
  Times=lapply(1:M,function(j){c()})
  Times_tot=Z=S=c();  
  
  
  ##############"" premier point de simulation
  i=1; 
  n=1
  grid_t = lapply(1:M,function(m){seq(0,Tmax,length=50000)})
  Lambda_star_i = sum(nu)
  S[i] =  - log(runif(1))/Lambda_star_i
  Z_n = sample(1:M,1,prob=nu) 
  Times[[Z_n]]=c(Times[[Z_n]],S[i])
  Times_tot[n] = S[i]
  Z[n] = Z_n; 
  
  ###################  ENSUITE
  while ((Times_tot[n]<Tmax)&(S[i]<Tmax)){
    #------------ ETAPE 2 of Chevallier
    i=i+1; 
    #----------- ETAPE 3 of Chevallier
    ############ majorations des lambdas sur[s_{i-1}, Tmax] ########### 
    #grid_t = lapply(1:M,function(m){seq(S[i-1],Tmax,length=50000)})
    #l = lambda_function(grid_t,h,nu,Times)
    #Lambda_star_i_1 = sum(vapply(1:M,function(m){max(l[[m]])},1))
    
    
    
    jumps = jumps_integral(Times,h$s,rep(S[i-1],M),rep(Tmax,M),h$smax)
    center_jumps_mat = center_jumps(jumps)
    
    len_Times = vapply(1:M,function(m){length(Times[[m]])},1)
    mat_centerjumps_utiles=calc_mat_absc_utiles(Times,center_jumps_mat,h$smax)
    l =  lambda_cond_multidim_mat(mat_centerjumps_utiles,h,nu)
    Lambda_star_i = sum(vapply(1:M,function(m){max(l[[m]])},1))
    
 
    #---------------- ETAPE 4 of Chevallier
    S[i] = S[i-1] - log(runif(1))/Lambda_star_i
    
    #---------------- ETAPE 5 of Chevallier
    U_i = runif(1)
    
    vect_lambda_S_i = unlist(lambda_cond_multidim(as.list(rep(S[i],M)),Times,h,nu),use.names = FALSE)

    if (U_i<= sum(vect_lambda_S_i)/Lambda_star_i){
      n = n+1; 
      Z_n = min(which(cumsum(vect_lambda_S_i)/Lambda_star_i>=U_i)) 
      Times[[Z_n]]=c(Times[[Z_n]],S[i])
      Times_tot[n] = S[i]
      if((op_affichage==1)&(n%%100==0)){print(c(n,Times_tot[n]))}
      Z[n] = Z_n; 
    }
  }
  res=list(); 
  w = which(Times_tot<Tmax)
  res$Times_tot=Times_tot[w]
  res$Z = Z[w]
  res$Times=lapply(1:M,function(m){res$Times_tot[res$Z==m]})
  res$S=S
  return(res)
}



###############################################
###############################################
# simulation_hawkes_multidim_allh=function(Tmax,h,M_h,nu,op_affichage=0){
#   
#   ### h in a decreasing function. 
#   ############"  h définie par M² alpha_{lm}exp{-beta_lm(u)} : h$alpha h$beta = de taille M*M
#   M = length(nu) # nombre de neurones
#   Times=lapply(1:M,function(j){c(0)})
#   Times_tot=Z=S=c(0); 
#   
#   ## calcul de la fonction d'intensites conditionnelle
#   i=1; 
#   n=1;
#   while ((Times_tot[n]<Tmax)&(S[i]<Tmax)){
#     #------------ ETAPE 2 of Chevallier
#     i=i+1; 
#     #----------- ETAPE 3 of Chevallier
#     ############ majorations des lambdas sur[s_{i-1}, Tmax] ########### 
#     
#     len_Times = vapply(1:M,function(m){length(Times[[m]])},1)
#     Lambda_star_i <- sum(vapply(1:M,function(m){nu[m]+sum(M_h[,m]*len_Times)},1))
#    
#     
#     
#     #---------------- ETAPE 4 of Chevallier
#     S[i] = S[i-1] - log(runif(1))/Lambda_star_i
#     
#     #---------------- ETAPE 5 of Chevallier
#     U_i = runif(1)
#     
#     vect_lambda_S_i = unlist(lambda_cond_multidim_allh(as.list(rep(S[i],M)),Times,h,nu),use.names = FALSE)
#     if (U_i<= sum(vect_lambda_S_i)/Lambda_star_i){
#       n = n+1; 
#       Z_n = min(which(cumsum(vect_lambda_S_i)/Lambda_star_i>=U_i)) 
#       Times[[Z_n]]=c(Times[[Z_n]],S[i])
#       Times_tot[n] = S[i]
#       if((op_affichage==1)&(n%%100==0)){print(c(n,Times_tot[n]))}
#       Z[n] = Z_n; 
#     }
#   }
#   res=list(); 
#   w = which(Times_tot<Tmax)
#   res$Times_tot=Times_tot[w]
#   res$Z = Z[w]
#   res$Times=lapply(1:M,function(m){res$Times_tot[res$Z==m]})
#   res$S=S
#   return(res)}
# 
# 
# ############################  
# plot_h = function(h_pop,Weights,M,smax,h_vrai,h_init_func){
#   NMC=length(Weights)
#   U = sample(1:NMC,50,replace=TRUE,prob = Weights)
#   par(mfrow=c(M,M))
#   absc=seq(0,smax[1],len=110)
#   for(p in 1:M^2){
#     ### estimations
#     h_moy=vapply(1:NMC,function(j){
#       h_chap=stepfun(h_pop[[j]]$s[[p]],c(0, h_pop[[j]]$alpha[[p]], 0),right=T)
#       estim_h_j=h_chap(absc)*Weights[j];
#       return(estim_h_j)
#     },absc)
#     h_estim=rowSums(h_moy)
#     h_med=vapply(1:NMC,function(j){
#       h_chap=stepfun(h_pop[[j]]$s[[p]],c(0, h_pop[[j]]$alpha[[p]], 0),right=T)
#       estim_h_j=h_chap(absc);
#       return(estim_h_j)
#     },absc)
#     h_estim_med=vapply(1:length(absc),function(k){u = h_med[k,]; o=order(u); ind = which(cumsum(Weights[o])>0.5)[1]; return(u[o[ind]])},1)
#     
#     #### plot
#     plot(absc,h_estim,type='l',main=p,ylim=c(min(h_vrai$alpha[[p]])*1.2,max(h_vrai$alpha[[p]])*1.2),lwd=2,xlab='')
#     for(j in U){
#       h_pop_j=stepfun(h_pop[[j]]$s[[p]],c(0, h_pop[[j]]$alpha[[p]], 0),right=T)
#       curve(h_pop_j,add=TRUE,col='grey')
#     }
#     h_p_vrai = stepfun(h_vrai$s[[p]],c(0,h_vrai$alpha[[p]],0),right=T)
#     curve(h_p_vrai,0,smax[p],col='red',add=TRUE,lty=2,lwd=2)
#     h_init_func_p = h_init_func[[p]]
#     curve(h_init_func_p,0,smax[p],col='magenta',add=TRUE,lty=3,lwd=2)
#     lines(absc,h_estim_med,lty=2,col='green',lwd=2)
#    
#   }
# }  
#   
############################  
plot_h_without_Approx = function(h_pop,Weights,M,smax,h_vrai){
  NMC=length(Weights)
  U = sample(1:NMC,50,replace=TRUE,prob = Weights)
  par(mfrow=c(M,M))
  absc=seq(0,smax[1],len=110)
  for(p in 1:M^2){
    ### estimations
    h_moy=vapply(1:NMC,function(j){
      h_chap=stepfun(h_pop[[j]]$s[[p]],c(0, h_pop[[j]]$alpha[[p]], 0),right=T)
      estim_h_j=h_chap(absc)*Weights[j];
      return(estim_h_j)
    },absc)
    h_estim=rowSums(h_moy)
    h_med=vapply(1:NMC,function(j){
      h_chap=stepfun(h_pop[[j]]$s[[p]],c(0, h_pop[[j]]$alpha[[p]], 0),right=T)
    estim_h_j=h_chap(absc);
      return(estim_h_j)
    },absc)
    h_estim_med=vapply(1:length(absc),function(k){u = h_med[k,]; o=order(u); ind = which(cumsum(Weights[o])>0.5)[1]; return(u[o[ind]])},1)
      
    #### plot
    plot(absc,h_estim,type='l',main=p,ylim=c(min(h_vrai$alpha[[p]])*1.2,max(h_vrai$alpha[[p]])*1.2),lwd=2,xlab='')
    for(j in U){
      h_pop_j=stepfun(h_pop[[j]]$s[[p]],c(0, h_pop[[j]]$alpha[[p]], 0),right=T)
      curve(h_pop_j,add=TRUE,col='grey')
    }
    h_p_vrai = stepfun(h_vrai$s[[p]],c(0,h_vrai$alpha[[p]],0),right=T)
    curve(h_p_vrai,0,smax[p],col='red',add=TRUE,lty=2,lwd=2)
    #h_init_func_p = h_init_func[[p]]
    #  curve(h_init_func_p,0,smax[p],col='magenta',add=TRUE,lty=3,lwd=2)
    lines(absc,h_estim_med,lty=2,col='green',lwd=2)
      
    }
  }   
  
  
  
  
############## 
dtriangle=function(x,a,b,c,log=FALSE){
  y = as.numeric((x>=a)&(x<c))* 2*(x-a)/((b-a)*(c-a)) + as.numeric((x<=b)&(x>=c))* 2*(b-x)/((b-a)*(b-c))
  if(log==TRUE){y=log(y)}
  return(y)
}

ptriangle=function(x,a,b,c){
  y = as.numeric((x>=a)&(x<c))* (x-a)^2/((b-a)*(c-a)) + as.numeric((x<=b)&(x>=c))* (1-(b-x)^2/((b-a)*(b-c)))+ as.numeric(x>b)
  return(y)
}


rtriangle=function(n,a,b,c){
  
  U=runif(n); 
  Fc = ptriangle(c,a,b,c)
  echan = (U<Fc)*(a+sqrt(U*(b-a)*(c-a))) + (U>=Fc)*( b - sqrt((1-U)*(b-a)*(b-c))) 
  return(echan)
}  
  
############################# 
ddirichlet=function(x,a,log=FALSE){
  if(is.vector(x)){x=matrix(x,nrow=1,ncol=length(x))}
  p=dim(x)[2]; 
  n=dim(x)[1]; 
  if(length(a)==1){a=rep(a,p)}
  a = matrix(a,ncol=p,nrow=n)
  B=prod(gamma(a))/gamma(sum(a))
  y = sum((a-1)*log(x))
  if(log==FALSE){y=exp(y)}
  return(y)
}
