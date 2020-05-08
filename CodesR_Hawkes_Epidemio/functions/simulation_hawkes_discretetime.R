
# Simulate data
simulate_data_multidim_dt = function(Tmax,h,nu){
  
  # Initialise
  M = length(nu)
  K = h$K
  Times=lapply(1:M,function(j){cbind(1:Tmax, rep(0,Tmax))})
  Intensities=lapply(1:M,function(j){rep(0,Tmax)})
  
  # First period of simulation
  i=1
  for (m in 1:M){
    Intensities[[m]][i] = sum(nu)
    Times[[m]][i,2] = rpois(1, lambda=Intensities[[m]][i])
  }
  
  for (i in 2:Tmax){
    absc_utiles = calc_absc_utiles_dt(Times,i,h$smax)
    l= lambda_cond_multidim_dt(absc_utiles,h,nu)
    for (m in 1:M){
      Intensities[[m]][i] = sum(unlist(l[[m]]))
      Times[[m]][i,2] = rpois(1, lambda=Intensities[[m]][i])
    }
  }
  
  return(Times)
  
}

# set.seed(3240)
# set.seed(46285)
# Times_multidim = simulate_data_multidim_dt(80, h_vrai, nu)


# Simulate data
simulate_data_dt = function(Tmax,h,nu){
  
  # Initialise
  K = h$K
  sim_data = cbind(seq(1,Tmax,1), rep(0,Tmax), rep(0,Tmax))

  # First period of simulation
  sim_data[1,2] = nu
  sim_data[1,3] = rpois(1, lambda=sim_data[1,2])
  
  for (i in 2:Tmax){
    sim_decay_times = as.numeric(i-sim_data[sim_data[,1]< i & (i-sim_data[,1]) <= K & sim_data[,3] > 0, 1])
    sim_decay_times_counts = as.numeric(sim_data[sim_data[,1] < i & (i-sim_data[,1]) <= K & sim_data[,3] >0, 3])
    if(length(sim_decay_times) > 0){
      sim_decay_vals = sapply(sim_decay_times, function(x) h$alpha[[1]][x])
      sim_decay_vals_counts = sim_decay_vals*sim_decay_times_counts
    }else{
      sim_decay_vals_counts = 0
    }
    
    sim_data[i,2] = nu + sum(sim_decay_vals_counts)
    sim_data[i,3] = rpois(1, lambda=sim_data[i,2])
  }
  
  return(list(Times = sim_data[,c(1,3)]))
  
}

# set.seed(3240)
# set.seed(46285)
# Times_unidim = simulate_data_dt(80, h_vrai, nu)

#####################################
#####################################

#--------------------
# dprior_delta <- function(delta,pi_0,log=TRUE){
#   d = dbinom(delta,1,1-pi_0,log=TRUE)
#   if(log==FALSE){d <-exp(d)}
#   return(d)
# }
#------------------------------
# dprior_K <- function(K,lambda_K,log=TRUE){
#   res=dpois(K-1,lambda_K,log=TRUE)
#   if(log==FALSE){res=exp(res)}
#   return(res)
# }
#----------------------
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
# ddirichlet=function(x,a,log=FALSE){
#   if(is.vector(x)){x=matrix(x,nrow=1,ncol=length(x))}
#   p=dim(x)[2]; 
#   n=dim(x)[1]; 
#   if(length(a)==1){a=rep(a,p)}
#   a = matrix(a,ncol=p,nrow=n)
#   B=prod(gamma(a))/gamma(sum(a))
#   y = sum((a-1)*log(x))
#   if(log==FALSE){y=exp(y)}
#   return(y)
# }
# 
# dprior_s=function(s,a,log=TRUE){ ### pour s matricielle
#   
#   if (is.vector(s)==TRUE){
#     K = length(s)-1
#     n=1
#     s = matrix(s,nrow=1,ncol=K+1)
#   }
#   if (is.matrix(s)==TRUE){
#     K = dim(s)[2]-1
#     n = dim(s)[1]
#   }
#   s_max = s[1,K+1]
#   u = t(apply(s,1,diff))/s_max
#   x = log(ddirichlet(u,matrix(a,n,K)))-(K-1)*log(s_max)
#   if (log!=TRUE){
#     x=exp(x) 
#   }
#   return(x)
# }                                         

dprior_h_regular_s=function(h, hyperParam_prior){
  M=hyperParam_prior$M
  log_prior=list()
  #log_prior$delta = dprior_delta(h$delta,hyperParam_prior$pi_0,log=TRUE)
  #log_prior$K = dprior_K(h$K,hyperParam_prior$lambda_K,log=TRUE)
  #log_prior$K[h$delta==0]=NA
  log_prior$lalpha = unlist(lapply(1:M^2,function(p){dprior_lalpha(h$lalpha[[p]],hyperParam_prior$mu_lalpha,hyperParam_prior$s_lalpha,hyperParam_prior$p_Z,log=TRUE)}))
  log_prior$lalpha[h$delta==0]=NA
  log_prior$h =  sum(log_prior$lalpha[h$delta==1])
  #log_prior$h =  sum(log_prior$delta)+sum(log_prior$lalpha[h$delta==1])+sum(log_prior$K[h$delta==1])
  return(log_prior)
}
