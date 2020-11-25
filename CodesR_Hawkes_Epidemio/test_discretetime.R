
rm(list = ls())


library(parallel)
library(mvtnorm)
library(ggplot2)
library(GGally)
#library(mcmcse)


#source('functions/simulation_hawkes_multidim.R')
#source('functions/functions_likelihood_hawkes.R')
#source('functions/functions_likelihood_hawkes_positive_intensity.R')
#source('functions/RJ_Kernel_regular_s_positive_intensity.R')
#source('functions/RJ_Kernel_moving_s_positive_intensity.R')
source('functions/simulation_hawkes_discretetime.R')
source('functions/functions_likelihood_discretetime.R')
source('functions/RJ_Kernel_regular_discretetime.R')
source('functions/RJ_Kernel_moving_discretetime.R')


op.plot= TRUE

#-----------------------------------
#-------- PARAM de simulation 
#-----------------------------------


M = 2
smax = rep(7,M^2)
nu = rep(1,M)
beta = rep(0.9,M)

s = alpha = lalpha = gamma = list()

#l = 1
#m = 1
for (m in 1:M){
  for (l in 1:M){
    p = (m-1)*M+l
    s[[p]] = c(0,3,5,smax[p]) 
    alpha[[p]] <- dnorm(s[[p]][-1],3,3)
    gamma[[p]] = c(1,alpha[[p]][2]/alpha[[p]][1], alpha[[p]][3]/alpha[[p]][1])
    #alpha[[p]] <- dnorm(s[[p]][-1],1.2,1)
    #gamma[[p]] = c(1,alpha[[p]][2]/alpha[[p]][1], alpha[[p]][3]/alpha[[p]][1])
    #alpha[[p]] <- dnorm(s[[p]][-1],2,1)/sum(dnorm(s[[p]][-1],2,1))
    #gamma[[p]] = c(1,alpha[[p]][2]/alpha[[p]][1], alpha[[p]][3]/alpha[[p]][1], alpha[[p]][4]/alpha[[p]][1])  }
  }
}

for (p in 1:M^2){
  K = sapply(1:M^2,function(p){length(alpha[[p]])});
  s = lapply(1:M^2,function(p){s[[p]]})
  s_diffs = lapply(1:M^2,function(p){diff(s[[p]])})
  delta = vapply(1:M^2,function(p){as.numeric(max(alpha[[p]])>0)},1)
}
  
h_vrai=list(alpha=alpha,s=s,smax=smax,s_diffs=s_diffs,K=K,delta=delta, beta=beta, gamma=gamma)
nu_vrai=nu;
theta_vrai = list(h_vrai=h_vrai,nu_vrai=nu_vrai)

#------- plot de $h_p$
par(mfrow=c(M,M))
I=matrix(0,M,M)
for (p in 1:M^2){
  l = (p-1)%%M+1
  m = (p-1)%/%M+1
  print(c(l,m))
    
  h_p  <- stepfun(h_vrai$s[[p]],c(0,h_vrai$alpha[[p]],0),right=T)
  curve(h_p,0,smax[p],col='red',lwd=2)
  I[l,m] =sum(diff(h_vrai$s[[p]])*(h_vrai$alpha[[p]]))
  curve(h_p,0,8)
}
eigen(I)

#-----------------------------------
#------ SIMULATION des données 
#-----------------------------------

sim_length = 20
Times = simulate_data_multidim_dt(sim_length, h_vrai, nu)


Tmax = sim_length;
Tinf = 7; 
Times_obs=lapply(1:M,function(m){u_m=Times[[m]]; u_m = u_m[u_m[,1] > Tinf ,]; u_m=u_m[u_m[,1] <= Tmax ,]})
Times_utiles=lapply(1:M,function(m){u_m=Times[[m]]; u_m=u_m[u_m[,1] > Tinf-max(h_vrai$smax) ,];  u_m=u_m[u_m[,1] <= Tmax ,]})
obs_utiles =calc_absc_utiles_dt(Times_utiles,Times_obs,h_vrai$smax)
data=list(obs_utiles=obs_utiles);
data$Times_utiles = Times_utiles  
data$Times_obs = Times_obs
data$Tinf = rep(Tinf,M)
data$Tmax = rep(Tmax,M)
data$smax = h_vrai$smax 


#--------------------------------
#------  PRIOR DISTRIBUTION
#-------------------------------
p_Z = c(0.5,0.5); p_Z=p_Z/sum(p_Z)
pi_0 = 0.5;
#mu_lalpha =2.5; s_lalpha = 1
#mu_lnu = rep(3.5,M); s_lnu=rep(1,M)
#mu_lbeta =2.5; s_lbeta = 1
s_lalpha = 0.5; mu_lalpha =log(0.5)-s_lalpha^2
#s_lgamma = rep(0.5,M); mu_lgamma =rep(log(0.5)-s_lgamma^2,M)
s_lgamma = 1; mu_lgamma =log(1)-s_lgamma^2
s_lnu=0.5; mu_lnu = log(1)-s_lalpha^2;
s_lbeta = 0.5; mu_lbeta =log(1)-s_lalpha^2

hyperParam_prior=list(M=M,mu_lalpha=mu_lalpha,s_lalpha=s_lalpha,mu_lnu=mu_lnu,s_lnu=s_lnu,p_Z=p_Z,a_s=2,pi_0=pi_0,
                      mu_lbeta=mu_lbeta,s_lbeta=s_lbeta,mu_lgamma=mu_lgamma,s_lgamma=s_lgamma)
#hyperParam_prior$a_lambda_K = 20;
#hyperParam_prior$b_lambda_K = 1
#lambda_K = rgamma(1,hyperParam_prior$a_lambda_K,hyperParam_prior$b_lambda_K)
#hyperParam_prior$lambda_K=lambda_K




#-----------------------------
#--------- REVERSIBLE JUMP
#------------------------------

#------ tuning 
N_MCMC = 60000
burn_in = 30000
#N_MCMC = 10000
#burn_in = 5000
#N_MCMC = 800000
#burn_in = 400000
op_echan = list(alpha = c(1:M^2),K = c(1:M^2),nu = c(1:M),delta = c(),s = c(1:M^2), beta=c(1:M), gamma = c(1:M^2))
#op_echan = list(alpha = c(1:M^2),K = c(1:M^2),nu = c(1:M),lambda_K = c(1),delta = c(),s = c(1:M^2))
par_algo_MCMC =list(op_echan = op_echan)
par_algo_MCMC$op_affichage = 1000
par_algo_MCMC$N_MCMC = N_MCMC
par_algo_MCMC$rho_lnu = 1
par_algo_MCMC$rho_lalpha = 1
par_algo_MCMC$rho_lbeta = 1
par_algo_MCMC$rho_lgamma = 1




#--------- Initialisation

h_init = list()
h_init$smax = smax
h_init$delta = rep(1,M^2)
h_init$K =   rep(1,M^2)
h_init$s =   list(); for (p in  1:(M^2)) {h_init$s[[p]] =   c(0,h_init$smax[p])}
h_init$s_diffs =   list(); for (p in  1:(M^2)) {h_init$s_diffs[[p]] =   diff(h_init$s[[p]])}
h_init$gamma =   lapply(1:M^2,function(p){c(1)})
h_init$lgamma =   lapply(1:M^2,function(p){log(h_init$gamma[[p]])})
h_init$lbeta = rnorm(M^2, mu_lbeta, s_lbeta)
h_init$beta =   exp(h_init$lbeta)

#h_init$K =   rep(smax,M^2)
#h_init$s =   list(); for (p in  1:(M^2)) {h_init$s[[p]] =   seq(0,h_init$smax[p],length =   h_init$K[p] + 1)}
#h_init$s_diffs =   list(); for (p in  1:(M^2)) {h_init$s_diffs[[p]] =   diff(s[[p]])}
#h_init$alpha =   lapply(1:M^2,function(p){dnorm(s[[p]][-1],1.2,1)})
# h_init$alpha =   lapply(1:M^2,function(p){dnorm(s[[p]][-1],2,1)})
# h_init$lalpha =   lapply(1:M^2,function(p){log(h_init$alpha[[1]])})
# h_init$alpha =   lapply(1:M^2,function(p){rep(0,smax)})
# h_init$lalpha =   lapply(1:M^2,function(p){rep(-10000,smax)})
#h_init$gamma =   lapply(1:M^2,function(p){c(1,alpha[[p]][2]/alpha[[p]][1], alpha[[p]][3]/alpha[[p]][1])})
#h_init$gamma =   lapply(1:M^2,function(p){c(1,alpha[[p]][2]/alpha[[p]][1], alpha[[p]][3]/alpha[[p]][1], alpha[[p]][4]/alpha[[p]][1])})
# h_init$beta = vapply(1:M,function(m){0.9},1)
# h_init$lbeta = vapply(1:M,function(m){log(h_init$beta)},1)
#nu_init  =    vapply(1:M,function(m){length(data$Times_obs[[m]][,1])/(data$Tmax[m] - data$Tinf[m])},1)
#h_init$gamma =   lapply(1:M^2,function(p){c(1,1,1)})
# h_init$gamma =   lapply(1:M^2,function(p){c(1,1,1,1)})
# h_init$lgamma =   lapply(1:M^2,function(p){log(h_init$gamma[[1]])})
# h_init$lbeta = vapply(1:M,function(m){rnorm(1, mu_lbeta, s_lbeta)},1)
# h_init$beta =  vapply(1:M,function(m){exp(h_init$lbeta)},1)

nu_init  =  exp(rnorm(M, mu_lnu, s_lnu))
theta_init =   list(h =   h_init,nu =   nu_init)
#theta_init =   list(h =   h_init,nu =   nu_init,lambda_K =   lambda_K)
INPUT =   list(theta =   theta_init)

#------- Running MCMC (Reversible Jump) 
#resMCMC = RJ_Kernel_regular_s_dt(data,INPUT,hyperParam_prior,par_algo_MCMC)
resMCMC = RJ_Kernel_moving_s_dt(data,INPUT,hyperParam_prior,par_algo_MCMC)
resMCMC$INPUT <- INPUT
resMCMC$hyperParam_prior <- hyperParam_prior
resMCMC$par_algo_MCMC <- par_algo_MCMC


# library(purrr)
# # Calculate multivariate ESS
# for (i in 1:resMCMC$h[[1]]$smax){
#   resMCMCess = resMCMC
#   ind = which(sapply(resMCMCess$h, function(x) x$K == i))
#   nu_ind = resMCMC$nu[ind]
#   h_ind = lapply(ind, function(x) list(resMCMCess$h[[x]]))
#   h_ind_matrix = sapply(1:length(ind), function(x) c(h_ind[[1]][[x]]$beta, h_ind[[1]][[x]]$gamma[[1]], h_ind[[1]][[x]]$s[[1]], h_ind[[1]][[x]]$K))
# }


#----------------------------------------------------
#-------------------------- POSTERIOR inference
#---------------------------------------------------- 
part=seq(burn_in,par_algo_MCMC$N_MCMC,by=5)
L=length(part)

p=4
absc=seq(0,theta_vrai$h_vrai$smax[p],len=100)
traj = matrix(0,L,length(absc))
# for (j in 1:L){
#   h_func_j = stepfun(resMCMC$h[[part[j]]]$s[[1]],c(0,resMCMC$h[[part[j]]]$alpha[[1]],0),right=TRUE)
#   traj[j,]=h_func_j(absc)
# }
# h_func_vrai = stepfun(theta_vrai$h_vrai$s[[1]],c(0,theta_vrai$h_vrai$alpha[[1]],0),right=TRUE)
for (j in 1:L){
  x= resMCMC$h[[part[j]]]$s[[p]]
  y = c(0,resMCMC$h[[part[j]]]$gamma[[p]]/sum(resMCMC$h[[part[j]]]$gamma[[p]]),0)
  h_func_j = stepfun(x, y, right=TRUE)
  traj[j,]=h_func_j(absc)
}
h_func_vrai = stepfun(theta_vrai$h_vrai$s[[p]],c(0,theta_vrai$h_vrai$gamma[[p]]/sum(theta_vrai$h_vrai$gamma[[p]]),0),right=TRUE)
est_h = apply(traj,2,mean)
h_vrai = h_func_vrai(absc)
h_05 = vapply(1:100,function(j){quantile(traj[,j],0.05)},1)
h_50 = vapply(1:100,function(j){quantile(traj[,j],0.5)},1)
h_95 = vapply(1:100,function(j){quantile(traj[,j],0.95)},1)
H = c(h_vrai,est_h,h_50)
var.names = rep(c("True h","Post Mean","Median"),each = length(absc))
Vec.absc = rep(absc,3)
estim.h = data.frame(H = H,var.names,est_h,Time = Vec.absc)

#  
rib <- data.frame(Time = absc)
rib$H_min <- h_05;
rib$H_max <- h_95

estim.h$Func = factor(estim.h$var.names , levels = c("True h","Post Mean","Median"))


p <-  ggplot()
p <- p +  geom_ribbon(data = rib,aes(x = Time, ymin = H_min, ymax = H_max),alpha = 0.1,fill="blue",colour='black')
p <- p +  geom_line(data = estim.h, aes(x=Time,y=H,group = Func,colour = Func,linetype = Func),lwd=1.01)
p <- p + scale_fill_brewer(palette =  "Dark2") + scale_color_brewer(palette="Dark2")
p <- p + ylab("Intensities") + xlab('Time')
p 



# #----------------------------------------------
# # ------------- MCMC diagnostics
# #---------------------------------------------- 
# all_samples = t(sapply(1:resMCMC$par_algo_MCMC$N_MCMC, function(j) cbind(resMCMC$nu[j], t(resMCMC$h[[j]]$alpha[[1]]))))
# apply(all_samples,2,median)
# par(mfrow=c(3,1))
# #k=1
# for (k in 1:ncol(all_samples)){
#   plot(all_samples[,k], type = "l", main=paste0(par_names[k]))
#   plot(density(all_samples[,k]), type = "l")
#   #acf(all_samples[,k])
# }

# thinned_samples = t(sapply(1:L, function(j) cbind(resMCMC$nu[part[j]], t(resMCMC$h[[part[j]]]$alpha[[1]]) , t(resMCMC$h[[part[j]]]$beta))))
# par_names = c("nu", "alpha_1", "alpha_2", "alpha_3", "beta")
if (M==1){
  thinned_samples = t(sapply(1:L, function(j) cbind(resMCMC$nu[part[j]], t(resMCMC$h[[part[j]]]$beta), t(resMCMC$h[[part[j]]]$K))))
} else {
  thinned_samples = t(sapply(1:L, function(j) c(resMCMC$nu[part[j],], resMCMC$h[[part[j]]]$beta, resMCMC$h[[part[j]]]$K)))
}
#par_names = c("nu", "beta", "K")
#colnames(thinned_samples) = par_names
ind = which(sapply(1:ncol(thinned_samples), function(x) length(unique(thinned_samples[,x])))>1)
apply(thinned_samples,2,median)


par(mfrow=c(4,3))
# par(mfrow=c(3,1))
# ind=4
for (k in ind){
  #plot(thinned_samples[,k], type = "l", main=paste0(par_names[k]))
  plot(thinned_samples[,k], type = "l")
}
for (k in ind){
  plot(density(thinned_samples[,k]), type = "l")
}
for (k in ind){
  acf(thinned_samples[,k])
}

par(mfrow=c(1,1))
ggpairs(as.data.frame(thinned_samples[,ind]))



# Plot likelihood
h_vrai=list(alpha=alpha,s=s,smax=smax,K=K,delta=delta, beta=beta, gamma=gamma)
nu_vrai=nu;

par(mfrow=c(3,2))

nu_vec = seq(0,3,0.01)
log_lik = sapply(nu_vec, function(n){nu_ll= n; log_likelihood_multidim_dt(data, h_vrai, nu_ll)  })
plot(x=nu_vec, y=log_lik, main="nu")

# alpha1_vec = seq(0,1,0.01)
# log_lik = sapply(alpha1_vec, function(a1){h_ll= h_vrai; h_ll$alpha[[1]][1]=a1;h_ll$lalpha[[1]][1]=log(a1); log_likelihood_multidim_dt(data, h_ll, nu_vrai)  })
# plot(x=alpha1_vec, y=log_lik, main="alpha1")
#
# alpha2_vec = seq(0,1,0.01)
# log_lik = sapply(alpha2_vec, function(a2){h_ll= h_vrai; h_ll$alpha[[1]][2]=a2;h_ll$lalpha[[1]][2]=log(a2); log_likelihood_multidim_dt(data, h_ll, nu_vrai)  })
# plot(x=alpha2_vec, y=log_lik, main="alpha2")
#
# alpha3_vec = seq(0,0.5,0.01)
# log_lik = sapply(alpha3_vec, function(a3){h_ll= h_vrai; h_ll$alpha[[1]][3]=a3;h_ll$lalpha[[1]][3]=log(a3); log_likelihood_multidim_dt(data, h_ll, nu_vrai)  })
# plot(x=alpha3_vec, y=log_lik, main="alpha3")

beta_vec = seq(0,2,0.01)
log_lik = sapply(beta_vec, function(b){h_ll= h_vrai; h_ll$beta=b;h_ll$lbeta=log(b); log_likelihood_multidim_dt(data, h_ll, nu_vrai)  })
plot(x=beta_vec, y=log_lik, main="beta")

gamma1_vec = seq(0,5,0.01)
log_lik = sapply(gamma1_vec, function(g1){h_ll= h_vrai; h_ll$gamma[[1]][2]=g1;h_ll$lgamma[[1]][2]=log(g1); log_likelihood_multidim_dt(data, h_ll, nu_vrai)  })
plot(x=gamma1_vec, y=log_lik, main="gamma1")

gamma2_vec = seq(0,4,0.01)
log_lik = sapply(gamma2_vec, function(g2){h_ll= h_vrai; h_ll$gamma[[1]][3]=g2;h_ll$lgamma[[1]][3]=log(g2); log_likelihood_multidim_dt(data, h_ll, nu_vrai)  })
plot(x=gamma2_vec, y=log_lik, main="gamma2")

gamma3_vec = seq(0,1,0.01)
log_lik = sapply(gamma3_vec, function(g3){h_ll= h_vrai; h_ll$gamma[[1]][4]=g3;h_ll$lgamma[[1]][4]=log(g3); log_likelihood_multidim_dt(data, h_ll, nu_vrai)  })
plot(x=gamma3_vec, y=log_lik, main="gamma3")

par(mfrow=c(1,1))

