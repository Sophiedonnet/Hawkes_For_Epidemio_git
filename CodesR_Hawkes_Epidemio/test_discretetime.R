rm(list = ls())


library(parallel)
library(mvtnorm)
library(ggplot2)
library(GGally)

#source('functions/simulation_hawkes_multidim.R')
#source('functions/functions_likelihood_hawkes.R')
#source('functions/functions_likelihood_hawkes_positive_intensity.R')
#source('functions/RJ_Kernel_regular_s_positive_intensity.R')
#source('functions/RJ_Kernel_moving_s_positive_intensity.R')
source('functions/simulation_hawkes_discretetime.R')
source('functions/functions_likelihood_discretetime.R')
source('functions/RJ_Kernel_regular_discretetime.R')


op.plot= TRUE

#-----------------------------------
#-------- PARAM de simulation 
#-----------------------------------


M = 1
smax = rep(3,M^2)
nu = rep(1,M)

s = alpha = lalpha = list()
l = 1
m = 1
p = (m-1)*M+l
s[[p]] = seq(0,smax[p],1) 
alpha[[p]] <- dnorm(s[[p]][-1],1.2,1)
K = sapply(1:M^2,function(p){length(alpha[[p]])});
s = lapply(1:M^2,function(p){seq(0,smax[p],len=K[p]+1)})
delta = vapply(1:M^2,function(p){as.numeric(max(alpha[[p]])>0)},1)

h_vrai=list(alpha=alpha,s=s,smax=smax,K=K,delta=delta)
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

Times = simulate_data_multidim_dt(80, h_vrai, nu)


Tmax = 80;
Tinf = 7; 
Times_obs=lapply(1:M,function(m){u_m=Times[[m]]; u_m = u_m[u_m[,1]>Tinf ,]; u_m=u_m[u_m[,1] <= Tmax ,]})
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
mu_lalpha =2.5; s_lalpha = 1
mu_lnu = rep(3.5,M); s_lnu=rep(1,M)
hyperParam_prior=list(M=M,mu_lalpha=mu_lalpha,s_lalpha=s_lalpha,mu_lnu=mu_lnu,s_lnu=s_lnu,p_Z=p_Z,a_s=2,pi_0=pi_0)
#hyperParam_prior$a_lambda_K = 20;
#hyperParam_prior$b_lambda_K = 1
#lambda_K = rgamma(1,hyperParam_prior$a_lambda_K,hyperParam_prior$b_lambda_K)
#hyperParam_prior$lambda_K=lambda_K




#-----------------------------
#--------- REVERSIBLE JUMP
#------------------------------

#------ tuning 
N_MCMC = 30000
burn_in = 10000
#N_MCMC = 50000
op_echan = list(alpha = c(1:M^2),K = c(1:M^2),nu = c(1:M),delta = c(),s = c(1:M^2))
#op_echan = list(alpha = c(1:M^2),K = c(1:M^2),nu = c(1:M),lambda_K = c(1),delta = c(),s = c(1:M^2))
par_algo_MCMC =list(op_echan = op_echan)
par_algo_MCMC$op_affichage = 1000
par_algo_MCMC$N_MCMC = N_MCMC
par_algo_MCMC$rho_lnu = 1
par_algo_MCMC$rho_lalpha = 1


#--------- Initialisation

h_init = list()
h_init$smax = rep(smax,M^2)
h_init$delta = rep(1,M^2)
h_init$K =   rep(smax,M^2)
h_init$s =   list(); for (p in  1:(M^2)) {h_init$s[[p]] =   seq(0,h_init$smax[p],length =   h_init$K[p] + 1)}
h_init$alpha =   lapply(1:M^2,function(p){dnorm(s[[p]][-1],1.2,1)}) 
h_init$lalpha =   lapply(1:M^2,function(p){log(h_init$alpha[[1]])})
#h_init$alpha =   lapply(1:M^2,function(p){rep(0,smax)})
#h_init$lalpha =   lapply(1:M^2,function(p){rep(-10000,smax)})
nu_init  =    vapply(1:M,function(m){length(data$Times_obs[[m]][,1])/(data$Tmax[m] - data$Tinf[m])},1)
theta_init =   list(h =   h_init,nu =   nu_init)
#theta_init =   list(h =   h_init,nu =   nu_init,lambda_K =   lambda_K)
INPUT =   list(theta =   theta_init)


#------- Running MCMC (Reversible Jump) 
resMCMC = RJ_Kernel_regular_s_dt(data,INPUT,hyperParam_prior,par_algo_MCMC)
resMCMC$INPUT <- INPUT
resMCMC$hyperParam_prior <- hyperParam_prior
resMCMC$par_algo_MCMC <- par_algo_MCMC


#----------------------------------------------------
#-------------------------- POSTERIOR inference
#---------------------------------------------------- 
part=seq(burn_in,par_algo_MCMC$N_MCMC,by=5)
L=length(part)



absc=seq(0,theta_vrai$h_vrai$smax[1],len=100)
traj = matrix(0,L,length(absc))
for (j in 1:L){
  h_func_j = stepfun(resMCMC$h[[part[j]]]$s[[1]],c(0,resMCMC$h[[part[j]]]$alpha[[1]],0),right=TRUE)
  traj[j,]=h_func_j(absc)
}
h_func_vrai = stepfun(theta_vrai$h_vrai$s[[1]],c(0,theta_vrai$h_vrai$alpha[[1]],0),right=TRUE)
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

#----------------------------------------------
# ------------- MCMC diagnostics
#---------------------------------------------- 
all_samples = t(sapply(burn_in:resMCMC$par_algo_MCMC$N_MCMC, function(j) cbind(resMCMC$nu[j], t(resMCMC$h[[j]]$alpha[[1]]))))
colMeans(all_samples)
par(mfrow=c(3,1))
k=1
#for (k in 1:ncol(all_samples)){
  plot(all_samples[,k], type = "l", main=paste0(par_names[k]))
  plot(density(all_samples[,k]), type = "l")
  acf(all_samples[,k])
#}

thinned_samples = t(sapply(1:L, function(j) cbind(resMCMC$nu[part[j]], t(resMCMC$h[[part[j]]]$alpha[[1]]))))
par_names = c("nu", "alpha_1", "alpha_2", "alpha_3")
#, "alpha_4", "alpha_5", "alpha_6", "alpha_7"
colnames(thinned_samples) = par_names
ind = which(colSums(thinned_samples)>0)
colMeans(thinned_samples)

par(mfrow=c(3,1))
k=1
#for (k in ind){
  plot(thinned_samples[,k], type = "l", main=paste0(par_names[k]))
  plot(density(thinned_samples[,k]), type = "l")
  acf(thinned_samples[,k])
#}

ggpairs(as.data.frame(thinned_samples[,ind]))

