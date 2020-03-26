#library(mvtnorm)

RJ_Kernel_moving_s_positive_intensity=function(data,INPUT,hyperParam_prior,par_algo_MCMC,op.plot=FALSE,h_vrai=c(),op.lambda=FALSE){
#op.plot=TRUE
#op_in_SMC=FALSE
## output  : theta,log_lik,acc_rate
###
M <- length(data$Times_obs)
op_echan=par_algo_MCMC$op_echan
op_affichage = par_algo_MCMC$op_affichage
N_MCMC = par_algo_MCMC$N_MCMC

compt_accept_heights =truc = 0
LL=c()
K_min=1;
K_max=100
pi_moves=matrix(0.5,K_max,2); pi_moves[1,]=c(1,0)  ### proba naissance / mort 


### if op.lambda==TRUE on estime un lambda par fonction hkl  

mat_obs_utiles = data$mat_obs_utiles
Times_utiles = data$Times_utiles
Times_obs= data$Times_obs


Times_obs=data$Times_obs;
mat_obs_utiles = data$mat_obs_utiles
Times_utiles = data$Times_utiles
Tmax = data$Tmax; if (length(Tmax)!=M){Tmax=rep(Tmax,M)}
Tinf = data$Tinf ; if (length(Tinf)!=M){Tinf=rep(Tinf,M)}




n_obs=vapply(1:M,function(m){length(Times_obs[[m]])},1)

smax=INPUT$theta$h$smax

pi_0 = hyperParam_prior$pi_0
p_Z= hyperParam_prior$p_Z # poids de 0 dans la loi a priori des h
mu_lalpha  = hyperParam_prior$mu_lalpha #
s_lalpha = hyperParam_prior$s_lalpha
mu_lnu  = hyperParam_prior$mu_lnu### nu : partie constante de l'intensité_
s_lnu  = hyperParam_prior$s_lnu

a_lambda_K = hyperParam_prior$a_lambda_K
b_lambda_K = hyperParam_prior$b_lambda_K
a_s <- hyperParam_prior$a_s
#pi_moves = prob_moves(K_min,K_max,p_K,grid_s=FALSE,op_echan_K=(length(op_echan$K)>0))


lambda_K = INPUT$theta$lambda_K

if(op.lambda==TRUE){lambda_K = rep(lambda_K,M^2)}





############################## INITIALISATION ########################## 
nu= INPUT$theta$nu; lnu=log(nu)
h = INPUT$theta$h
K = vapply(1:M^2,function(p){length(h$alpha[[p]])},2)
delta=h$delta; 
s = h$s; 
h$K = K; 




##########################################################################



rho_lnu = par_algo_MCMC$rho_lnu
rho_lalpha = par_algo_MCMC$rho_lalpha
if(op.plot==TRUE){par(mfrow=c(M,M))}
#####################

seq_MCMC=list()
seq_MCMC$nu = matrix(0,N_MCMC,M); seq_MCMC$nu[1,] = nu  
seq_MCMC$h=vector(mode = "list", length = N_MCMC); seq_MCMC$h[[1]]=h
if(op.lambda==FALSE){seq_MCMC$lambda_K = rep(0,N_MCMC); seq_MCMC$lambda_K[1] <- lambda_K; }
if (op.lambda==TRUE){seq_MCMC$lambda_K = matrix(0,N_MCMC,M^2);seq_MCMC$lambda_K[1, ] <- lambda_K; }

  
idx <- func_idx(mat_obs_utiles,h$s)
J  <- calc_J(Times_utiles,h$s,Tinf,Tmax)

log_lik= log_likelihood_hawkes_multidim_positive_intensity(data,idx,J,h,nu,op='vec')
log_prior_h = dprior_h_moving_s(h,hyperParam_prior)
log_prior_lnu = dnorm(lnu,mu_lnu,s_lnu,log=TRUE)

  

absc=seq(0,INPUT$theta$h$smax[1],len=110)
absc=seq(0,INPUT$theta$h$smax[1],len=110)
list_mat_h=list(); for (p in 1:(M^2)){list_mat_h[[p]]= matrix(0,1000,length(absc))}




taux_accept_birth = taux_accept_death <-0


########################################################
########################################################
########################################################
  for (it in 2:N_MCMC){

  ## print(it)
  #if (it==2){print(c(it,nu,h$alpha[[1]]))}
    if (it%%par_algo_MCMC$op_affichage==0){print(c(it,nu,lambda_K))}
    
    
    ######################### nu | Times,aires

     indic_simu='nu'
     for (m in op_echan$nu){
       rho_lnu=c(0.01,0.1,0.5)*par_algo_MCMC$rho_lnu; S = sample(1:3,1)
       lnu_c= lnu;
       lnu_c[m] = lnu[m]+rnorm(1,0,rho_lnu[S]);
       nu_c = exp(lnu_c)

       log_lik_c  = log_likelihood_hawkes_multidim_positive_intensity(data,idx,J,h,nu_c,op='vec')
       #log_lik_c = vapply(1:M,function(m){length(Times_obs[[m]])*log(nu_c[m])-(Tmax[m]-Tinf[m])*nu_c[m]},1)

       log_prior_lnu_c_m = dnorm(lnu_c[m],mu_lnu[m],s_lnu[m],log=TRUE)
       p1 =log_lik_c[m] - log_lik[m]
       p2 =log_prior_lnu_c_m- log_prior_lnu[m]

       test_accept = log(runif(1))<=(p1+p2);

       if (test_accept==TRUE){
         nu = nu_c
         lnu = lnu_c
         log_lik =log_lik_c
         log_prior_lnu[m] =log_prior_lnu_c_m
       }
     }
  
    ######## Adjusted Langevin Metropolis on log nu

       if(length(op_echan$nu)==M){
         indic_simu <- 'nu_langevin'
          rho_lnu <-0.2
          eval_lambda <-lambda_cond_multidim_mat_positive_intensity(mat_obs_utiles,idx,h,nu)
         grad_log_post <- vapply(1:M,function(m){nu[m]*sum(1/eval_lambda[[m]])},1) -1/s_lnu^2*(lnu-mu_lnu)-nu*(Tmax-Tinf)
          lnu_c <- lnu + rho_lnu*grad_log_post + sqrt(2*rho_lnu)*rnorm(M); 
         nu_c=exp(lnu_c)
         log_lik_c <- log_likelihood_hawkes_multidim_positive_intensity(data,idx,J,h,nu_c,op="vec")
         log_prior_lnu_c <- dnorm(lnu_c,mu_lnu,s_lnu,log=TRUE)
         q_nu_c_nu <- sum(dnorm(lnu_c,lnu + rho_lnu*grad_log_post,2*rho_lnu,log=TRUE))
       
         eval_lambda_c<- lambda_cond_multidim_mat_positive_intensity(mat_obs_utiles,idx,h,nu_c)
         grad_log_post_c <- vapply(1:M,function(m){sum(nu_c[m]/eval_lambda_c[[m]])},1) -1/s_lnu^2*(lnu_c-mu_lnu)-nu_c*(Tmax-Tinf)
       
         q_nu_nu_c <-sum(dnorm(lnu,lnu_c + rho_lnu*grad_log_post_c,2*rho_lnu,log=TRUE))
         proba_accept<- sum(log_lik_c- log_lik) +sum(log_prior_lnu_c -log_prior_lnu) - (q_nu_c_nu -q_nu_nu_c);
       
         test_accept <- log(runif(1))<=proba_accept;
         if (test_accept==TRUE){
           
           lnu <- lnu_c; nu <-exp(lnu)
           log_lik <-log_lik_c
           log_prior_lnu <-log_prior_lnu_c
       }
    }
  #----------------------- move on delta
  indic_simu='delta'
  for (p in op_echan$delta){
    m <- (p-1)%/%M+1
    h_c <- h; 
      
    
    if(h$delta[p]==1){
      move='decrease'
      h_c$delta[p] <- 0; 
      h_c$alpha[[p]]<- c(0)
      h_c$lalpha[[p]]<- c(-10000)
      h_c$K[p] <- 1; 
      h_c$s[[p]] <- c(0,smax[p])
    }else{
        move='increase'
        h_c$delta[p] <- 1;
        if(op.lambda==TRUE){h_c$K[p] <- rpois(1,lambda_K[p])+1;}else{h_c$K[p] <- rpois(1,lambda_K)+1;}
        h_c$lalpha[[p]]<- rprior_lalpha(h_c$K[p],mu_lalpha,s_lalpha,p_Z)
        h_c$alpha[[p]] = exp(h_c$lalpha[[p]])
        h_c$s[[p]] <- rprior_s(h_c$K[p],smax[p],a_s)
      }
      
      idx_c <- func_idx(mat_obs_utiles,h_c$s)
      J_c  <- calc_J(Times_utiles,h_c$s,Tinf,Tmax)
      log_lik_c <- log_lik 
      log_lik_c[m] = log_likelihood_hawkes_multidim_positive_intensity_m(data,idx_c,J_c,h_c,nu,m)
      
      log_prior_h_c = dprior_h_moving_s(h_c,hyperParam_prior)
      if(move=='increase'){
        h_small=h; 
        h_big=h_c; 
        if(op.lambda==FALSE){
        q_small_to_big = dprior_s(h_big$s[[p]],hyperParam_prior$a_s,log=TRUE)+dpois(h_big$K[p]-1,lambda_K,log=TRUE)+dprior_lalpha(h_big$lalpha[[p]],mu_lalpha,s_lalpha,p_Z,log=TRUE)
        }
        if(op.lambda==TRUE){
          q_small_to_big = dprior_s(h_big$s[[p]],hyperParam_prior$a_s,log=TRUE)+dpois(h_big$K[p]-1,lambda_K[p],log=TRUE)+dprior_lalpha(h_big$lalpha[[p]],mu_lalpha,s_lalpha,p_Z,log=TRUE)
        }
        
        proba_accept = sum(log_lik_c-log_lik)+log_prior_h_c$h-log_prior_h$h-q_small_to_big
      }
      if(move=='decrease'){
        h_small=h_c; 
        h_big=h; 
        if(op.lambda==FALSE){
        q_small_to_big = dprior_s(h_big$s[[p]],hyperParam_prior$a_s,log=TRUE)+dpois(h_big$K[p]-1,lambda_K,log=TRUE)+dprior_lalpha(h_big$lalpha[[p]],mu_lalpha,s_lalpha,p_Z,log=TRUE)
        }else{
          q_small_to_big = dprior_s(h_big$s[[p]],hyperParam_prior$a_s,log=TRUE)+dpois(h_big$K[p]-1,lambda_K[p],log=TRUE)+dprior_lalpha(h_big$lalpha[[p]],mu_lalpha,s_lalpha,p_Z,log=TRUE)
          
          
        }
        
        proba_accept = sum(log_lik_c-log_lik)+log_prior_h_c$h-log_prior_h$h+q_small_to_big
        }
      
    
      test_accept = log(runif(1))<proba_accept
      if (test_accept==TRUE){ h = h_c; J = J_c; idx = idx_c; log_lik = log_lik_c; log_prior_h  = log_prior_h_c}
  }
#     
    
  
  #----------------------move heigths  ----------------
     indic_simu='heights'
     for (p in intersect(op_echan$alpha,which(h$delta==1))){
#       for (p in 2:2){
       m <- (p-1)%/%M+1
       for (k in 1:(h$K[p])){
       #  for (k in 1:1){
           
           
         h_c <- h;          
         rho_lalpha <- c(0.01,0.1,1)*par_algo_MCMC$rho_lalpha;
         if(h$alpha[[p]][k]==0){h_c$lalpha[[p]][k] <- rnorm(1,mu_lalpha,s_lalpha); p_Z_move = 1;
         }else{
           p_Z_move =0.5#1-exp(-h$alpha[[p]][k]/10) 
           Z = rbinom(1,1,p_Z_move); 
           if(Z==0){h_c$lalpha[[p]][k]=-10000}else{h_c$lalpha[[p]][k]=rnorm(1,h$lalpha[[p]][k],rho_lalpha[sample(1:3,1)])}
         }
         h_c$alpha[[p]]=exp(h_c$lalpha[[p]])
       #------------ proba of move  
         
         
        if((sum(h_c$alpha[[p]]==0)==h_c$K[p])){proba_accept=-Inf}else{ 
          p_Z_move_c = 0.5#1-exp(-h_c$alpha[[p]][k]/10) 
          d1 = sum(vapply(1:3,function(j){dnorm(h_c$lalpha[[p]][k],h$lalpha[[p]][k],rho_lalpha[j])},1))/3
          d1c = sum(vapply(1:3,function(j){dnorm(h$lalpha[[p]][k],h_c$lalpha[[p]][k],rho_lalpha[j])},1))/3
          p3_1 = log((h$alpha[[p]][k]==0)*dnorm(h_c$lalpha[[p]][k],mu_lalpha,s_lalpha) + (h$alpha[[p]][k]>0)*((1-p_Z_move) * (h_c$alpha[[p]][k]==0)+ p_Z_move*(h_c$alpha[[p]][k]>0)*d1))
          p3_2 = log((h_c$alpha[[p]][k]==0)*dnorm(h$lalpha[[p]][k],mu_lalpha,s_lalpha) + (h_c$alpha[[p]][k]>0)*((1-p_Z_move_c) * (h$alpha[[p]][k]==0)+ p_Z_move_c*(h$alpha[[p]][k]>0)*d1c))
          log_p_passage = p3_1-p3_2
        #log_p_passage=0
#       # ------------- mise ? jour 
#       
          log_prior_h_c = dprior_h_moving_s(h_c,hyperParam_prior)
          log_lik_c = log_likelihood_hawkes_multidim_positive_intensity(data,idx,J,h_c,nu,op='vec')            
          proba_accept = sum(log_lik_c-log_lik) -log_p_passage +log_prior_h_c$h-log_prior_h$h
        }
        
        
        test_accept = log(runif(1))<proba_accept
      
        if (test_accept==TRUE){ 
          h = h_c; 
          log_lik = log_lik_c; 
          log_prior_h  = log_prior_h_c;  compt_accept_heights = compt_accept_heights+1
        }
       }
      }
#print(compt_accept_heights/truc*100)
#   #----------------------move steps ----------------
  indic_simu = 's'
  for (p in intersect(intersect(op_echan$s,which(h$delta==1)),which(h$K>1))){
      
      m <- (p-1)%/%M+1
      for (k in 2:h$K[p]){
        #---- initialisation
        h_c = h; 
        s_c_p =s_p= h$s[[p]]
        #--- proposition
        rho_s = c(2,4); S=sample(c(1,2),1)
        beta_s = (s_p[k+1]-s_p[k-1])/(s_p[k]-s_p[k-1])*(rho_s[S]-1)-(rho_s[S]-2)
        s_c_p[k] = s_p[k-1] + (s_p[k+1]-s_p[k-1])*rbeta(1,rho_s[S],beta_s)
        h_c$s[[p]]=s_c_p
      
        #-----  probability of move
        beta_1s = (s_p[k+1]-s_p[k-1])/(s_p[k]-s_p[k-1])*(rho_s[1]-1)-(rho_s[1]-2)
        beta_2s = (s_p[k+1]-s_p[k-1])/(s_p[k]-s_p[k-1])*(rho_s[2]-1)-(rho_s[2]-2)
        z = (s_c_p[k]-s_p[k-1])/(s_p[k+1]-s_p[k-1])
        
        p3_sc_s = log(0.5/(s_p[k+1]-s_p[k-1])*(dbeta(z,rho_s[1],beta_1s)+dbeta(z,rho_s[2],beta_2s)))

        beta_1s_c = (s_c_p[k+1]-s_c_p[k-1])/(s_c_p[k]-s_c_p[k-1])*(rho_s[1]-1)-(rho_s[1]-2)
        beta_2s_c = (s_c_p[k+1]-s_c_p[k-1])/(s_c_p[k]-s_c_p[k-1])*(rho_s[2]-1)-(rho_s[2]-2)
        z_c = (s_p[k]-s_c_p[k-1])/(s_c_p[k+1]-s_c_p[k-1])
        p3_s_sc = log(0.5/(s_c_p[k+1]-s_c_p[k-1])*(dbeta(z_c,rho_s[1],beta_1s_c)+dbeta(z_c,rho_s[2],beta_2s_c)))
        #--------- updates
       
        idx_c <- func_idx(mat_obs_utiles,h_c$s)
        J_c  <- calc_J(Times_utiles,h_c$s,Tinf,Tmax)
        
        log_prior_h_c = dprior_h_moving_s(h_c,hyperParam_prior)
        log_lik_c = log_lik
        log_lik_c[m] = log_likelihood_hawkes_multidim_positive_intensity_m(data,idx_c,J_c,h_c,nu,m)            
        
        proba_accept = sum(log_lik_c-log_lik)+sum(log_prior_h_c$h-log_prior_h$h)-(p3_sc_s-p3_s_sc)

        test_accept = log(runif(1))<proba_accept
        
        if (test_accept==TRUE){ h = h_c; J = J_c; idx = idx_c; log_lik = log_lik_c; log_prior_h  = log_prior_h_c}
      }
    }
#     
#     
#   #-----------------------------------------
#   #----------------------birth and death  ----------------
#   
    for (p in intersect(op_echan$K,which(h$delta==1))){
       K_p=h$K[p]; 
       log_lik_c = log_lik
       m <- (p-1)%/%M+1
       h_c <- h; 
       move=sample(c('birth','death'),1,prob=pi_moves[K_p,])
       if(move=='birth'){
         indic_simu='birth'
         h_small <- h; 
         h_big <- h;  
         h_big$K[p] =  h_small$K[p] + 1
         prob_k_star = diff(h_small$s[[p]])/h_small$smax[p]
         k_star =sample(1:(h_small$K[p]),1,prob=prob_k_star) #### on aura tendance à proposer de couper les longs paliers. 
         h_big$s[[p]] = sort(c(h_small$s[[p]],runif(1,h_small$s[[p]][k_star],h_small$s[[p]][k_star+1]))) #### au final s* est simulée uniformement sur [0, smax]
         if(k_star>1){mu_k_star =0.5*(h_small$lalpha[[p]][k_star-1]+h_small$lalpha[[p]][k_star]); sigma = 0.1}
         if(k_star==1){mu_k_star =h_small$lalpha[[p]][1]; sigma = 0.1}
         lalpha_k_star = mu_k_star + sigma*rnorm(1);
         if(k_star==1){h_big$lalpha[[p]]=c(lalpha_k_star,h_small$lalpha[[p]])}
         if(k_star>1){h_big$lalpha[[p]]=c(h_small$lalpha[[p]][1:(k_star-1)],lalpha_k_star,h_small$lalpha[[p]][(k_star):h_small$K[p]])}
         Jacob = 1;  
         log_prob_proposal =dnorm(lalpha_k_star,mu_k_star,sigma,log=TRUE)+log(1/h_small$smax[p])+log(pi_moves[h_small$K[p],1])-log(pi_moves[h_big$K[p],2]);  
         h_big$alpha=lapply(1:M^2,function(p){exp(h_big$lalpha[[p]])})
         h_c=h_big; 
       }
   
       if (move=='death'){
         indic_simu='death'
         h_big <- h 
         h_small <- h ;
         h_small$K[p] <- h_big$K[p]-1 
         k_star <- sample(1:(h_big$K[p]-1),1) 
         h_small$s[[p]] <- h_big$s[[p]][-(k_star+1)]; 
         h_small$lalpha[[p]] <- h_big$lalpha[[p]][-k_star]
         lalpha_k_star <-  h_big$lalpha[[p]][k_star]
         if(k_star>1){mu_k_star =0.5*(h_small$lalpha[[p]][k_star-1]+h_small$lalpha[[p]][k_star]); sigma =0.1};# max(abs(h_small$lalpha[[p]][k_star-1]-h_small$lalpha[[p]][k_star]),0.1)}
         if(k_star==1){mu_k_star =h_small$lalpha[[p]][1]; sigma = 0.1}#max(h_small$lalpha[[p]][k_star],0.1)} 
         
         log_prob_proposal =dnorm(lalpha_k_star,mu_k_star,sigma,log=TRUE)+log(1/h_small$smax[p])+log(pi_moves[h_small$K[p],1])-log(pi_moves[h_big$K[p],2]);  
         log_prob_proposal = -log_prob_proposal
         h_small$alpha=lapply(1:M^2,function(p){exp(h_small$lalpha[[p]])})
         h_c=h_small;
        }
       
        if(prod(h_c$alpha[[p]]==0)==1){
          h_c$alpha[[p]]=c(0)
          h_c$lalpha[[p]]=c(-10000)
          h_c$delta[p]=0; 
          h_c$K[p]=1
          h_c$s[[p]]=c(0,smax[p])
        }
        idx_c <- func_idx(mat_obs_utiles,h_c$s)
        J_c  <- calc_J(Times_utiles,h_c$s,Tinf,Tmax)
        log_lik_c= log_lik
        log_prior_h_c = dprior_h_moving_s(h_c,hyperParam_prior)
        log_lik_c[m] = log_likelihood_hawkes_multidim_positive_intensity_m(data,idx_c,J_c,h_c,nu,m)             
        proba_accept = log_lik_c[m]-log_lik[m] +log_prior_h_c$h-log_prior_h$h-log_prob_proposal; 
        test_accept = log(runif(1))<proba_accept
        if((test_accept==TRUE)&(move=="birth")){taux_accept_birth = taux_accept_birth+1 }
        if((test_accept==TRUE)&(move=="death")){taux_accept_death = taux_accept_death+1 }  #       
        if (test_accept==TRUE){ h = h_c; J = J_c; idx = idx_c; log_lik = log_lik_c; log_prior_h  = log_prior_h_c}
       }
#       
    ########### SPECIAL MOVE : passer de 2K ? K  
#     indic_simu = "Special move K/2"
#     for (p in intersect(which(h$K>2),intersect(op_echan$K,which(h$delta==1)))){
#       m <- (p-1)%/%M+1
#       
#       if (h$K[p]%%2==0){
#         h_c <- h; 
#         h_c$K[p] <- h$K[p]/2;
#         
#         h_c$s[[p]] <- h$s[[p]][seq(1,h$K[p]+1,by=2)]
#         mid <- (h$alpha[[p]][seq(2,h$K[p],by=2)] +h$alpha[[p]][seq(1,h$K[p],by=2)])*0.5
#         h_c$alpha[[p]] <- mid; 
#         u <- 2*(h$alpha[[p]][seq(2,h$K[p],by=2)]-mid)
#         
#         #%% likelihood
#         idx_c <- func_idx(mat_obs_utiles,h_c$s)
#         J_c  <- calc_J(Times_utiles,h_c$s,Tinf,Tmax)
#         log_lik_c <- log_lik
#         log_lik_c[m] = log_likelihood_hawkes_multidim_positive_intensity_m(data,idx_c,J_c,h_c,nu,m)
#         
#         log_prior_h_c=dprior_h_moving_s(h_c,hyperParam_prior)
#         q_c = sum(dnorm(u,0,10,log=TRUE))-sum(log(diff(h$s[[p]])))
#         
#         proba_accept = sum(log_lik_c-log_lik)+log_prior_h_c$h-log_prior_h$h+ sum(q_c)
#         test_accept = log(runif(1))<proba_accept
#         if (test_accept==TRUE){
#           #print(c(p,'accept special move K/2'))
#           h = h_c
#           log_lik = log_lik_c
#           log_prior_h= log_prior_h_c
#           J = J_c
#           idx=idx_c
#           }
#       }
#     }
#     
#     indic_simu = "Special move *2"
#     for (p in intersect(op_echan$K,which(h$delta==1))){
#       m <- (p-1)%/%M+1
#       log_lik_c = log_lik
#       h_c <- h; 
#       h_c$K[p] <- h$K[p]*absc=seq(0,INPUT$theta$h$smax[1],len=110)
     #list_mat_h=list()
     #if(op.plot==TRUE){for (p in 1:(M^2)){list_mat_h[[p]]= matrix(0,1000,length(absc))}}2; 
#       h_c$s[[p]] <- sort(c(h$s[[p]],vapply(1:(h$K[p]),function(k){runif(1,h$s[[p]][k],h$s[[p]][k+1])},1)))
#       
#       s_dep =h$alpha[[p]]/2
#       u <- s_dep*rnorm(h$K[p],0,1)
#       
#       h_c$alpha[[p]] <-c(vapply(1:(h$K[p]),function(k){c(h$alpha[[p]][k]-u[k]/2,h$alpha[[p]][k]+u[k]/2)},c(1,1)))
#       #%% likelihood
#       if(min(h_c$alpha[[p]])<0){test_accept<-FALSE}else{
#         idx_c <- func_idx(mat_obs_utiles,h_c$s)
#         J_c  <- calc_J(Times_utiles,h_c$s,Tinf,Tmax)
#         log_lik_c[m] = log_likelihood_hawkes_multidim_positive_intensity_m(data,idx_c,J_c,h_c,nu,m)
#         log_prior_h_c =dprior_h_moving_s(h_c,hyperParam_prior)
#         q_c = sum(dnorm(u[s_dep>0],0,s_dep[s_dep>0],log=TRUE)) +sum(dunif(h_c$s[[p]][seq(2,h_c$K[p],by=2)],h$s[[p]][1:(h$K[p])],h$s[[p]][2:(h$K[p]+1)],log=TRUE )); 
#         
#         proba_accept = sum(log_lik_c-log_lik)+log_prior_h_c$h-log_prior_h$h - q_c
#         test_accept = log(runif(1))<proba_accept
#       }
#       if (test_accept==TRUE){
#         h = h_c
#         log_lik = log_lik_c
#         log_prior_h= log_prior_h_c
#         J = J_c
#         idx=idx_c
#        }
#     }
#     
#     ################################### 
#     
#     indic_simu = 'Special move 1 -> K'
#     
#     for (p in intersect(intersect(op_echan$K,which(h$delta==1)),which(h$K==1))){
#       m <- (p-1)%/%M+1
#       K = 1
#       log_lik_c = log_lik
#       h_c <- h; 
#       h_c$K[p] <- rpois(1,2*lambda_K)+2; 
#       h_c$s[[p]] <- rprior_s(h_c$K[p],smax[p],a_s)
#       
#       d <- rdirichlet(1,rep(4,h_c$K[p]))
#       h_c$alpha[[p]] <- 1/K*d*h$alpha[[p]]
#       
#       
#       #%% likelihood
#       idx_c <- func_idx(mat_obs_utiles,h_c$s)
#       J_c  <- calc_J(Times_utiles,h_c$s,Tinf,Tmax)
#       log_lik_c[m] = log_likelihood_hawkes_multidim_positive_intensity_m(data,idx_c,J_c,h_c,nu,m)
#       log_prior_h_c = dprior_h_moving_s(h_c,hyperParam_prior)
#       
#       q_c = ddirichlet(d,rep(4,h_c$K[p]),log=TRUE) + dpois(h_c$K[p]-2,2*lambda_K,log=TRUE) + dprior_s(h_c$s[[p]],a_s,log=TRUE)
#       log_Jac = h_c$K[p]*log(h_c$K[p])+(h_c$K[p]-1)*log(h$alpha[[p]])
#       
#       
#       proba_accept = sum(log_lik_c-log_lik)+log_prior_h_c$h-log_prior_h$h  - sum(q_c)-log_Jac
#       test_accept = log(runif(1))<proba_accept
#       if(min(h_c$alpha[[p]])<0){test_accept<-FALSE}
#       if (test_accept==TRUE){
#         #print(c(p,'accept special move +'))
#         h = h_c
#         log_lik = log_lik_c
#         log_prior_h_c = log_prior_h
#         J = J_c
#         idx=idx_c
#       }
#     }
#     
#     
#         
    indic_simu='lambda_K' 
    if(op_echan$lambda_K==1){  
      if(op.lambda==FALSE){
       lambda_K=rgamma(1,hyperParam_prior$a_lambda_K+sum(h$K[h$delta==1]-1),hyperParam_prior$b_lambda_K+sum(h$delta==1))}
      if(op.lambda==TRUE){
         lambda_K=rgamma(M^2,hyperParam_prior$a_lambda_K+h$K-1,hyperParam_prior$b_lambda_K+1)}
         
       hyperParam_prior$lambda_K <- lambda_K
     }
#   
  
  
  
#   #===============   STOCKAGE des variables
#   
  indic_simu='stock'   
  for (p in 1:(M^2)){ 
  h_p = stepfun(h$s[[p]],c(0,h$alpha[[p]],0),right=T)
  list_mat_h[[p]][it%%1000+1,] =h_p(absc)
  }
  
  indic_simu='plot'
  if ((op.plot==TRUE)&(it%%op_affichage==0)){
      par(mfrow=c(M,M))
      #plot(LL,type='l') 
      for(p in 1:M^2){
        absc=seq(0,smax[p],len=110)
        l = (p-1)%%M+1
        m =(p-1)%/%M+1
        h_p = stepfun(h$s[[p]],c(0,h$alpha[[p]],0),right=T)
        #h_c_p = stepfun(h_c$s[[p]],c(0,h_c$alpha[[p]],0),right=T)
        curve(h_p,0,smax[p],col='red',lty=2,lwd=2,ylim=c(0,40),main=paste(l,'sur',m,'; iter',it,sep=" "))
        if(length(h_vrai)>0){
          h_vrai_p = stepfun(h_vrai$s[[p]],c(0,h_vrai$alpha[[p]],0),right=T)
          curve(h_vrai_p,0,smax[p],col='orange',lwd=2,add=TRUE)
        }
        #curve(h_c_p,0,smax[p],col='blue',lwd=2,add=TRUE)
        lines(absc,colSums(list_mat_h[[p]])/1000,type='s',col='magenta',lwd=2)
        
      }
  }
    
  seq_MCMC$nu[it,] = nu;
  seq_MCMC$h[[it]] = h;
  if(op.lambda==FALSE){seq_MCMC$lambda_K[it] =  lambda_K}else{seq_MCMC$lambda_K[it,]=lambda_K} 
  
  OUTPUT = seq_MCMC
  
  }
return(OUTPUT)  }






