#library(mvtnorm)

RJ_Kernel_regular_s_dt = function(data, INPUT, hyperParam_prior, par_algo_MCMC,op.plot = FALSE, h_vrai = c(), op.test.on.rho = FALSE) { 
  #op_in_SMC=FALSE
  ## output  : theta,log_lik,acc_rate 
  ###
  M <- length(data$Times_obs)
  op_echan <- par_algo_MCMC$op_echan
  op_affichage <- par_algo_MCMC$op_affichage
  N_MCMC <- par_algo_MCMC$N_MCMC
  
  
  K_max <- 100
  pi_moves <- matrix(0.5,K_max,2); 
  pi_moves[1,] <- c(1,0)  ### proba naissance / mort 
  
  
  
  obs_utiles = data$obs_utiles
  Times_utiles = data$Times_utiles
  #Times_obs = data$Times_obs
  Tmax = data$Tmax; 
  if (length(Tmax) != M) { 
    Tmax = rep(Tmax,M)
  }
  
  Tinf = data$Tinf; 
  if (length(Tinf) != M) {
    Tinf = rep(Tinf,M)
  }
  
  #n_obs = vapply(1:M,function(m) { length(Times_obs[[m]])},1)
  smax = INPUT$theta$h$smax 
  
  #pi_0 = hyperParam_prior$pi_0
  p_Z =  hyperParam_prior$p_Z # poids de 0 dans la loi a priori des h
  mu_lalpha  = hyperParam_prior$mu_lalpha #
  s_lalpha = hyperParam_prior$s_lalpha
  mu_lnu  = hyperParam_prior$mu_lnu### nu : partie constante de l'intensité_
  s_lnu  = hyperParam_prior$s_lnu
  
  # a_lambda_K = hyperParam_prior$a_lambda_K
  # b_lambda_K = hyperParam_prior$b_lambda_K
  # a_s <- hyperParam_prior$a_s
  #pi_moves = prob_moves(K_min,K_max,p_K,grid_s=FALSE,op_echan_K=(length(op_echan$K)>0))
  
  
  #lambda_K = INPUT$theta$lambda_K
  
  
  
  
  
  
  
  ############################## INITIALISATION ########################## 
  nu = INPUT$theta$nu; 
  lnu = log(nu)
  h = INPUT$theta$h
  
  K = vapply(1:M^2,function(p) { length(h$alpha[[p]])},2)
  # delta = h$delta; 
  # s = h$s; 
  h$K = K; 
  
  
  
  
  ##########################################################################
  
  
  
  rho_lnu = par_algo_MCMC$rho_lnu
  rho_lalpha = par_algo_MCMC$rho_lalpha
  
  if (op.plot == TRUE) {
    par(mfrow = c(M,M))
  }
  #####################
  
  seq_MCMC <- list()
  seq_MCMC$nu = matrix(0,N_MCMC,M); 
  seq_MCMC$nu[1,] = nu  
  seq_MCMC$h = vector(mode = "list", length = N_MCMC); 
  seq_MCMC$h[[1]] <- h
  #seq_MCMC$lambda_K = rep(0,N_MCMC); 
  #seq_MCMC$lambda_K[1] <- lambda_K; 
  
  
  #idx <- func_idx(mat_obs_utiles,h$s)
  #J  <- calc_J(Times_utiles,h$s,Tinf,Tmax)
  
  log_lik = log_likelihood_multidim_dt(data,h,nu,op = 'vec')
  log_prior_h = dprior_h_regular_s(h,hyperParam_prior)
  log_prior_lnu = dnorm(lnu,mu_lnu,s_lnu,log = TRUE)

  # absc = seq(0,INPUT$theta$h$smax[1],len = 110)
  # list_mat_h = list(); 
  # for (p in 1:(M^2)) { 
  #   list_mat_h[[p]] = matrix(0,1000,length(absc))
  # }
  
  taux_accept_birth <-  taux_accept_death <- 0
  compt_accept_heights <-   0
  
  ########################################################
  ########################################################
  ########################################################
  for (it in 2:N_MCMC) { 
    
    if (it %% par_algo_MCMC$op_affichage == 0) { 
      #print(c(it,nu,lambda_K))
      print(c(it,nu))
    }
    
    
    ######################### nu | Times,aires
    
     for (m in op_echan$nu) { 
      rho_lnu = c(0.01,0.1,0.5)*par_algo_MCMC$rho_lnu; 
      S = sample(1:3,1)
      lnu_c = lnu;
      lnu_c[m] = lnu[m] + rnorm(1,0,rho_lnu[S]);
      nu_c = exp(lnu_c)
      
      log_lik_c  = log_likelihood_multidim_dt(data,h,nu_c,op = 'vec')
      log_prior_lnu_c_m = dnorm(lnu_c[m],mu_lnu[m],s_lnu[m],log = TRUE)
      p1 = log_lik_c[m] - log_lik[m]
      p2 = log_prior_lnu_c_m - log_prior_lnu[m]
      
      test_accept = log(runif(1)) <= (p1 + p2);
      
      if (test_accept == TRUE) { 
        nu = nu_c
        lnu = lnu_c
        log_lik = log_lik_c
        log_prior_lnu[m] = log_prior_lnu_c_m
      }
    }
    
    ######## Adjusted Langevin Metropolis on log nu
    
    if (length(op_echan$nu)  ==  M) {
      rho_lnu <- 0.2
      eval_lambda <- lambda_cond_multidim_dt(obs_utiles,h,nu)
      grad_log_post <- vapply(1:M,function(m) { nu[m]*sum(data$Times_obs[[m]][,2]/eval_lambda[[m]])},1) - nu[m]*(Tmax - Tinf) - 1/s_lnu^2*(lnu - mu_lnu) 
      lnu_c <- lnu  +  rho_lnu*grad_log_post  +  sqrt(2*rho_lnu)*rnorm(M); 
      nu_c <- exp(lnu_c)
      log_lik_c <- log_likelihood_multidim_dt(data,h,nu_c,op = "vec")
      log_prior_lnu_c <- dnorm(lnu_c,mu_lnu,s_lnu,log = TRUE)
      q_nu_c_nu <- sum(dnorm(lnu_c,lnu  +  rho_lnu * grad_log_post,2 * rho_lnu,log = TRUE))
      
      eval_lambda_c <- lambda_cond_multidim_dt(obs_utiles,h,nu_c)
      grad_log_post_c <- vapply(1:M,function(m) { nu_c[m]*sum(data$Times_obs[[m]][,2]/eval_lambda[[m]])},1) - nu_c[m]*(Tmax - Tinf) - 1/s_lnu^2*(lnu_c - mu_lnu) 
      
      q_nu_nu_c <- sum(dnorm(lnu,lnu_c  +  rho_lnu*grad_log_post_c,2 * rho_lnu,log = TRUE))
      proba_accept <- sum(log_lik_c - log_lik)  + sum(log_prior_lnu_c - log_prior_lnu) - (q_nu_c_nu - q_nu_nu_c);
      
      test_accept <- log(runif(1)) <= proba_accept;

      if (test_accept == TRUE) { 
        lnu <- lnu_c; nu <- exp(lnu)
        log_lik <- log_lik_c
        log_prior_lnu <- log_prior_lnu_c
      }
    }
    ###########################################
    #----------------------- move on lalpha
    #########################################
    
    for (p in op_echan$alpha) {
      
      m <- (p - 1) %/% M + 1
      
      for (k in 1:(h$K[p])) { 
        
        h_c <- h;          
        
        rho_lalpha <- c(0.01,0.1,1)*par_algo_MCMC$rho_lalpha;
        
        if (h$alpha[[p]][k] == 0) {
          h_c$lalpha[[p]][k] <- rnorm(1,mu_lalpha,s_lalpha); 
          p_Z_move = 1;
          
        }else{
          p_Z_move = 0.5
          Z = rbinom(1,1,p_Z_move); 
          
          if (Z == 0) { 
            h_c$lalpha[[p]][k] = -10000
          }else{
            h_c$lalpha[[p]][k] = rnorm(1,h$lalpha[[p]][k],rho_lalpha[sample(1:3,1)])
          }
        }
        h_c$alpha[[p]] = exp(h_c$lalpha[[p]])
        #------------ proba of move  



        if ((sum(h_c$alpha[[p]] == 0) == h_c$K[p])) { 
          proba_accept = -Inf 
        }else{
          p_Z_move_c = 0.5
          d1 = sum(vapply(1:3,function(j) { dnorm(h_c$lalpha[[p]][k],h$lalpha[[p]][k],rho_lalpha[j])},1))/3
          d1c = sum(vapply(1:3,function(j) { dnorm(h$lalpha[[p]][k],h_c$lalpha[[p]][k],rho_lalpha[j])},1))/3
          p3_1 = log((h$alpha[[p]][k] == 0) * dnorm(h_c$lalpha[[p]][k],mu_lalpha,s_lalpha)  +  
                       (h$alpha[[p]][k] > 0)*((1 - p_Z_move) * (h_c$alpha[[p]][k] == 0) +  
                        p_Z_move*(h_c$alpha[[p]][k] > 0)*d1))
          p3_2 = log((h_c$alpha[[p]][k] == 0) * dnorm(h$lalpha[[p]][k],mu_lalpha,s_lalpha)  +  
                       (h_c$alpha[[p]][k] > 0)*((1 - p_Z_move_c) * (h$alpha[[p]][k] == 0) +  
                        p_Z_move_c*(h$alpha[[p]][k] > 0)*d1c))
          log_p_passage = p3_1 - p3_2
          log_prior_h_c = dprior_h_regular_s(h_c,hyperParam_prior)
          log_lik_c = log_likelihood_multidim_dt(data,h_c,nu,op = 'vec')            
          proba_accept = sum(log_lik_c - log_lik) - log_p_passage  + log_prior_h_c$h - log_prior_h$h
        }
        test_accept = log(runif(1)) < proba_accept
        
        if (test_accept == TRUE) {  
          h = h_c; 
          log_lik = log_lik_c; 
          log_prior_h  = log_prior_h_c;  
          compt_accept_heights = compt_accept_heights + 1
        }
      }
    }
    

    #   ###########################################
    #----------------------- birth and death
    #########################################

    for (p in op_echan$K) {
      K_p = h$K[p];
      log_lik_c = log_lik
      m <- (p - 1) %/% M + 1
      h_c <- h;
      move = sample(c('birth','death'),1,prob = pi_moves[K_p,])
      if (move == 'birth') {
        indic_simu = 'birth'
        h_small <- h;
        h_big <- h;
        h_big$K[p] =  h_small$K[p]  +  1
        prob_k_star = diff(h_small$s[[p]])/h_small$smax[p]
        k_star = sample(1:(h_small$K[p]),1,prob = prob_k_star) #### on aura tendance à proposer de couper les longs paliers.
        h_big$s[[p]] = seq(0,smax[p], length.out = h_big$K[p] + 1)
        if (k_star > 1) { mu_k_star = 0.5*(h_small$lalpha[[p]][k_star - 1] + h_small$lalpha[[p]][k_star]); sigma = 0.1}
        if (k_star == 1) { mu_k_star = h_small$lalpha[[p]][1]; sigma = 0.1}
        lalpha_k_star = mu_k_star  +  sigma*rnorm(1);
        if (k_star == 1) { h_big$lalpha[[p]] = c(lalpha_k_star,h_small$lalpha[[p]])}
        if (k_star > 1) { h_big$lalpha[[p]] = c(h_small$lalpha[[p]][1:(k_star - 1)],lalpha_k_star,h_small$lalpha[[p]][(k_star):h_small$K[p]])}
        # Jacob = 1;
        log_prob_proposal = dnorm(lalpha_k_star,mu_k_star,sigma,log = TRUE) + log(1/h_small$K[p]) + log(pi_moves[h_small$K[p],1]) - log(pi_moves[h_big$K[p],2]);
        h_big$alpha = lapply(1:M^2,function(p) { exp(h_big$lalpha[[p]])})
        h_c = h_big;
      }

      if (move == 'death') {
        indic_simu = 'death'
        h_big <- h
        h_small <- h ;
        h_small$K[p] <- h_big$K[p] - 1
        k_star <- sample(1:(h_big$K[p] - 1),1)
        h_small$s[[p]] <- seq(0,smax[p],len = h_small$K[p] + 1);
        h_small$lalpha[[p]] <- h_big$lalpha[[p]][-k_star]
        lalpha_k_star <-  h_big$lalpha[[p]][k_star]
        if (k_star > 1) { mu_k_star = 0.5*(h_small$lalpha[[p]][k_star - 1] + h_small$lalpha[[p]][k_star]); sigma = 0.1};
        if (k_star == 1) { mu_k_star = h_small$lalpha[[p]][1]; sigma = 0.1}

        log_prob_proposal = dnorm(lalpha_k_star,mu_k_star,sigma,log = TRUE) + log(1 / h_small$K[p]) + log(pi_moves[h_small$K[p],1]) - log(pi_moves[h_big$K[p],2])
        log_prob_proposal = -log_prob_proposal
        h_small$alpha = lapply(1:M^2,function(p) { exp(h_small$lalpha[[p]])})
        h_c = h_small;
      }

      if (prod(h_c$alpha[[p]] == 0) == 1) {
        h_c$alpha[[p]] = c(0)
        h_c$lalpha[[p]] = c(-10000)
        h_c$delta[p] = 0;
        h_c$K[p] = 1
        h_c$s[[p]] = c(0,smax[p])
      }

      log_lik_c = log_lik
      log_prior_h_c = dprior_h_regular_s(h_c,hyperParam_prior)
      log_lik_c[m] = log_likelihood_multidim_dt(data,h_c,nu,m)
      proba_accept = log_lik_c[m] - log_lik[m]  + log_prior_h_c$h - log_prior_h$h - log_prob_proposal;
      test_accept = log(runif(1)) < proba_accept


      if ((test_accept == TRUE) & (move == "birth")) { taux_accept_birth = taux_accept_birth + 1 }
      if ((test_accept == TRUE) & (move == "death")) { taux_accept_death = taux_accept_death + 1 }  #
      if (test_accept == TRUE) {  h = h_c; log_lik = log_lik_c; log_prior_h  = log_prior_h_c}
    }

    ########################################### 
    ########### lambda_K
    ###########################################
    # indic_simu = 'lambda_K' 
    # if (op_echan$lambda_K == 1) {   
    #   lambda_K = rgamma(1,hyperParam_prior$a_lambda_K + sum(h$K[h$delta == 1] - 1),hyperParam_prior$b_lambda_K + sum(h$delta == 1))
    #   hyperParam_prior$lambda_K <- lambda_K
    # }
    # #   
    
    
    
    #   # ==  ==  ==  ==  ==  ==  == =   STOCKAGE des variables
    # #   
    # indic_simu = 'stock'   
    # for (p in 1:(M^2)) {  
    #   h_p = stepfun(h$s[[p]],c(0,h$alpha[[p]],0),right = T)
    #   list_mat_h[[p]][it %% 1000 + 1,] = h_p(absc)
    # }
    # 
    # indic_simu = 'plot'
    # if ((op.plot == TRUE) & (it %% op_affichage == 0)) { 
    #   par(mfrow = c(M,M))
    #   #plot(LL,type='l') 
    #   for (p in 1:M^2) { 
    #     absc = seq(0,smax[p],len = 110)
    #     l = (p - 1) %% M + 1
    #     m = (p - 1) %/% M + 1
    #     h_p = stepfun(h$s[[p]],c(0,h$alpha[[p]],0),right = T)
    #     curve(h_p,0,smax[p],col = 'red',lty = 2,lwd = 2,ylim = c(0,40),main = paste(l,'sur',m,'; iter',it,sep  = " "))
    #     lines(absc,colSums(list_mat_h[[p]])/1000,type = 's',col = 'magenta',lwd = 2)
    #     
    #   }
    # }
    
    seq_MCMC$nu[it,] = nu;
    seq_MCMC$h[[it]] = h;
    #seq_MCMC$lambda_K[[it]] =  lambda_K
    OUTPUT = seq_MCMC
    
    
  }
  return(OUTPUT)  }






