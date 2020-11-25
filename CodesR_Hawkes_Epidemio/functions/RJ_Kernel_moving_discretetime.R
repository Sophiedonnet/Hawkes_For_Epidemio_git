#library(mvtnorm)

RJ_Kernel_moving_s_dt = function(data, INPUT, hyperParam_prior, par_algo_MCMC,op.plot = FALSE, h_vrai = c(), op.test.on.rho = FALSE) { 
  #op_in_SMC=FALSE
  ## output  : theta,log_lik,acc_rate 
  ###
  M <- length(data$Times_obs)
  op_echan <- par_algo_MCMC$op_echan
  op_affichage <- par_algo_MCMC$op_affichage
  N_MCMC <- par_algo_MCMC$N_MCMC
  
  
  
  
  obs_utiles = data$obs_utiles
  #Times_utiles = data$Times_utiles
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
  
  pi_moves <- matrix(0.5,smax[1],2); 
  pi_moves[1,] <- c(1,0)  ### proba naissance / mort 
  pi_moves[smax[1],] <- c(0,1)  ### proba naissance / mort 
  
  
  
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
  
  K = vapply(1:M^2,function(p) { length(h$gamma[[p]])},2)
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
  #log_prior_h = dprior_h_regular_s(h,hyperParam_prior)
  log_prior_lgamma = vapply(1:M^2, function(i) sum(dnorm(h$lgamma[[i]],mu_lgamma,s_lgamma, log=TRUE)), 1)
  log_prior_lnu = vapply(1:M, function(i) dnorm(lnu[i],mu_lnu,s_lnu,log = TRUE), 1)
  log_prior_lbeta = vapply(1:M^2, function(i) dnorm(h$lbeta[i],mu_lbeta,s_lbeta,log = TRUE), 1)
  log_prior_s = (factorial(h$K)*factorial(smax - h$K))/(factorial(smax))
  
  # absc = seq(0,INPUT$theta$h$smax[1],len = 110)
  # list_mat_h = list(); 
  # for (p in 1:(M^2)) { 
  #   list_mat_h[[p]] = matrix(0,1000,length(absc))
  # }
  
  taux_accept_birth <-  taux_accept_death <- 0
  taux_accept_split <-  taux_accept_merge <- 0
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
      log_prior_lnu_c_m = dnorm(lnu_c[m],mu_lnu,s_lnu,log = TRUE)
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
    
    # ######## Adjusted Langevin Metropolis on log nu
    # 
    # if (length(op_echan$nu)  ==  M) {
    #   rho_lnu <- 0.1
    #   eval_lambda <- lambda_cond_multidim_dt(obs_utiles,h,nu)
    #   grad_log_post <- vapply(1:M,function(m) { nu[m]*sum(data$Times_obs[[m]][,2]/eval_lambda[[m]])},1) - nu[m]*(Tmax - Tinf) - 1/s_lnu^2*(lnu - mu_lnu)
    #   lnu_c <- lnu  +  rho_lnu*grad_log_post  +  sqrt(2*rho_lnu)*rnorm(M);
    #   nu_c <- exp(lnu_c)
    #   log_lik_c <- log_likelihood_multidim_dt(data,h,nu_c,op = "vec")
    #   log_prior_lnu_c <- dnorm(lnu_c,mu_lnu,s_lnu,log = TRUE)
    #   q_nu_c_nu <- sum(dnorm(lnu_c,lnu  +  rho_lnu * grad_log_post,2 * rho_lnu,log = TRUE))
    # 
    #   eval_lambda_c <- lambda_cond_multidim_dt(obs_utiles,h,nu_c)
    #   grad_log_post_c <- vapply(1:M,function(m) { nu_c[m]*sum(data$Times_obs[[m]][,2]/eval_lambda[[m]])},1) - nu_c[m]*(Tmax - Tinf) - 1/s_lnu^2*(lnu_c - mu_lnu)
    # 
    #   q_nu_nu_c <- sum(dnorm(lnu,lnu_c  +  rho_lnu*grad_log_post_c,2 * rho_lnu,log = TRUE))
    #   proba_accept <- sum(log_lik_c - log_lik)  + sum(log_prior_lnu_c - log_prior_lnu) - (q_nu_c_nu - q_nu_nu_c);
    # 
    #   test_accept <- log(runif(1)) <= proba_accept;
    # 
    #   if (test_accept == TRUE) {
    #     lnu <- lnu_c; nu <- exp(lnu)
    #     log_lik <- log_lik_c
    #     log_prior_lnu <- log_prior_lnu_c
    #   }
    # }
    # 
    
    ######################### beta

    for (m in op_echan$beta) {
      rho_lbeta = c(0.01,0.1,0.5)*par_algo_MCMC$rho_lbeta;
      S = sample(1:3,1)
      h_c = h;
      h_c$lbeta[m] = h$lbeta[m] + rnorm(1,0,rho_lbeta[S]);
      h_c$beta = exp(h_c$lbeta)

      log_lik_c  = log_likelihood_multidim_dt(data,h_c,nu,op = 'vec')
      log_prior_lbeta_c_m = dnorm(h_c$lbeta[m],mu_lbeta,s_lbeta,log = TRUE)
      p1 = log_lik_c[m] - log_lik[m]
      p2 = log_prior_lbeta_c_m - log_prior_lbeta[m]

      test_accept = log(runif(1)) <= (p1 + p2);

      if (test_accept == TRUE) {
        h = h_c
        log_lik = log_lik_c
        log_prior_lbeta[m] = log_prior_lbeta_c_m
      }
    }

    # ###########################################
    # #----------------------- move on lgamma
    # #########################################

    for (p in intersect(op_echan$gamma, which(h$K > 1))) {

      m <- (p - 1) %/% M + 1

      for (k in 2:(h$K[p])){
        #for (k in c(3)){
        
        h_c <- h;
        rho_lgamma <- c(0.01,0.1,0.5)*par_algo_MCMC$rho_lgamma;
        h_c$lgamma[[p]][k] = rnorm(1,h$lgamma[[p]][k],rho_lgamma[sample(1:3,1)])
        h_c$gamma[[p]][k] = exp(h_c$lgamma[[p]][k])

        log_lik_c = log_likelihood_multidim_dt(data,h_c,nu,op = 'vec')
        log_prior_lgamma_c = sum(dnorm(h_c$lgamma[[p]],mu_lgamma,s_lgamma,log = TRUE))

        proba_accept = sum(log_lik_c - log_lik) + log_prior_lgamma_c - log_prior_lgamma[p]
        #proba_accept = sum(log_lik_c - log_lik)

        test_accept = log(runif(1)) < proba_accept
        
        # print(proba_accept)
        # print(test_accept)
        # print(it)
        # print(h_c$gamma[[1]])
        # print(log_lik_c)
        # print(log_lik)
        # print(log_prior_lgamma_c)
        # print(log_prior_lgamma)

        if (test_accept == TRUE) {
          h = h_c;
          log_lik = log_lik_c;
          log_prior_lgamma[p] = log_prior_lgamma_c
          compt_accept_heights = compt_accept_heights + 1
        }
      }
    }
    
   
    # ###########################################
    # #----------------------- move on s
    # #########################################
    
    for (p in intersect(op_echan$s, which(h$K > 1))) {
      
      m <- (p - 1) %/% M + 1
      if (length(2:(h$K[p])) ==1){
        rand_s = 2
      } else {
        rand_s = sample(2:(h$K[p]), replace=FALSE)
      }
      for (k in rand_s){
        
        h_c <- h;
        
        # Determine possible moves
        vacant_s = setdiff(seq(h$s[[p]][(k-1)],h$s[[p]][(k+1)]), c(h$s[[p]][(k-1)], h$s[[p]][k], h$s[[p]][(k+1)]))
        
        if (length(vacant_s)>0){
          
          # Propose move
          if (length(vacant_s) ==1){
            h_c$s[[p]][k] = as.numeric(vacant_s)
          } else {
            h_c$s[[p]][k] = sample(vacant_s,1)
          }
          h_c$s_diffs = lapply(1:M^2,function(p) { diff(h_c$s[[p]])})
          
          log_lik_c = log_likelihood_multidim_dt(data,h_c,nu,op = 'vec')
            
          proba_accept = sum(log_lik_c - log_lik) # prior not dependent on s, only K which is constant here
          
          test_accept = log(runif(1)) < proba_accept

          if (test_accept == TRUE) {
            h = h_c;
            log_lik = log_lik_c;
          }
        }
      }
    }
    
    
    
    ###########################################
    #----------------------- move on lalpha
    #########################################

    # for (p in op_echan$alpha) {
    # 
    #   m <- (p - 1) %/% M + 1
    # 
    #   for (k in 1:(h$K[p])) {
    #   #for (k in c(2)) {
    # 
    #     h_c <- h;
    # 
    #     # if (k!=3){
    #        rho_lalpha <- c(0.01,0.1,0.5)*par_algo_MCMC$rho_lalpha;
    #     # } else {
    #     #  rho_lalpha <- c(0.01,0.05,0.1)*par_algo_MCMC$rho_lalpha;
    #     #}
    # 
    #     # if (h$alpha[[p]][k] == 0) {
    #     #   h_c$lalpha[[p]][k] <- rnorm(1,mu_lalpha,s_lalpha);
    #     #   p_Z_move = 1;
    #     #
    #     # }else{
    #       #p_Z_move = 0.5
    #       # p_Z_move = 1
    #       # Z = rbinom(1,1,p_Z_move);
    # 
    #       # if (Z == 0) {
    #       #   h_c$lalpha[[p]][k] = -10000
    #       # }else{
    #         h_c$lalpha[[p]][k] = rnorm(1,h$lalpha[[p]][k],rho_lalpha[sample(1:3,1)])
    #     #   }
    #     # }
    #     h_c$alpha[[p]] = exp(h_c$lalpha[[p]])
    #     #------------ proba of move
    # 
    # 
    # 
    #     if ((sum(h_c$alpha[[p]] == 0) == h_c$K[p])) {
    #       proba_accept = -Inf
    #     }else{
    #       # p_Z_move_c = 0.5
    #       # d1 = sum(vapply(1:3,function(j) { dnorm(h_c$lalpha[[p]][k],h$lalpha[[p]][k],rho_lalpha[j])},1))/3
    #       # d1c = sum(vapply(1:3,function(j) { dnorm(h$lalpha[[p]][k],h_c$lalpha[[p]][k],rho_lalpha[j])},1))/3
    #       # p3_1 = log((h$alpha[[p]][k] == 0) * dnorm(h_c$lalpha[[p]][k],mu_lalpha,s_lalpha)  +
    #       #              (h$alpha[[p]][k] > 0)*((1 - p_Z_move) * (h_c$alpha[[p]][k] == 0) +
    #       #                                       p_Z_move*(h_c$alpha[[p]][k] > 0)*d1))
    #       # p3_2 = log((h_c$alpha[[p]][k] == 0) * dnorm(h$lalpha[[p]][k],mu_lalpha,s_lalpha)  +
    #       #              (h_c$alpha[[p]][k] > 0)*((1 - p_Z_move_c) * (h$alpha[[p]][k] == 0) +
    #       #                                         p_Z_move_c*(h$alpha[[p]][k] > 0)*d1c))
    #       #
    #       #
    #       # log_p_passage = p3_1 - p3_2
    #       # log_prior_h_c = dprior_h_regular_s(h_c,hyperParam_prior)
    #       log_p_passage=0
    #       log_prior_h_c = sum(dnorm(h_c$lalpha[[p]],mu_lalpha,s_lalpha, log=TRUE))
    #       log_lik_c = log_likelihood_multidim_dt(data,h_c,nu,op = 'vec')
    # 
    #       #proba_accept = sum(log_lik_c - log_lik) - log_p_passage  + log_prior_h_c$h - log_prior_h$h
    #       proba_accept = sum(log_lik_c - log_lik) - log_p_passage  + log_prior_h_c - log_prior_h
    #     }
    #     test_accept = log(runif(1)) < proba_accept
    # 
    #     if (test_accept == TRUE) {
    #       h = h_c;
    #       log_lik = log_lik_c;
    #       log_prior_h  = log_prior_h_c;
    #       compt_accept_heights = compt_accept_heights + 1
    #     }
    #   }
    # }


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
        
        vacant_s = setdiff(seq(0,smax[p]), h_small$s[[p]])
        if (length(vacant_s)>0){
          
          if (length(vacant_s) ==1){
            k_star_val = as.numeric(vacant_s)
          } else {
            k_star_val = as.numeric(sample(vacant_s,1))
          }
          
          h_big$s[[p]] = sort(c(h_small$s[[p]], k_star_val))
          k_star = which(h_big$s[[p]]==k_star_val)-1
          
          # if (k_star < h_small$K[p]) { 
          #   mu_k_star = 0.5*(h_small$lgamma[[p]][k_star] + h_small$lgamma[[p]][k_star+1])
          #   sigma = 0.1
          # }
          # if (k_star == h_small$K[p]) { 
          #   mu_k_star = h_small$lgamma[[p]][h_small$K[p]]
          #   sigma = 0.1
          # }
          if (h_small$K[p] > 1){
            mu_k_star = mean(h_small$lgamma[[p]][-1])
          } else {
            mu_k_star = 0
          }
          sigma = 0.1
          
          lgamma_k_star = mu_k_star + sigma*rnorm(1);
          if (k_star < h_small$K[p]) { 
            h_big$lgamma[[p]] = c(h_small$lgamma[[p]][1:k_star],lgamma_k_star,h_small$lgamma[[p]][(k_star+1):h_small$K[p]])
          }
          if (k_star == h_small$K[p]) { 
            h_big$lgamma[[p]] = c(h_small$lgamma[[p]], lgamma_k_star)
          }
          
          # Jacob = 1;
          log_prob_proposal = dnorm(lgamma_k_star,mu_k_star,sigma,log = TRUE) + log(1/(smax[p]-h_small$K[p])) + log(pi_moves[h_small$K[p],1]) - 
                              log(pi_moves[h_big$K[p],2]) - log(1/h_small$K[p]);
          h_big$gamma = lapply(1:M^2,function(p) { exp(h_big$lgamma[[p]])})
          h_big$s_diffs = lapply(1:M^2,function(p) { diff(h_big$s[[p]])})
          h_c = h_big;
        }
      }

      if (move == 'death') {
        indic_simu = 'death'
        h_big <- h
        h_small <- h ;
        h_small$K[p] <- h_big$K[p] - 1
        
        k_star <- as.numeric(sample(1:(h_small$K[p]),1))
        h_small$s[[p]] <- h_big$s[[p]][-(k_star + 1)]
        
        h_small$lgamma[[p]] <- h_big$lgamma[[p]][-(k_star+1)]
        lgamma_k_star <-  h_big$lgamma[[p]][(k_star+1)]
        
        # if (k_star < h_small$K[p]) { 
        #   mu_k_star = 0.5*(h_small$lgamma[[p]][k_star] + h_small$lgamma[[p]][k_star+1])
        #   sigma = 0.1
        # }
        # if (k_star == h_small$K[p]) { 
        #   mu_k_sta r = h_small$lgamma[[p]][h_small$K[p]]
        #   sigma = 0.1
        # }
        
        if (h_small$K[p] > 1){
          mu_k_star = mean(h_small$lgamma[[p]][-1])
        } else {
          mu_k_star = 0
        }
        sigma = 0.1
        
        log_prob_proposal = dnorm(lgamma_k_star,mu_k_star,sigma,log = TRUE) + log(1/(smax[p]-h_small$K[p])) + log(pi_moves[h_small$K[p],1]) - 
          log(pi_moves[h_big$K[p],2]) - log(1/h_small$K[p]);
        log_prob_proposal = -log_prob_proposal
        h_small$gamma = lapply(1:M^2,function(p) { exp(h_small$lgamma[[p]])})
        h_small$s_diffs = lapply(1:M^2,function(p) { diff(h_small$s[[p]])})
        h_c = h_small;
      }

      # if (prod(h_c$alpha[[p]] == 0) == 1) {
      #   h_c$alpha[[p]] = c(0)
      #   h_c$lalpha[[p]] = c(-10000)
      #   h_c$delta[p] = 0;
      #   h_c$K[p] = 1
      #   h_c$s[[p]] = c(0,smax[p])
      # }

      log_lik_c = log_likelihood_multidim_dt(data,h_c,nu,op = 'vec')
      log_prior_lgamma_c = sum(dnorm(h_c$lgamma[[p]],mu_lgamma,s_lgamma, log=TRUE))
      log_prior_s_c = (factorial(h_c$K)*factorial(smax - h_c$K))/(factorial(smax))
      # assume all options for K are equally likely
      
      proba_accept = log_lik_c[m] - log_lik[m]  + log_prior_lgamma_c - log_prior_lgamma[p] + log_prior_s_c[p] - log_prior_s[p] - log_prob_proposal;
      test_accept = log(runif(1)) < proba_accept
      
      if ((test_accept == TRUE) & (move == "birth")) { 
        taux_accept_birth = taux_accept_birth + 1 
      }
      if ((test_accept == TRUE) & (move == "death")) { 
        taux_accept_death = taux_accept_death + 1 
      }  
      if (test_accept == TRUE) {  
        h = h_c
        log_lik = log_lik_c
        log_prior_lgamma[p]  = log_prior_lgamma_c
        log_prior_s  = log_prior_s_c
      }
      
    }

    
    #   ###########################################
    #----------------------- split and merge
    #########################################
    
    for (p in op_echan$K) {
      K_p = h$K[p];
      log_lik_c = log_lik
      m <- (p - 1) %/% M + 1
      h_c <- h;
      move = sample(c('split','merge'),1,prob = pi_moves[K_p,])
      if (move == 'split') {
        indic_simu = 'split'
        h_small <- h;
        h_big <- h;
        h_big$K[p] =  h_small$K[p]  +  1
        
        vacant_s = setdiff(seq(0,smax[p]), h_small$s[[p]])
        if (length(vacant_s)>0){
          
          if (length(vacant_s) ==1){
            k_star_val = as.numeric(vacant_s)
          } else {
            k_star_val = as.numeric(sample(vacant_s,1))
          }
          
          h_big$s[[p]] = sort(c(h_small$s[[p]], k_star_val))
          k_star = which(h_big$s[[p]]==k_star_val)-1
          
          v = runif(1,0,1)
          if (k_star != 1){
            lgamma_k_star1 = log(v * h_small$gamma[[p]][k_star] * (h_big$s[[p]][(k_star+2)] - h_big$s[[p]][(k_star)]) / 
                                   (h_big$s[[p]][(k_star+1)] - h_big$s[[p]][(k_star)]))
          } else {
            lgamma_k_star1 = log(1)
          }
          lgamma_k_star2 = log((1-v) * h_small$gamma[[p]][k_star] * (h_big$s[[p]][(k_star+2)] - h_big$s[[p]][(k_star)]) / 
                          (h_big$s[[p]][(k_star+2)] - h_big$s[[p]][(k_star+1)]))
            
          if (k_star == h_small$K[p]){ 
            if (h_small$K[p] == 1){
              h_big$lgamma[[p]] = c(lgamma_k_star1, lgamma_k_star2)
            } else {
              h_big$lgamma[[p]] = c(h_small$lgamma[[p]][1:(k_star-1)], lgamma_k_star1, lgamma_k_star2)
            }
          } else if (k_star < h_small$K[p]) { 
            if (k_star == 1){
              h_big$lgamma[[p]] = c(lgamma_k_star1, lgamma_k_star2, h_small$lgamma[[p]][(k_star+1):h_small$K[p]])
            } else{
              h_big$lgamma[[p]] = c(h_small$lgamma[[p]][1:(k_star-1)],lgamma_k_star1, lgamma_k_star2, h_small$lgamma[[p]][(k_star+1):h_small$K[p]])
            }
          } 
          stopifnot(length(h_big$lgamma[[p]]) == h_big$K[p])
          
          log_prob_proposal = log(1/(smax[p]-h_small$K[p])) + log(pi_moves[h_small$K[p],1]) - 
            log(pi_moves[h_big$K[p],2]) - log(1/h_small$K[p]) - 
            (h_small$gamma[[p]][k_star] * (h_big$s[[p]][(k_star+2)] - h_big$s[[p]][(k_star)])^2 / 
            ((h_big$s[[p]][(k_star+1)] - h_big$s[[p]][(k_star)])*(h_big$s[[p]][(k_star+2)] - h_big$s[[p]][(k_star+1)]))) # Jacobian
          h_big$s_diffs = lapply(1:M^2,function(p) { diff(h_big$s[[p]])})
          h_big$gamma = lapply(1:M^2,function(p) { exp(h_big$lgamma[[p]])})
          h_c = h_big;
        }
      }
      
      if (move == 'merge') {
        indic_simu = 'merge'
        h_big <- h
        h_small <- h ;
        h_small$K[p] <- h_big$K[p] - 1
        
        k_star <- as.numeric(sample(1:(h_small$K[p]),1))
        h_small$s[[p]] <- h_big$s[[p]][-(k_star + 1)]
        
        lgamma_k_star <- log((h_big$s[[p]][(k_star+1)] - h_big$s[[p]][(k_star)])*h_big$gamma[[p]][k_star] + 
                               (h_big$s[[p]][(k_star+2)] - h_big$s[[p]][(k_star+1)])*h_big$gamma[[p]][(k_star+1)])
          
        
        if (k_star == h_small$K[p]){ 
          if (h_small$K[p] == 1){
            h_small$lgamma[[p]] = c(h_big$lgamma[[p]][1])
          } else {
            h_small$lgamma[[p]] = c(h_big$lgamma[[p]][1:(k_star-1)], lgamma_k_star)
          }
        } else if (k_star < h_small$K[p]) { 
          if (k_star == 1){
            h_small$lgamma[[p]] = c(h_big$lgamma[[p]][1], h_big$lgamma[[p]][(k_star+2):h_big$K[p]])
          } else{
            h_small$lgamma[[p]] = c(h_big$lgamma[[p]][1:(k_star-1)],lgamma_k_star, h_big$lgamma[[p]][(k_star+2):h_big$K[p]])
          }
        } 
        stopifnot(length(h_small$lgamma[[p]]) == h_small$K[p])
        
        h_small$gamma = lapply(1:M^2,function(p) { exp(h_small$lgamma[[p]])})
        h_small$s_diffs = lapply(1:M^2,function(p) { diff(h_small$s[[p]])})
        
        log_prob_proposal = log(1/(smax[p]-h_small$K[p])) + log(pi_moves[h_small$K[p],1]) - 
          log(pi_moves[h_big$K[p],2]) - log(1/h_small$K[p]) - 
          (h_small$gamma[[p]][k_star] * (h_big$s[[p]][(k_star+2)] - h_big$s[[p]][(k_star)])^2 / 
             ((h_big$s[[p]][(k_star+1)] - h_big$s[[p]][(k_star)])*(h_big$s[[p]][(k_star+2)] - h_big$s[[p]][(k_star+1)]))) # Jacobian
        log_prob_proposal = -log_prob_proposal
        h_c = h_small;
      }
      
      # if (prod(h_c$alpha[[p]] == 0) == 1) {
      #   h_c$alpha[[p]] = c(0)
      #   h_c$lalpha[[p]] = c(-10000)
      #   h_c$delta[p] = 0;
      #   h_c$K[p] = 1
      #   h_c$s[[p]] = c(0,smax[p])
      # }
      
      log_lik_c = log_likelihood_multidim_dt(data,h_c,nu,op='vec')
      log_prior_lgamma_c = sum(dnorm(h_c$lgamma[[p]],mu_lgamma,s_lgamma, log=TRUE))
      log_prior_s_c = (factorial(h_c$K)*factorial(smax - h_c$K))/(factorial(smax))
      # assume all options for K are equally likely
      
      proba_accept = log_lik_c[m] - log_lik[m]  + log_prior_lgamma_c - log_prior_lgamma[p] + log_prior_s_c[p] - log_prior_s[p] - log_prob_proposal;
      test_accept = log(runif(1)) < proba_accept
      
      if ((test_accept == TRUE) & (move == "split")) { 
        taux_accept_split = taux_accept_split + 1 
      }
      if ((test_accept == TRUE) & (move == "merge")) { 
        taux_accept_merge = taux_accept_merge + 1 
      }  
      if (test_accept == TRUE) {  
        h = h_c
        log_lik = log_lik_c
        log_prior_lgamma[p]  = log_prior_lgamma_c
        log_prior_s  = log_prior_s_c
      }
      
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
  print(taux_accept_birth)
  print(taux_accept_death)
  print(taux_accept_split)
  print(taux_accept_merge)
  
  return(OUTPUT)  }






