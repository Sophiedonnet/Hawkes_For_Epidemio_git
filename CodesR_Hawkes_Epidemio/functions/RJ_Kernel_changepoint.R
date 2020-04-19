#library(mvtnorm)

RJ_Kernel_changepoint = function(data,INPUT,hyperParam_prior,par_algo_MCMC) {
  #op.plot=TRUE
  #op_in_SMC=FALSE
  ## output  : theta,log_lik,acc_rate
  ###
  M <- length(data$Times_obs)
  op_echan <- par_algo_MCMC$op_echan
  op_affichage <-  par_algo_MCMC$op_affichage
  N_MCMC <-  par_algo_MCMC$N_MCMC
  
  compt_accept_heights <-  0
  
  K_max <- 50
  pi_moves <- matrix(0.5,K_max,2);
  pi_moves[1,] <- c(1,0)  ### proba naissance / mort 
  
  #Times_utiles <- data$Times_utiles
  #Times_obs <- data$Times_obs
  
  tau = data$tau; if (length(tau) != M) { tau = rep(tau,M)}
  
  Tmax <- data$Tmax; if (length(Tmax) != M) {Tmax = rep(Tmax, M)}
  Tinf <- data$Tinf; if (length(Tinf) != M) {Tinf = rep(Tinf, M)}
  
  n_obs <- vapply(1:M, function(m) {length(Times_obs[[m]])}, 1)
  
  smax <- INPUT$theta$h$smax
  
  pi_0 = hyperParam_prior$pi_0
  p_Z = hyperParam_prior$p_Z # poids de 0 dans la loi a priori des h
  mu_lalpha  = hyperParam_prior$mu_lalpha #
  s_lalpha = hyperParam_prior$s_lalpha
  mu_lnu  = hyperParam_prior$mu_lnu### nu : partie constante de l'intensité_
  s_lnu  = hyperParam_prior$s_lnu
  
  ## Need to add hyperparameter prior on the log of the lag L and beta
  mu_llag  = hyperParam_prior$mu_llag 
  s_llag = hyperParam_prior$s_llag
  
  mu_lbeta  = hyperParam_prior$mu_lbeta 
  s_lbeta = hyperParam_prior$s_lbeta
  
  #a_lambda_K = hyperParam_prior$a_lambda_K
  #b_lambda_K = hyperParam_prior$b_lambda_K
  #a_s <- hyperParam_prior$a_s
  #pi_moves = prob_moves(K_min,K_max,p_K,grid_s=FALSE,op_echan_K=(length(op_echan$K)>0))
  
  
  lambda_K = INPUT$theta$lambda_K
  
  ############################## INITIALISATION ##########################
  nu = INPUT$theta$nu
  lnu = log(nu)
  h = INPUT$theta$h

  K = vapply(1:M^2,function(p) { length(h$alpha[[p]])},2)
  # delta = h$delta; 
  # s = h$s; 
  h$K = K
  
  
  ### add L and beta entries to INPUT
  Lag = INPUT$theta$Lag
  llag = log(Lag)
  
  beta = INPUT$theta$beta
  lbeta = log(beta)

  
  
  
  
  
  ##########################################################################
  
  
  
  rho_lnu = par_algo_MCMC$rho_lnu
  rho_lalpha = par_algo_MCMC$rho_lalpha
  
  ### Need to add sizes of move for L and beta
  rho_llag = par_algo_MCMC$rho_llag
  rho_lbeta = par_algo_MCMC$rho_lbeta
  
  if (op.plot == TRUE) {
    par(mfrow = c(M, M))
  }
  #####################
  
  seq_MCMC = list()
  seq_MCMC$nu = matrix(0, N_MCMC, M)
  seq_MCMC$nu[1, ] = nu
  seq_MCMC$h = vector(mode = "list", length = N_MCMC)
  seq_MCMC$h[[1]] = h
  seq_MCMC$lambda_K = rep(0, N_MCMC)
  seq_MCMC$lambda_K[1] <- lambda_K
  
  ### New entries
  seq_MCMC$Lag = matrix(0, N_MCMC, M)
  seq_MCMC$Lag[1, ] = Lag
  
  seq_MCMC$beta = matrix(0, N_MCMC, M)
  seq_MCMC$beta[1, ] = beta
  
  ### Need to duplicate all following computations
  Times_obs_1=lapply(1:M,function(m){u_m=Times[[m]]; u_m = u_m[u_m>Tinf]; u_m=u_m[u_m <= tau+Lag]})
  Times_obs_2=lapply(1:M,function(m){u_m=Times[[m]]; u_m = u_m[u_m>min(Tmax,tau+Lag)]; u_m=u_m[u_m <= Tmax]})
  
  Times_utiles_1=lapply(1:M,function(m){u_m=Times[[m]]; u_m=u_m[u_m > Tinf-max(smax)];  u_m=u_m[u_m <= tau+Lag]})
  mat_obs_utiles_1 =calc_mat_absc_utiles(Times_utiles_1,Times_obs_1,smax)
  
  Times_utiles_2=lapply(1:M,function(m){u_m=Times[[m]]; u_m=u_m[u_m > min(Tmax,tau+Lag)-max(smax)];  u_m=u_m[u_m <= Tmax]})
  mat_obs_utiles_2 =calc_mat_absc_utiles(Times_utiles_2,Times_obs_2,smax)
                                       
  idx_1 <- func_idx(mat_obs_utiles_1, h$s)
  J_1  <- calc_J(Times_utiles_1, h$s, Tinf, tau+Lag)
  idx_2 <- func_idx(mat_obs_utiles_2, h$s)
  J_2  <- calc_J(Times_utiles_2, h$s, tau+Lag, Tmax)
  
  log_lik = log_likelihood_hawkes_multidim_positive_intensity_2(data,  mat_obs_utiles_1,  mat_obs_utiles_2,
                                                                idx_1, idx_2, J_1, J_2, h, nu, beta, Lag, op ='vec')
  log_prior_h = dprior_h_moving_s(h, hyperParam_prior)
  log_prior_lnu = dnorm(lnu, mu_lnu, s_lnu, log = TRUE)
  
  ### Add prior for beta and L
  log_prior_llag = dnorm(llag, mu_llag, s_llag, log = TRUE)
  log_prior_lbeta = dnorm(lbeta, mu_lbeta, s_lbeta, log = TRUE)
  
  
  taux_accept_birth = taux_accept_death <- 0
  compt_accept_heights <-   0
  
  ########################################################
  ########################################################
  ########################################################
  for (it in 2:N_MCMC) {
    ## print(it)
    #if (it==2){print(c(it,nu,h$alpha[[1]]))}
    if (it %% par_algo_MCMC$op_affichage == 0) {
      print(c(it, nu, lambda_K))
    }
    
    
    ######################### nu | Times,aires
    
    
    for (m in op_echan$nu) {
      rho_lnu = c(0.01, 0.1, 0.5) * par_algo_MCMC$rho_lnu
      S = sample(1:3, 1)
      lnu_c = lnu
      
      lnu_c[m] = lnu[m] + rnorm(1, 0, rho_lnu[S])
      
      nu_c = exp(lnu_c)
      
      # new log likelihood
      log_lik_c = log_likelihood_hawkes_multidim_positive_intensity_2(data,  mat_obs_utiles_1,  mat_obs_utiles_2, idx_1, idx_2, J_1, J_2, h, nu_c, beta,Lag,
                                                                      op ='vec')
      #log_lik_c = vapply(1:M,function(m){length(Times_obs[[m]])*log(nu_c[m])-(Tmax[m]-Tinf[m])*nu_c[m]},1)
      
      log_prior_lnu_c_m = dnorm(lnu_c[m], mu_lnu[m], s_lnu[m], log = TRUE)
      p1 = log_lik_c[m] - log_lik[m]
      p2 = log_prior_lnu_c_m - log_prior_lnu[m]
      
      #print(c("nu",p1,p2))
      
      test_accept = log(runif(1)) <= (p1 + p2)
      
      
      if (test_accept == TRUE) {
        nu <- nu_c
        lnu <- lnu_c
        log_lik <- log_lik_c
        log_prior_lnu[m] <- log_prior_lnu_c_m
      }
    }
    
    ######## Adjusted Langevin Metropolis on log nu
    
    if (length(op_echan$nu) == M) {
      rho_lnu <- 0.2
      
      ## Separates t_i < tau + L and t_i > tau + L
      
      h_beta = duplicate(h)
      for (p in  1:M) {h_beta$alpha[[p]] =  h$alpha[[p]] * beta}
      
      
      eval_lambda_1 <- lambda_cond_multidim_mat_positive_intensity(mat_obs_utiles_1,idx_1,h,nu)
      eval_lambda_2 <-lambda_cond_multidim_mat_positive_intensity(mat_obs_utiles_2, idx_2, h_beta, nu)
      
      grad_log_post <-vapply(1:M, function(m) {nu[m] * (sum(1 / eval_lambda_1[[m]]) + sum(1 / eval_lambda_2[[m]]))}, 1) - 1 / s_lnu ^ 2 * (lnu - mu_lnu) - nu * (Tmax - Tinf)
      lnu_c <-lnu + rho_lnu * grad_log_post + sqrt(2 * rho_lnu) * rnorm(M)
      
      nu_c = exp(lnu_c)
      
      # New log likelihood
      log_lik_c = log_likelihood_hawkes_multidim_positive_intensity_2(data, mat_obs_utiles_1,  mat_obs_utiles_2, idx_1, idx_2, J_1, J_2, h, nu_c, beta,Lag,
                                                                      op ='vec')
      log_prior_lnu_c <- dnorm(lnu_c, mu_lnu, s_lnu, log = TRUE)
      q_nu_c_nu <- sum(dnorm(lnu_c, lnu + rho_lnu * grad_log_post, 2 * rho_lnu, log = TRUE))
      
      ## Same operation
      eval_lambda_c_1 <-lambda_cond_multidim_mat_positive_intensity(mat_obs_utiles_1, idx_1, h, nu_c)
      eval_lambda_c_2 <-lambda_cond_multidim_mat_positive_intensity(mat_obs_utiles_2, idx_2, h_beta, nu_c)
      grad_log_post_c <-vapply(1:M, function(m) {nu[m] * (sum(1 / eval_lambda_1[[m]]) + sum(1 / eval_lambda_2[[m]]))}, 1) 
                                    - 1 / s_lnu ^ 2 * (lnu_c - mu_lnu) - nu_c * (Tmax - Tinf)
      
      q_nu_nu_c <-
        sum(dnorm(lnu, lnu_c + rho_lnu * grad_log_post_c, 2 * rho_lnu, log = TRUE))
      proba_accept <-
        sum(log_lik_c - log_lik) + sum(log_prior_lnu_c - log_prior_lnu) - (q_nu_c_nu -
                                                                             q_nu_nu_c)
      
      #print(c("nu langevin",proba_accept))
      
      test_accept <- log(runif(1)) <= proba_accept
      
      if (test_accept == TRUE) {
        lnu <- lnu_c
        nu <- exp(lnu)
        log_lik <- log_lik_c
        log_prior_lnu <- log_prior_lnu_c
      }
    }
 
    
    #----------------------move heigths  ----------------
    
    for (p in op_echan$alpha) {
      #       for (p in 2:2){
      m <- (p - 1) %/% M + 1
      for (k in 1:(h$K[p])) {
        #  for (k in 1:1){
        
        
        h_c <- h
        
        rho_lalpha <- c(0.01, 0.1, 1) * par_algo_MCMC$rho_lalpha
        
        if (h$alpha[[p]][k] == 0) {
          h_c$lalpha[[p]][k] <- rnorm(1, mu_lalpha, s_lalpha)
          p_Z_move = 1
          
        } else{
          p_Z_move = 0.5#1-exp(-h$alpha[[p]][k]/10)
          Z = rbinom(1, 1, p_Z_move)
          
          if (Z == 0) {
            h_c$lalpha[[p]][k] = -10000
          } else{
            h_c$lalpha[[p]][k] = rnorm(1, h$lalpha[[p]][k], rho_lalpha[sample(1:3, 1)])
          }
        }
        h_c$alpha[[p]] = exp(h_c$lalpha[[p]])
        #------------ proba of move
        
        
        if ((sum(h_c$alpha[[p]] == 0) == h_c$K[p])) {
          proba_accept = -Inf
        } else{
          p_Z_move_c = 0.5#1-exp(-h_c$alpha[[p]][k]/10)
          d1 = sum(vapply(1:3, function(j) {
            dnorm(h_c$lalpha[[p]][k], h$lalpha[[p]][k], rho_lalpha[j])
          }, 1)) / 3
          d1c = sum(vapply(1:3, function(j) {
            dnorm(h$lalpha[[p]][k], h_c$lalpha[[p]][k], rho_lalpha[j])
          }, 1)) / 3
          p3_1 = log((h$alpha[[p]][k] == 0) * dnorm(h_c$lalpha[[p]][k], mu_lalpha, s_lalpha) + (h$alpha[[p]][k] >
                                                                                                  0) * ((1 - p_Z_move) * (h_c$alpha[[p]][k] == 0) + p_Z_move * (h_c$alpha[[p]][k] >
                                                                                                                                                                  0) * d1
                                                                                                  )
          )
          p3_2 = log((h_c$alpha[[p]][k] == 0) * dnorm(h$lalpha[[p]][k], mu_lalpha, s_lalpha) + (h_c$alpha[[p]][k] >
                                                                                                  0) * ((1 - p_Z_move_c) * (h$alpha[[p]][k] == 0) + p_Z_move_c * (h$alpha[[p]][k] >
                                                                                                                                                                    0) * d1c
                                                                                                  )
          )
          log_p_passage = p3_1 - p3_2
          #log_p_passage=0
          #       # ------------- mise ? jour
          #
          log_prior_h_c = dprior_h_moving_s(h_c, hyperParam_prior)
          
          ## New log likelihood
          log_lik_c = log_likelihood_hawkes_multidim_positive_intensity_2(data, mat_obs_utiles_1,  mat_obs_utiles_2, idx_1, idx_2, J_1, J_2, h_c, nu, beta,Lag,
                                                                          op ='vec')
 
          proba_accept = sum(log_lik_c - log_lik) - log_p_passage + log_prior_h_c$h -
            log_prior_h$h
        }
        
        #print(c("h",proba_accept))
        
        test_accept = log(runif(1)) < proba_accept
        
        if (test_accept == TRUE) {
          h = h_c
          
          log_lik = log_lik_c
          
          log_prior_h  = log_prior_h_c
          compt_accept_heights = compt_accept_heights + 1
        }
      }
    }
    #print(compt_accept_heights/truc*100)
    # #   #----------------------move steps ----------------
    # 
    # for (p in intersect(op_echan$s, which(h$K > 1))) {
    #   m <- (p - 1) %/% M + 1
    #   for (k in 2:h$K[p]) {
    #     #---- initialisation
    #     h_c = h
    #     
    #     s_c_p = s_p = h$s[[p]]
    #     #--- proposition
    #     rho_s = c(2, 4)
    #     S = sample(c(1, 2), 1)
    #     beta_s = (s_p[k + 1] - s_p[k - 1]) / (s_p[k] - s_p[k - 1]) * (rho_s[S] -
    #                                                                     1) - (rho_s[S] - 2)
    #     s_c_p[k] = s_p[k - 1] + (s_p[k + 1] - s_p[k - 1]) * rbeta(1, rho_s[S], beta_s)
    #     h_c$s[[p]] = s_c_p
    #     
    #     #-----  probability of move
    #     beta_1s = (s_p[k + 1] - s_p[k - 1]) / (s_p[k] - s_p[k - 1]) * (rho_s[1] -
    #                                                                      1) - (rho_s[1] - 2)
    #     beta_2s = (s_p[k + 1] - s_p[k - 1]) / (s_p[k] - s_p[k - 1]) * (rho_s[2] -
    #                                                                      1) - (rho_s[2] - 2)
    #     z = (s_c_p[k] - s_p[k - 1]) / (s_p[k + 1] - s_p[k - 1])
    #     
    #     p3_sc_s = log(0.5 / (s_p[k + 1] - s_p[k - 1]) * (dbeta(z, rho_s[1], beta_1s) +
    #                                                        dbeta(z, rho_s[2], beta_2s)))
    #     
    #     beta_1s_c = (s_c_p[k + 1] - s_c_p[k - 1]) / (s_c_p[k] - s_c_p[k -
    #                                                                     1]) * (rho_s[1] - 1) - (rho_s[1] - 2)
    #     beta_2s_c = (s_c_p[k + 1] - s_c_p[k - 1]) / (s_c_p[k] - s_c_p[k -
    #                                                                     1]) * (rho_s[2] - 1) - (rho_s[2] - 2)
    #     z_c = (s_p[k] - s_c_p[k - 1]) / (s_c_p[k + 1] - s_c_p[k - 1])
    #     p3_s_sc = log(0.5 / (s_c_p[k + 1] - s_c_p[k - 1]) * (
    #       dbeta(z_c, rho_s[1], beta_1s_c) + dbeta(z_c, rho_s[2], beta_2s_c)
    #     ))
    #     #--------- updates
    #     
    #     idx_c <- func_idx(mat_obs_utiles, h_c$s)
    #     J_c  <- calc_J(Times_utiles, h_c$s, Tinf, Tmax)
    #     
    #     log_prior_h_c = dprior_h_moving_s(h_c, hyperParam_prior)
    #     log_lik_c = log_lik
    #     log_lik_c[m] = log_likelihood_hawkes_multidim_positive_intensity_m(data, idx_c, J_c, h_c, nu, m)
    #     
    #     proba_accept = sum(log_lik_c - log_lik) + sum(log_prior_h_c$h -
    #                                                     log_prior_h$h) - (p3_sc_s - p3_s_sc)
    #     
    #     test_accept = log(runif(1)) < proba_accept
    #     
    #     if (test_accept == TRUE) {
    #       h = h_c
    #       J = J_c
    #       idx = idx_c
    #       log_lik = log_lik_c
    #       log_prior_h  = log_prior_h_c
    #     }
    #   }
    # }
    # #
    #
    #   #-----------------------------------------
    
    
    
    #   #----------------------birth and death  ----------------  ###### Do I need this part?
    #
    for (p in op_echan$K) {
      K_p = h$K[p]
      
      log_lik_c = log_lik
      m <- (p - 1) %/% M + 1
      h_c <- h
      
      move = sample(c('birth', 'death'), 1, prob = pi_moves[K_p, ])
      if (move == 'birth') {
        indic_simu = 'birth'
        h_small <- h
        
        h_big <- h
        
        h_big$K[p] =  h_small$K[p] + 1
        prob_k_star = diff(h_small$s[[p]]) / h_small$smax[p]
        k_star = sample(1:(h_small$K[p]), 1, prob = prob_k_star) #### on aura tendance à proposer de couper les longs paliers.
        h_big$s[[p]] = sort(c(h_small$s[[p]], runif(1, h_small$s[[p]][k_star], h_small$s[[p]][k_star +
                                                                                                1]))) #### au final s* est simulée uniformement sur [0, smax]
        if (k_star > 1) {
          mu_k_star = 0.5 * (h_small$lalpha[[p]][k_star - 1] + h_small$lalpha[[p]][k_star])
          sigma = 0.1
        }
        if (k_star == 1) {
          mu_k_star = h_small$lalpha[[p]][1]
          sigma = 0.1
        }
        lalpha_k_star = mu_k_star + sigma * rnorm(1)
        
        if (k_star == 1) {
          h_big$lalpha[[p]] = c(lalpha_k_star, h_small$lalpha[[p]])
        }
        if (k_star > 1) {
          h_big$lalpha[[p]] = c(h_small$lalpha[[p]][1:(k_star - 1)], lalpha_k_star, h_small$lalpha[[p]][(k_star):h_small$K[p]])
        }
        
        
        log_prob_proposal = dnorm(lalpha_k_star, mu_k_star, sigma, log =
                                    TRUE) + log(1 / h_small$smax[p]) + log(pi_moves[h_small$K[p], 1]) - log(pi_moves[h_big$K[p], 2])
        
        h_big$alpha = lapply(1:M ^ 2, function(p) {
          exp(h_big$lalpha[[p]])
        })
        h_c = h_big
        
      }
      
      if (move == 'death') {
        indic_simu = 'death'
        h_big <- h
        h_small <- h
        
        h_small$K[p] <- h_big$K[p] - 1
        k_star <- sample(1:(h_big$K[p] - 1), 1)
        h_small$s[[p]] <- h_big$s[[p]][-(k_star + 1)]
        
        h_small$lalpha[[p]] <- h_big$lalpha[[p]][-k_star]
        lalpha_k_star <-  h_big$lalpha[[p]][k_star]
        if (k_star > 1) {
          mu_k_star = 0.5 * (h_small$lalpha[[p]][k_star - 1] + h_small$lalpha[[p]][k_star])
          sigma = 0.1
        }
        # max(abs(h_small$lalpha[[p]][k_star-1]-h_small$lalpha[[p]][k_star]),0.1)}
        if (k_star == 1) {
          mu_k_star = h_small$lalpha[[p]][1]
          sigma = 0.1
        }#max(h_small$lalpha[[p]][k_star],0.1)}
        
        log_prob_proposal = dnorm(lalpha_k_star, mu_k_star, sigma, log =
                                    TRUE) + log(1 / h_small$smax[p]) + log(pi_moves[h_small$K[p], 1]) - log(pi_moves[h_big$K[p], 2])
        
        log_prob_proposal = -log_prob_proposal
        h_small$alpha = lapply(1:M ^ 2, function(p) {
          exp(h_small$lalpha[[p]])
        })
        h_c = h_small
        
      }
      
      if (prod(h_c$alpha[[p]] == 0) == 1) {
        h_c$alpha[[p]] = c(0)
        h_c$lalpha[[p]] = c(-10000)
        h_c$delta[p] = 0
        
        h_c$K[p] = 1
        h_c$s[[p]] = c(0, smax[p])
      }
      
      idx_c_1 <- func_idx(mat_obs_utiles_1, h_c$s)
      J_c_1  <- calc_J(Times_utiles_1, h_c$s, Tinf, tau+Lag)
      idx_c_2 <- func_idx(mat_obs_utiles_2, h_c$s) 
      J_c_2  <- calc_J(Times_utiles_2, h_c$s, tau+Lag, Tmax)
      log_lik_c = log_lik
      log_prior_h_c = dprior_h_moving_s(h_c, hyperParam_prior)
      

      ## New log likelihood
      log_lik_c = log_likelihood_hawkes_multidim_positive_intensity_2(data, mat_obs_utiles_1,  mat_obs_utiles_2,
                                                                      idx_c_1, idx_c_2, J_c_1, J_c_2, h_c, nu, beta,Lag,
                                                                      op ='vec')
      #print(c("h 2", log_lik_c[m], log_lik[m] , log_prior_h_c$h, log_prior_h$h, log_prob_proposal))

      proba_accept = log_lik_c[m] - log_lik[m] + log_prior_h_c$h - log_prior_h$h - log_prob_proposal
      
      test_accept = log(runif(1)) < proba_accept
      if ((test_accept == TRUE) &
          (move == "birth")) {
        taux_accept_birth = taux_accept_birth + 1
      }
      if ((test_accept == TRUE) &
          (move == "death")) {
        taux_accept_death = taux_accept_death + 1
      }  #
      if (test_accept == TRUE) {
        h = h_c
        J_1 = J_c_1
        idx_1 = idx_c_1
        J_2 = J_c_2
        idx_2 = idx_c_2
        
        log_lik = log_lik_c
        log_prior_h  = log_prior_h_c
      }
    }
    
    # 
    # if (op_echan$lambda_K == 1) {
    #   lambda_K = rgamma(
    #     1,
    #     hyperParam_prior$a_lambda_K + sum(h$K[h$delta == 1] - 1),
    #     hyperParam_prior$b_lambda_K + sum(h$delta == 1)
    #   )
    #   
    #   hyperParam_prior$lambda_K <- lambda_K
    # }
    #
    
    ###################################### Moves on beta ###########################################################
    
    for (m in op_echan$beta) {
      rho_lbeta = c(0.01, 0.1, 0.5) * par_algo_MCMC$rho_lbeta
      S = sample(1:3, 1)
      lbeta_c = lbeta
      
      lbeta_c[m] = lbeta[m] + rnorm(1, 0, rho_lbeta[S])
      
      beta_c = exp(lbeta_c)
      
      # only second part of log likelihood
      log_lik_c = log_likelihood_hawkes_multidim_positive_intensity_2(data, mat_obs_utiles_1,  mat_obs_utiles_2,idx_1,
                                                                    idx_2, J_1, J_2, h, nu, beta_c,Lag,
                                                                      op ='vec')

      log_prior_lbeta_c_m = dnorm(lbeta_c[m], mu_lbeta[m], s_lbeta[m], log = TRUE)
      p1 = log_lik_c[m] - log_lik[m]
      p2 = log_prior_lbeta_c_m - log_prior_lbeta[m]
      
      test_accept = log(runif(1)) <= (p1 + p2)
      
      #print(c("beta",p1,p2))
      
      
      if (test_accept == TRUE) {
        beta = beta_c
        lbeta = lbeta_c
        log_lik = log_lik_c
        log_prior_lbeta[m] = log_prior_lbeta_c_m
      }
    }
    
    ###################################### Moves on lag L ###########################################################
    
    ### Moves the lag only 
    
    for (m in op_echan$lag) {
      rho_llag = c(0.01, 0.1, 0.5) * par_algo_MCMC$rho_llag
      S = sample(1:3, 1)
      llag_c = llag
      
      llag_c[m] = llag[m] + rnorm(1, 0, rho_llag[S])
      
      lag_c = exp(llag_c)
      
      Times_obs_c_1=lapply(1:M,function(m){u_m=Times[[m]]; u_m = u_m[u_m>Tinf]; u_m=u_m[u_m <= tau+lag_c]})
      Times_obs_c_2=lapply(1:M,function(m){u_m=Times[[m]]; u_m = u_m[u_m>min(Tmax,tau+lag_c)]; u_m=u_m[u_m <= Tmax]})
      
      Times_utiles_c_1=lapply(1:M,function(m){u_m=Times[[m]]; u_m=u_m[u_m > Tinf-max(smax)];  u_m=u_m[u_m <= tau+lag_c]})
      mat_obs_utiles_c_1 =calc_mat_absc_utiles(Times_utiles_x_1,Times_obs_c_1,smax)
      
      Times_utiles_c_2=lapply(1:M,function(m){u_m=Times[[m]]; u_m=u_m[u_m > min(Tmax,tau+lag_c)-max(smax)];  u_m=u_m[u_m <= Tmax]})
      mat_obs_utiles_c_2 =calc_mat_absc_utiles(Times_utiles_c_2,Times_obs_c_2,smax)
      
      idx_c_1 <- func_idx(mat_obs_utiles_c_1, h$s)
      J_c_1  <- calc_J(Times_utiles_c_1, h$s, Tinf, tau+lag_c)
      idx_c_2 <- func_idx(mat_obs_utiles_c_2, h$s)
      J_c_2  <- calc_J(Times_utiles_c_2, h$s, tau+lag_c, Tmax)
      
      
      
      ## New log likelihood
      log_lik_c = log_likelihood_hawkes_multidim_positive_intensity_2(data, mat_obs_utiles_c_1,  mat_obs_utiles_c_2,
                                                                      idx_c_1, idx_c_2, J_c_1, J_c_2, h, nu, beta,lag_c,
                                                                      op ='vec')
      
      log_prior_llag_c_m = dnorm(llag_c[m], mu_llag[m], s_llag[m], log = TRUE)
      p1 = log_lik_c[m] - log_lik[m]
      p2 = log_prior_llag_c_m - log_prior_llag[m]
      
      #print(c("lag",p1,p2))
      
      test_accept = log(runif(1)) <= (p1 + p2)
      
      
      if (test_accept == TRUE) {
        Lag = lag_c
        llag = llag_c
        log_lik = log_lik_c
        log_prior_llag[m] = log_prior_llag_c_m
        
        Times_obs_1 = Times_obs_c_1
        Times_obs_2 = Times_obs_c_2
        mat_obs_utiles_1 =  mat_obs_utiles_c_1
        mat_obs_utiles_2 = mat_obs_utiles_c_2
        idx_1 = idx_c_1
        idx_2 =  idx_c_2
        J_1 = J_c_1
        J_2 = J_c_2
      }
    }
    
    
    
    #   #===============   STOCKAGE des variables
    #
    
    
    seq_MCMC$nu[it, ] = nu
    seq_MCMC$h[[it]] = h
    seq_MCMC$lambda_K[it] =  lambda_K
    seq_MCMC$Lag[it,] =  Lag
    seq_MCMC$beta[it,] =  beta
    
    
    
  }
  OUTPUT = seq_MCMC
  return(OUTPUT)
}
