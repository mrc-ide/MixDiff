RunMCMC <- function(obs_dat, 
                    MCMC_settings,
                    hyperpriors,
                    index_dates,
                    index_dates_order) ### CHANGE THIS SO index_dates_order is computed automatically from index_dates
{
  
  ###############################################
  ### define augmented data to be used for initialisation of the chain ###
  ###############################################
  
  aug_dat <- initialise_aug_data(obs_dat, index_dates_order)
  
  ###############################################
  ### define parameters to be used for initialisation of the chain ###
  ###############################################
  
  index_dates <- list(matrix(c(1, 2), nrow=2), cbind(c(1, 2), c(1, 3)), cbind(c(1, 2), c(2, 3), c(1, 4)), cbind(c(1, 2), c(2, 3), c(1, 4)) )
  theta <- initialise_theta_from_aug_dat(aug_dat, index_dates)
  
  ###############################################
  ### Initalise the MCMC chains ###
  ###############################################
  
  range_dates <- find_range(obs_dat)
  
  # to store param values
  curr_theta <- theta
  theta_chain <- list()
  theta_chain[[1]] <- curr_theta
  
  # to store augmented data values
  curr_aug_dat <- aug_dat
  aug_dat_chain <- list()
  aug_dat_chain[[1]] <- curr_aug_dat
  
  logpost_chain <- rep(NA, (MCMC_settings$chain_properties$n_iter - MCMC_settings$chain_properties$burnin) / MCMC_settings$chain_properties$record_every)
  logpost_chain[1] <- lposterior_total(curr_aug_dat, curr_theta, obs_dat, hyperpriors, index_dates, range_dates)
  
  n_accepted_D_moves <- 0
  n_proposed_D_moves <- 0
  
  n_accepted_mu_moves <- 0
  n_proposed_mu_moves <- 0
  
  n_accepted_CV_moves <- 0
  n_proposed_CV_moves <- 0
  
  ###############################################
  ### Run the MCMC ###
  ###############################################
  
  print("... Burnin ...")
  
    for(k in seq_len(MCMC_settings$chain_properties$n_iter-1))
    {
      output_stuff <- (k>=MCMC_settings$chain_properties$burnin) & (k %% MCMC_settings$chain_properties$record_every)==0
      
      if(output_stuff)
      {
        print(sprintf("... %d / %d ...", k, MCMC_settings$chain_properties$n_iter))
      }
      
      # move some of the D_i
      if(MCMC_settings$moves_switch$D_on)
      {
        for(g in seq_len(n_groups))
        {
          for(j in seq_len(ncol(curr_aug_dat$D[[g]])))
          {
            to_update <- sample(1:nrow(obs_dat[[g]]), round(nrow(obs_dat[[g]])*MCMC_settings$moves_options$fraction_Di_to_update)) # proposing moves for only a certain fraction of dates
            n_10_to_update <- floor(length(to_update) / MCMC_settings$moves_options$move_D_by_groups_of_size)
            for(i in seq_len(length(n_10_to_update)))
            {
              tmp <- move_Di (to_update[MCMC_settings$moves_options$move_D_by_groups_of_size*(i-1)+(seq_len(MCMC_settings$moves_options$move_D_by_groups_of_size))], g, j, 
                              curr_aug_dat,
                              curr_theta, 
                              obs_dat, 
                              hyperpriors, 
                              index_dates,
                              range_dates) 
              n_proposed_D_moves <- n_proposed_D_moves + 1
              n_accepted_D_moves <- n_accepted_D_moves + tmp$accept
              if(tmp$accept==1) curr_aug_dat <- tmp$new_aug_dat # if accepted move, update accordingly
            }
          }
        }
      }
      
      # move zeta using Gibbs sampler
      if(MCMC_settings$moves_switch$zeta_on)
      {
        tmp <- move_zeta_gibbs(curr_aug_dat,
                               curr_theta, 
                               obs_dat, 
                               hyperpriors) 
        curr_theta <- tmp$new_theta # always update with new theta (Gibbs sampler)
      }
      
     # move mu
      if(MCMC_settings$moves_switch$mu_on)
      {
        for(g in seq_len(n_groups))
        {
          for(j in seq(2,ncol(curr_aug_dat$D[[g]]),1))
          {
            tmp <- move_lognormal(what="mu", g, j-1, MCMC_settings$moves_options$sdlog_mu, 
                                  curr_aug_dat,
                                  curr_theta, 
                                  obs_dat, 
                                  hyperpriors,
                                  index_dates)
            n_proposed_mu_moves <- n_proposed_mu_moves + 1
            n_accepted_mu_moves <- n_accepted_mu_moves + tmp$accept
            if(tmp$accept==1) curr_theta <- tmp$new_theta # if accepted move, update accordingly
          }
        }
      }
      
      # move CV
      if(MCMC_settings$moves_switch$CV_on)
      {
        for(g in seq_len(n_groups))
        {
          for(j in seq(2,ncol(curr_aug_dat$D[[g]]),1))
          {
            tmp <- move_lognormal(what="CV", g, j-1, MCMC_settings$moves_options$sdlog_CV, 
                                  curr_aug_dat,
                                  curr_theta, 
                                  obs_dat, 
                                  hyperpriors,
                                  index_dates)
            n_proposed_CV_moves <- n_proposed_CV_moves + 1
            n_accepted_CV_moves <- n_accepted_CV_moves + tmp$accept
            if(tmp$accept==1) curr_theta <- tmp$new_theta # if accepted move, update accordingly
          }
        }
      }
      
      # recording value of parameters and corresponding posterior after all moves 
      if( output_stuff )
      {
        idx <- (k - MCMC_settings$chain_properties$burnin) / MCMC_settings$chain_properties$record_every+1
        theta_chain[[idx]] <- curr_theta
        aug_dat_chain[[idx]] <- curr_aug_dat
        logpost_chain[idx] <- lposterior_total(curr_aug_dat, curr_theta, obs_dat, hyperpriors, index_dates, range_dates) #### CONSIDER DOING THIS USING SAPPLY AFTER THE WHOLE THING
      }
    }
  
  ###############################################
  ### Compute acceptance probabilities ###
  ###############################################
  
  accept_prob <- list(
    D_moves=n_accepted_D_moves / n_proposed_D_moves,
    mu_moves=n_accepted_mu_moves / n_proposed_mu_moves,
    CV_moves=n_accepted_CV_moves / n_proposed_CV_moves,
    zeta_moves=1)
  
  ###############################################
  ### Return lisst of outputs of interest ###
  ###############################################
  
  res <- list(theta_chain=theta_chain, aug_dat_chain=aug_dat_chain, logpost_chain=logpost_chain, accept_prob=accept_prob)
  
  return(res)
}